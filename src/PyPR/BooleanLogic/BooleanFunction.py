import json
from collections.abc import Callable, Collection, Iterator, Mapping
from typing import Any, Protocol, Self

from numba import njit  # noqa: F401 -- used by the source compile() execs

import PyPR.JSON_Serialization

from PyPR.BooleanLogic.Latex import LatexStyle, LatexTerm, fill_template

# stack markers for the DAG walks: `postorder` uses _FINISHED to mean "the node
# beneath me is finished"; the tree printers use _CLOSE and _SEPARATOR for the
# closing bracket and the comma between arguments
_FINISHED = object()
_CLOSE = object()
_SEPARATOR = object()


class _RenderedValues(dict):
    """Rendered values keyed by node, falling back to override strings on a miss.

    Generation hooks look their arguments up with `values[arg]`. Keeping this a
    real dict keeps that lookup at C speed; only an override node -- which is
    never rendered, so never stored -- reaches `__missing__`.
    """

    overrides: Mapping[Any, str]

    def __missing__(self, node: Any) -> str:
        return self.overrides[node]


class _CopiesWithMemo(dict):
    """Copies keyed by original node, falling back to the deepcopy memo on a miss.

    `__deepcopy__` stops its walk at nodes copied earlier in the same deepcopy
    call; their copies live in the memo (keyed by id), so a lookup of one of them
    misses here and is answered from there.
    """

    memo: dict[int, Any]

    def __missing__(self, node: Any) -> Any:
        return self.memo[id(node)]


class IndexableContainer[K,V](Protocol):
    def __getitem__(self, key: K, /) -> V: ...

class BooleanFunction:
    args: tuple["BooleanFunction", ...]
    arg_limit: int | None

    def __init__(self,
        *args: "BooleanFunction",
        arg_limit: int | None = None
    ):
        self.args = args
        self.arg_limit = arg_limit

    def postorder(self,
        stop: Callable[["BooleanFunction"], bool] | None = None
    ) -> Iterator["BooleanFunction"]:
        """Iterate over the nodes of the DAG, each once, every node after all of its arguments.

        This is the walk every whole-function method is built on. The method keeps
        its results in a dict keyed by node and fills it in the order yielded, so
        by the time a node arrives the entries for its arguments are already there::

            values = {}
            for node in fn.postorder():
                values[node] = node._eval(values, array)

        What differs between methods -- what to compute at a leaf versus a gate,
        where to store it, where to stop -- stays in the loop body and the `stop`
        predicate, so none of it has to be expressed as a hook on the walk.

        Nodes come in the order a left-to-right depth-first search finishes them,
        and a node reached along several paths is yielded the first time only.
        Nodes are identified by object identity, so two equal but distinct `VAR(0)`
        objects are two nodes. The walk is iterative, so deep functions do not run
        into Python's recursion limit.

        :param stop: A predicate marking nodes to treat as already handled: a node
            for which it returns True is neither yielded nor descended into. It
            is used to resume from an earlier call (nodes that already have an id,
            a SAT label or a copy) and to cut the DAG at boundaries (subfunctions
            and overrides in code generation). The caller supplies the values of
            stopped nodes itself. Defaults to None, which walks the whole DAG.
        :type stop: Callable[[BooleanFunction], bool] | None
        :yield: Each node of the DAG not cut off by `stop`, arguments first.
        :rtype: Iterator[BooleanFunction]
        """
        # This loop runs under every eval, so it is written for speed: it
        # allocates nothing per node and has no hook calls when `stop` is None.
        # Measured against the per-method loops it replaced, it is no slower.
        seen: set[BooleanFunction] = set()
        stack: list[Any] = [self]
        pop = stack.pop

        while stack:
            node = pop()

            # a node entered earlier sits beneath this marker, and every
            # argument pushed above the marker has now been yielded
            if node is _FINISHED:
                yield pop()
                continue

            if node in seen:
                continue
            if stop is not None and stop(node):
                continue
            seen.add(node)

            args = node.args
            if args:
                # arguments are pushed reversed so they are popped left to right
                stack += (node, _FINISHED)
                stack += args[::-1]
            else:
                # nothing beneath it, so a leaf is finished as soon as it is reached
                yield node

    def _copy(self,
        child_copies: dict["BooleanFunction","BooleanFunction"]
    ) -> Self:
        """Rebuild this node over copies of its children.

        The tree walks behind `__copy__`, `compose` and `_merge_redundant` copy
        a function bottom-up: by the time a node is reached, each of its
        children has a copy in `child_copies`, and this builds the node's own
        copy on top of them. It is an instance method because a node is more
        than its children -- a gate carries `arg_limit`, a fusion node its
        template -- and that additional data can only come from the node being
        copied. A subclass with data of its own overrides this to carry it
        over, copying it so that the result shares nothing with the original.

        :param child_copies: Copies of this node's children, keyed by original.
        :type child_copies: dict[BooleanFunction, BooleanFunction]
        :return: A copy of this node over the copied children.
        :rtype: Self
        """
        return type(self)(
            *(child_copies[arg] for arg in self.args),
            arg_limit = self.arg_limit
        )

    def copy(self) -> Self:
        """An alias of `__copy__()` which creates a copy of a BooleanFunction.

        Creates a copy of the boolean function on which it is called.
        the output DAG structure is identical to the DAG structure
        of the function on which it is called, and the returned function
        is the same subclass as the input.       
        
        :return copy: a copy of the input function
        """
        return self.__copy__()

    def __copy__(self) -> Self:
        """Creates a copy of a BooleanFunction.

        Creates a copy of the boolean function on which it is called.
        the output DAG structure is identical to the DAG structure
        of the function on which it is called, and the returned function
        is the same subclass as the input. The copy is always deep, so
        `copy.copy`, `copy.deepcopy` and `.copy()` all give the same result.

        :return copy: a copy of the input fn.
        """
        return self.__deepcopy__({})

    def __deepcopy__(self,
        memo: dict[int, Any]
    ) -> Self:
        """Creates a copy of a BooleanFunction which shares no nodes with the original.

        The DAG is rebuilt bottom-up, so a node reached along several paths is
        copied once and the copy has the same sharing as the original. `memo` is
        the `copy.deepcopy` memo, keyed by `id`: a node already in it is reused
        rather than copied again. That carries the sharing across everything
        copied in one `deepcopy` call -- for example the bits of a feedback
        function which reference common subfunctions.

        :param memo: The `copy.deepcopy` memo, mapping `id` of an original to its copy.
        :type memo: dict[int, Any]
        :return: A copy of the input function.
        :rtype: Self
        """
        copies = _CopiesWithMemo()
        copies.memo = memo

        # Nodes copied earlier in the same deepcopy call are cut off the walk, and
        # their copies come from the memo when a parent looks them up. A top-level
        # copy starts from an empty memo, so it skips the per-node check entirely.
        stop = (lambda n: id(n) in memo) if memo else None

        for node in self.postorder(stop=stop):
            if node.is_leaf():
                # Overwritten in Inputs.py
                new_node = node.__copy__()
            else:
                new_node = node._copy(copies)

            copies[node] = new_node
            memo[id(node)] = new_node

        return memo[id(self)]

    def add_arguments(
        self,
        *new_args: "BooleanFunction"
    ) -> None:
        """Adds one or more arguments to a BooleanFunction.

        :param tuple[BooleanFunction] new_args: A variable length list of arguments to add   
        
        :raises ValueError: if the number of arguments would make len(self.args)
            greater than the allowed number of arguments (self.arg_limit).
        """
        if (not self.arg_limit) or (len(self.args) + len(new_args) <= self.arg_limit):
            self.args = tuple(list(self.args) + list(new_args))
        else:
            raise ValueError(f"{type(self)} object supports at most {self.arg_limit} arguments")

    def remove_arguments(
            self,
            *remove_args: "BooleanFunction"
        ) -> None:
        """Removes one or more arguments from a BooleanFunction.

        Removes the specified arguments from the function. If one or more of the arguments are 
        are not present, they are skipped (no error is raised). Additionally, if no arguments are
        passed, all arguments are removed.

        :param tuple[BooleanFunction] remove_args: A variable length list of arguments to remove  
        """
        if not remove_args:
            self.args = ()
        else:
            self.args = tuple(x for x in self.args if x not in remove_args)

    def subfunctions(self) -> list["BooleanFunction"]:
        """Returns a list of subfunctions

        A subfunction is a BooleanFunction object which is referenced 
        more than one time by BooleanFunctions above it in the DAG. The
        subfunctions are returned in the order they are first encountered
        by a DFS (postorder) traveral. This means the returned array is
        topologically sorted.

        :returns subfunctions: A topologically sorted list of BooleanFunction 
            which are referenced by multiple parents
        """
        # a node is a subfunction when more than one argument slot points at it.
        # Each node is yielded once, so scanning every node's arguments visits
        # every incoming edge exactly once -- including both edges when a gate
        # takes the same argument twice.
        nodes = []
        referenced: set[BooleanFunction] = set()
        shared: set[BooleanFunction] = set()
        for node in self.postorder():
            nodes.append(node)
            for arg in node.args:
                if arg in referenced:
                    shared.add(arg)
                else:
                    referenced.add(arg)

        return [node for node in nodes if node in shared and not node.is_leaf()]

    def inputs(self) -> list["BooleanFunction"]:
        """Returns a list of the input nodes

        The list returns all leaf BooleanFunction objects (e.g. `VAR` or `CONST`), 
        in the order they are first encountered by a DFS (postorder) traveral.
        The function does not merge or exclude any objects (even if they are semantically 
        equivalent). In other words, two VAR objects might both be included if they 
        are distinct objects, even if they reference the same variable.

        :returns inputs: A list of all the leaf nodes of the DAG, which serve as 
            inputs to the function
        """
        return [node for node in self.postorder() if node.is_leaf()]


    # string generation:
    def pretty_str(self) -> str:
        """Returns a nicely formatted string for pretty printing.

        This returns a string with functions opened and closed on new lines. Subfunctions
        are broken out and printed in a topologically sorted order, so that the DAG structure 
        can be seen. Each function is printed with opening and closing parenthesis, inside which
        Arguments are printed on new lines and indented. this takes up much more vertical space,
        but is easier to parse.

        :returns pretty_string: a nicely formatted string for pretty printing
        """
        subfuncs = self.subfunctions()
        labels = {root: f"(subfunction {i+1})" for i, root in enumerate(subfuncs)}
        pretty_strings = []

        # the indent for each depth, built once; depths are reached one level at
        # a time, so a depth not yet cached is always the next one
        prefixes: list[str] = []

        # a tree walk rather than a fold: a gate's line opens before its
        # arguments and its bracket closes after them, so every line is emitted
        # exactly once, in order, and the cost is linear in the output
        for root in subfuncs + [self]:
            name = "(Main Function)" if root is self else labels[root]
            lines = [f"{name} = (\n"]

            stack: list[tuple[Any, int]] = [(root, 0)]
            while stack:
                node, depth = stack.pop()
                if depth < len(prefixes):
                    prefix = prefixes[depth]
                else:
                    prefix = "   |" * depth + "   "
                    prefixes.append(prefix)

                if node is _CLOSE:
                    lines.append(prefix + ")\n")
                elif node is not root and node in labels:
                    lines.append(prefix + labels[node] + "\n")
                elif node.is_leaf():
                    lines.append(prefix + node.dense_str() + "\n")
                else:
                    lines.append(prefix + f"{type(node).__name__} (\n")
                    stack.append((_CLOSE, depth))
                    stack.extend((arg, depth + 1) for arg in reversed(node.args))

            lines.append(")\n")
            pretty_strings.append("".join(lines))

        return "\n\n".join(pretty_strings)

    def dense_str(self) -> str:
        """Returns a densely formatted string for debugging and variable inspection.

        Subfunctions are abbreviated, but not printed. Functions are printed inline,
        inside nested parenthesis. Although these choices make it harder to parse, they 
        make it safer to print large functions, and in almost all cases they still reveal 
        enough of the structure to debug or explore further. Due to the increased safety, this is
        the default behavior for `__str__()`, and `pretty_str` must be called explicitly.
         
         
        :returns dense_string: a densely formatted string for debugging and variable inspection.
        """
        subfuncs = self.subfunctions()
        labels = {root: f"(subfunction {i+1})" for i, root in enumerate(subfuncs)}

        # a tree walk rather than a fold: a gate opens before its arguments and
        # closes after them, so the pieces are appended in output order and
        # joined once, which keeps the cost linear in the output
        pieces: list[str] = []
        stack: list[Any] = [self]
        while stack:
            node = stack.pop()

            if node is _CLOSE:
                pieces.append(")")
            elif node is _SEPARATOR:
                pieces.append(",")
            elif node in labels:
                pieces.append(labels[node])
            elif node.is_leaf():
                pieces.append(node.dense_str())
            else:
                pieces.append(f"{type(node).__name__}(")
                stack.append(_CLOSE)

                # pushed last-first, with a separator between each pair, so they
                # pop as arg0 , arg1 , ... argN
                args = node.args
                for i in range(len(args) - 1, -1, -1):
                    stack.append(args[i])
                    if i:
                        stack.append(_SEPARATOR)

        return "".join(pieces)

    def __str__(self):
        """An alias for `dense_str`, which returns a densely \
        formatted string for debugging and variable inspection.

        Subfunctions are abbreviated, but not printed. Functions are printed inline,
        inside nested parenthesis. Although these choices make it harder to parse, they 
        make it safer to print large functions, and in almost all cases they still reveal 
        enough of the structure to debug or explore further. Due to the increased safety, this is
        the default behavior for `__str__()`, and `pretty_str` must be called explicitly.
         
        :returns dense_string: a densely formatted string for debugging and variable inspection.
        """
        return self.dense_str()

    # code generation:
    def _generate_assignments(self,
        hook: str,
        context: Any,
        output_name: str,
        subfunction_name: Callable[[int], str],
        overrides: Mapping["BooleanFunction", str] | None = None,
        inline_subfunctions: bool = False,
    ) -> list[tuple[str, Any]]:
        """Split the function into named assignments and render each one bottom-up.

        This is the driver shared by every `generate_*` method; they differ only
        in how one node is rendered and how an assignment is laid out. Each
        subfunction gets its own assignment, in topological order, followed by
        one for the function itself. An assignment's expression is rendered over
        its *region*: the nodes reachable from it without passing through another
        subfunction or an override. Within a region, an assigned subfunction
        stands for its name and an override node for its override string, so
        the hook sees a plain string for both.

        :param hook: The name of the per-node hook, called on every rendered node
            as `node.<hook>(values, context)`, e.g. `"_generate_c"`.
        :type hook: str
        :param context: The hook's second argument: the state array's name for
            code, the `LatexStyle` for LaTeX.
        :type context: Any
        :param output_name: The name assigned the function itself.
        :type output_name: str
        :param subfunction_name: The name of the subfunction with the given
            1-based index.
        :type subfunction_name: Callable[[int], str]
        :param overrides: Nodes to write as the given string instead of rendering
            them. An override node gets no assignment of its own.
        :type overrides: Mapping[BooleanFunction, str] | None
        :param inline_subfunctions: If True, no subfunctions are split out, and the
            function is rendered as a single assignment.
        :type inline_subfunctions: bool
        :return: One (name, rendered value) pair per assignment, dependencies first.
        :rtype: list[tuple[str, Any]]
        """
        if overrides is None:
            overrides = {}

        roots = [] if inline_subfunctions else self.subfunctions()
        names = {root: subfunction_name(i + 1) for i, root in enumerate(roots)}
        roots.append(self)
        names[self] = output_name

        # Values already known are not re-rendered, so the walk stops at them:
        # overrides, assigned subfunctions (by name), and leaves shared with an
        # earlier region. Override nodes are never stored, so looking one up
        # misses and falls through to its string. Without overrides the stop
        # test is the dict's own C-level __contains__, with no Python call per
        # node; with a single region and no overrides nothing is known in
        # advance, so there is nothing to stop at and no test at all.
        # (no __init__ on these dict subclasses: attributes are set after
        # construction, which keeps construction a C-level call)
        values = _RenderedValues()
        values.overrides = overrides
        if overrides:
            def stop(node: BooleanFunction) -> bool:
                return node in values or node in overrides
        elif len(roots) > 1:
            stop = values.__contains__
        else:
            stop = None

        assignments = []
        for root in roots:
            # an override replaces its node everywhere, including its assignment
            if root in overrides:
                continue

            for node in root.postorder(stop=stop):
                values[node] = getattr(node, hook)(values, context)

            assignments.append((names[root], values[root]))
            values[root] = names[root]

        return assignments

    def _generate_c(self,
        c_strings: Mapping["BooleanFunction", str],
        array_name: str
    ) -> str:
        """Hook which specifies how to generate C for a BooleanFunction subclass.

        Given the C expressions of the node's arguments and the name of the state
        array, return the C expression for this node. `generate_c` calls this on
        every node, so a custom subclass generates valid C as long as it
        implements it.

        :param c_strings: The C expression of each argument, keyed by node.
        :type c_strings: Mapping[BooleanFunction, str]
        :param array_name: The name of the state array, which `VAR` nodes index.
        :type array_name: str
        :raises NotImplementedError: If not overridden.
        :return: A C expression computing this node.
        :rtype: str
        """
        raise NotImplementedError  # implemented for each subclass

    def generate_c(self,
        output_name: str = 'output',
        subfunction_prefix: str = 'fn',
        array_name: str = 'array',
        overrides: Mapping["BooleanFunction", str] | None = None
    ) -> list[str]:
        """Generates C statements computing the function.

        Each subfunction is assigned to its own variable, in topological order,
        and the last statement assigns the function itself to `output_name`. The
        statements are bare assignments: declaring the variables and the state
        array is left to the caller.

        :param output_name: The variable the function is assigned to. Defaults to 'output'.
        :type output_name: str
        :param subfunction_prefix: Subfunction variables are named
            `{subfunction_prefix}_{index}`. Defaults to 'fn'.
        :type subfunction_prefix: str
        :param array_name: The name of the array holding the state. Defaults to 'array'.
        :type array_name: str
        :param overrides: Nodes to write as the given string instead of generating
            them; an overridden node gets no statement of its own. Used to refer to
            expressions already computed elsewhere. Defaults to None.
        :type overrides: Mapping[BooleanFunction, str] | None
        :return: One C statement per line, dependencies first.
        :rtype: list[str]
        """
        assignments = self._generate_assignments(
            "_generate_c",
            array_name,
            output_name,
            lambda index: f"{subfunction_prefix}_{index}",
            overrides,
        )
        return [f"{name} = {expr};" for name, expr in assignments]

    def _generate_VHDL(self,
        vhdl_strings: Mapping["BooleanFunction", str],
        array_name: str
    ) -> str:
        """Hook which specifies how to generate VHDL for a BooleanFunction subclass.

        Given the VHDL expressions of the node's arguments and the name of the state
        signal, return the VHDL expression for this node. `generate_VHDL` calls this
        on every node, so a custom subclass generates valid VHDL as long as it
        implements it.

        :param vhdl_strings: The VHDL expression of each argument, keyed by node.
        :type vhdl_strings: Mapping[BooleanFunction, str]
        :param array_name: The name of the state signal, which `VAR` nodes index.
        :type array_name: str
        :raises NotImplementedError: If not overridden.
        :return: A VHDL expression computing this node.
        :rtype: str
        """
        raise NotImplementedError  # implemented for each subclass

    def generate_VHDL(self,
        output_name: str = 'output',
        subfunction_prefix: str = 'fn',
        array_name: str = 'array',
        overrides: Mapping["BooleanFunction", str] | None = None
    ) -> list[str]:
        """Generates VHDL signal assignments computing the function.

        Each subfunction is assigned to its own signal, in topological order, and
        the last assignment drives `output_name` with the function itself. This is
        mostly used by `FeedbackFunction.write_VHDL`, which declares the signals.

        :param output_name: The signal the function drives. Defaults to 'output'.
        :type output_name: str
        :param subfunction_prefix: Subfunction signals are named
            `{subfunction_prefix}_{index}`. Defaults to 'fn'.
        :type subfunction_prefix: str
        :param array_name: The name of the state signal. Defaults to 'array'.
        :type array_name: str
        :param overrides: Nodes to write as the given string instead of generating
            them; an overridden node gets no assignment of its own. Defaults to None.
        :type overrides: Mapping[BooleanFunction, str] | None
        :return: One VHDL assignment per line, dependencies first.
        :rtype: list[str]
        """
        assignments = self._generate_assignments(
            "_generate_VHDL",
            array_name,
            output_name,
            lambda index: f"{subfunction_prefix}_{index}",
            overrides,
        )
        return [f"{name} <= {expr};" for name, expr in assignments]

    def _generate_python(self,
        python_strings: Mapping["BooleanFunction", str],
        array_name: str
    ) -> str:
        """Hook which specifies how to generate Python for a BooleanFunction subclass.

        Given the Python expressions of the node's arguments and the name of the
        state array, return the Python expression for this node. `generate_python`
        calls this on every node, so a custom subclass generates valid Python as
        long as it implements it.

        :param python_strings: The Python expression of each argument, keyed by node.
        :type python_strings: Mapping[BooleanFunction, str]
        :param array_name: The name of the state array, which `VAR` nodes index.
        :type array_name: str
        :raises NotImplementedError: If not overridden.
        :return: A Python expression computing this node.
        :rtype: str
        """
        raise NotImplementedError  # implemented for each subclass

    def generate_python(self,
        output_name: str = 'output',
        subfunction_prefix: str = 'fn',
        array_name: str = 'array',
        overrides: Mapping["BooleanFunction", str] | None = None
    ) -> list[str]:
        """Generates Python statements computing the function.

        Each subfunction is assigned to its own variable, in topological order, and
        the last statement assigns the function itself to `output_name`. The
        expressions use only `^`, `&`, `|` and `1 - x`, so they are valid numba
        code over integer bits; `compile` and `FeedbackFunction.compile` are built
        on this.

        :param output_name: The variable the function is assigned to. Defaults to 'output'.
        :type output_name: str
        :param subfunction_prefix: Subfunction variables are named
            `{subfunction_prefix}_{index}`. Defaults to 'fn'.
        :type subfunction_prefix: str
        :param array_name: The name of the array holding the state. Defaults to 'array'.
        :type array_name: str
        :param overrides: Nodes to write as the given string instead of generating
            them; an overridden node gets no statement of its own. Defaults to None.
        :type overrides: Mapping[BooleanFunction, str] | None
        :return: One Python statement per line, dependencies first.
        :rtype: list[str]
        """
        assignments = self._generate_assignments(
            "_generate_python",
            array_name,
            output_name,
            lambda index: f"{subfunction_prefix}_{index}",
            overrides,
        )
        return [f"{name} = {expr}" for name, expr in assignments]

    def _generate_latex(self,
        terms: Mapping["BooleanFunction", LatexTerm | str],
        style: LatexStyle
    ) -> LatexTerm:
        """Hook which specifies how to render a BooleanFunction subclass in LaTeX.

        The default looks the node's class name up in `style.operators`, so a new
        gate class only needs an entry there. Override this for a node that is not
        an operator joining its arguments (as the leaves and fused nodes do).

        :param terms: The rendering of each argument, keyed by node. A plain string
            is a subfunction name or override, and is treated as an atom.
        :type terms: Mapping[BooleanFunction, LatexTerm | str]
        :param style: The notation to render in.
        :type style: LatexStyle
        :return: The rendered node.
        :rtype: LatexTerm
        """
        return style.render_gate(type(self).__name__, [terms[arg] for arg in self.args])

    def generate_latex(self,
        output_name: str = 'f',
        subfunction_name: str | Callable[..., str] = 'g_{$index}',
        style: LatexStyle | None = None,
        overrides: Mapping["BooleanFunction", str] | None = None,
        inline_subfunctions: bool = False
    ) -> list[str]:
        """Generates LaTeX equations for the function.

        Each subfunction gets its own equation, in topological order, followed by
        one defining the function itself, so a DAG with shared structure reads as
        a short system of equations rather than one expression with repeats. Set
        `inline_subfunctions` to write the whole function as one expression.

        Parentheses are placed by precedence rather than around every gate, so an
        ANF reads `x_{0} x_{1} \\oplus x_{2}`. The notation itself -- operator
        symbols, precedence, negation, variables, constants, delimiters and the
        layout of each line -- comes from `style`; see `LatexStyle`. The default
        line layout is `name &= expr`, for an `align` environment::

            lines = fn.generate_latex()
            tex = "\\\\begin{align*}\\n" + " \\\\\\\\\\n".join(lines) + "\\n\\\\end{align*}"

        :param output_name: The left-hand side of the final equation. Defaults to 'f'.
        :type output_name: str
        :param subfunction_name: Template naming a subfunction, with field `$index`
            (1-based), or a callable taking `index`. Defaults to 'g_{$index}'.
        :type subfunction_name: str | Callable[..., str]
        :param style: The notation to write in. Defaults to `LatexStyle()`.
        :type style: LatexStyle | None
        :param overrides: Nodes to write as the given LaTeX instead of expanding
            them; an overridden node is treated as an atom and gets no equation of
            its own. Defaults to None.
        :type overrides: Mapping[BooleanFunction, str] | None
        :param inline_subfunctions: If True, write the function as one equation
            with no subfunctions split out. Defaults to False.
        :type inline_subfunctions: bool
        :return: One equation per line, dependencies first.
        :rtype: list[str]
        """
        if style is None:
            style = LatexStyle()

        assignments = self._generate_assignments(
            "_generate_latex",
            style,
            output_name,
            lambda index: fill_template(subfunction_name, index=index),
            overrides,
            inline_subfunctions,
        )
        return [
            style.format_line(name, term if isinstance(term, str) else term.text)
            for name, term in assignments
        ]


    # convenient node manipulations
    def _binarize(self,
        cache: dict["BooleanFunction", "BooleanFunction"]
        ) -> Self:
        """Helper function which specifies how to binarize a BooleanFunction subclass

        Creates an equivalent version of the function in which all gates have at most 2 inputs.

        :param cache: A dict which maps child nodes to their binarized output.
        :type cache: dict[BooleanFunction, BooleanFunction]
        :raises NotImplementedError: If not overriden
        :return: A new BooleanFunction which is the binarized version of the input node.
        :rtype: BooleanFunction
        """
        raise NotImplementedError

    def binarize(self) -> Self:
        """Creates an equivalent version of the function in which all gates have at most 2 inputs.

        For the associative gates (`XOR`, `AND`, and `OR`) any gate which has more than two
        args will be split into a tree of binary gates of the same type. The negated versions
        (`XNOR`, `NAND`, and `NOR`) are translated into a tree of the corresponding associative
        gate with a negated gate at the root. `NOT` gates are unaffected. This is useful for 
        analyses which can only handle binary gates (such as the tseytin transform).

        :return: A new BooleanFunction which is the binarized version of the input function.
        :rtype: BooleanFunction
        """
        new_nodes = {}
        for node in self.postorder():
            if node.is_leaf():
                # Overwritten in Inputs.py
                new_nodes[node] = node.__copy__()
            else:
                new_nodes[node] = node._binarize(new_nodes)
        return new_nodes[self]

    def _remap_indices(self,
        index_map: Any
    ) -> None:
        """Helper function which remaps the index of a particular leaf node.

        Modifies the node in-place to remap the index of a particular leaf node (usually `VAR`).
        Does not affect `CONST` nodes.

        :param index_map: Any container which supports indexing via integers 
            (e.g. `index_map[old_index] = new_index`). This is very frequently 
            a dict or list, but other types are accepted as well.
        :type index_map: Any
        :raises NotImplementedError: If not overriden
        """
        raise NotImplementedError

    def remap_indices(self,
        index_map: IndexableContainer[int,int],
        in_place: bool = False
    ) -> Self:
        """Remap the input variable indices.

        :param index_map: Any container which supports indexing via integers 
            (e.g. `index_map[old_index] = new_index`). This is very frequently 
            a dict or list, but other types are accepted as well.
        :type index_map: Any
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to `False`
        :type in_place: bool, optional
        :return: A BooleanFunction with the indices remapped.
        :rtype: BooleanFunction
        """
        if in_place: fn = self
        else: fn = self.copy()

        for leaf in fn.inputs():
            leaf._remap_indices(index_map)
        return fn

    def _remap_constants(self,
        const_map: list[tuple[Any,Any]]
    ) -> None:
        """Helper function which remaps the constant of a particular leaf node.

        Modifies the node in-place to remap the constant of a particular leaf node 
        (usually `CONST`). Does not affect `VAR` nodes.

        :param const_map: A container which contains (old_const, new_const) pairs. Every instance of
            old_const will be replaced with the first instance of old_const found (so to prevent 
            errors ensure a stable iteration order on the container and sane equality checks on
            the constants). This pair structure allows for keys/constants which are not hashable. 
        :type const_map: Any
        :raises NotImplementedError: If not overriden
        """
        raise NotImplementedError

    def remap_constants(self,
        constant_map: list[tuple[Any,Any]],
        in_place: bool = False
    ) -> Self:
        """Remap the input constants.

        :param const_map: A container which contains (old_const, new_const) pairs. Every instance of
            old_const will be replaced with the first instance of old_const found (so to prevent 
            errors ensure a stable iteration order on the container and sane equality checks on
            the constants). This pair structure allows for keys/constants which are not hashable. 
        :type const_map: Any
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to `False`
        :type in_place: bool, optional
        :return: A BooleanFunction with the constants remapped.
        :rtype: BooleanFunction
        """
        if in_place: fn = self
        else: fn = self.copy()

        for leaf in fn.inputs():
            leaf._remap_constants(constant_map)
        return fn

    def _shift_indices(self,
        shift_amount: int
    ) -> None:
        """Helper function which remaps the index of a particular leaf node.

        Modifies the node in-place to remap the index of a particular leaf node (usually `VAR`).
        Does not affect `CONST` nodes. This is similar to a simplified version of`_remap_indices` 
        where every index is mapped from `i` to `i+shift_amount` rather than according to a lookup

        :param shift_amount: The amount to shift each index (e.g. `i` becomes `i + shift amount`)
        :type index_map: int
        :raises NotImplementedError: If not overriden
        """
        raise NotImplementedError

    def shift_indices(self,
        shift_amount: int,
        in_place: bool = False
    ) -> Self:
        """Shift the input variable indices by a fixed amount.

        :param shift_amount: The amount to shift each index (e.g. `i` becomes `i + shift amount`)
        :type index_map: int
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to `False`
        :type in_place: bool, optional
        :return: A BooleanFunction with the indices shifted.
        :rtype: BooleanFunction
        """
        if in_place: fn = self
        else: fn = self.copy()

        for leaf in fn.inputs():
            leaf._shift_indices(shift_amount)
        return fn

    def condense_idxs(self,
        in_place: bool = False
    ) -> Self:
        """Reduce all variable indices to be consecutive 

        Reduces all variable indices as much as possible while maintaining the relative
        order. This also makes all variables consecutive, which reduces the size of the
        list needed to pass inputs to a function.
        
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to`False`
        :type in_place: bool, optional
        :return: A BooleanFunction with consecutive variables
        :rtype: BooleanFunction
        """
        return self.remap_indices(
            {v:i for i,v in enumerate(sorted(self.idxs_used()))},
            in_place
        )

    def _compose(self,
        input_map: IndexableContainer[int,"BooleanFunction"],
        in_place: bool = False
    ) -> Self:
        """Helper function which helps compose functions.

        Returns a function where each input node has been replaced with the function
        at the corresponding index in input_map. The inplace flag can be set to do this
        modification in place, rather than constructing a new function.
        No simplifications are made to the newly composed function.

        :param input_map: Any container which supports indexing via integers 
            (e.g. `index_map[old_index] = function`). This is very frequently 
            a dict or list, but other types are accepted as well.
        :type input_map: Any
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to `False`
        :type in_place: bool, optional
        :raises NotImplementedError: If not overriden
        """
        raise NotImplementedError

    def compose(self,
        input_map: IndexableContainer[int,"BooleanFunction"],
        in_place: bool = False
    ) -> Self:
        """Compose a BooleanFunction with a container mapping input variables to other BooleanFunctions

        Use the given functions to replace leaf nodes which have their inndices in the container.
        Not all indices have to be mapped (and indices which aren't will remain unaffection)

        :param input_map: Any container which supports indexing via integers 
            (e.g. `index_map[old_index] = function`). This is very frequently 
            a dict or list, but other types are accepted as well.
        :type input_map: Any
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to `False`
        :type in_place: bool, optional
        :return: A BooleanFunction with the indices remapped.
        :rtype: BooleanFunction
        """
        new_nodes = {}
        for node in self.postorder():
            if node.is_leaf():
                # Overwritten in Inputs.py
                new_nodes[node] = node._compose(input_map, in_place)
            elif in_place:
                node.args = tuple(new_nodes[arg] for arg in node.args)
                new_nodes[node] = node
            else:
                new_nodes[node] = node._copy(new_nodes)
        return new_nodes[self]

    def _merge_redundant(self,
        cache: dict["BooleanFunction","BooleanFunction"],
        subfunctions: Collection["BooleanFunction"],
        in_place: bool = False
    ) -> "BooleanFunction":
        """Helper function which determines how to simplify a node for `merge_redundant`

        determines how to reduce or simplify functions, and is expected to be overridden'
        for nodes which can be simplified better than the default. For example, for the
        associative gates (`XOR`,`AND`, and `OR`), this method pulls out arguments of the
        same type, and raises them up, effectively merging nodes as much as possible. Other
        simplifications might be possible for your given node.

        :param cache: a dictionary mapping child nodes to their corresponding output.
        :type cache: dict[BooleanFunction,BooleanFunction]
        :param subfunctions: The subfunctions of the root node on which
            `merge_redundant` was called. Only membership is meaningful: it is
            passed as a set, so testing it is constant time.
        :type subfunctions: Collection[BooleanFunction]
        :param in_place: If `True`, modify the function in place and return self,
        instead of returning a new function. Defaults to `False`
        :type in_place: bool, optional
        :return: A reduced or simplified version of this node.
        :rtype: BooleanFunction
        """
        if len(self.args) == 1:
            return cache[self.args[0]]
        elif in_place:
            self.args = tuple([cache[arg] for arg in self.args])
            return self
        else:
            return self._copy(cache)

    def merge_redundant(self,
        in_place: bool = False
    ) -> "BooleanFunction":
        """Performs some basic heuristic simplifications on a BooleanFunction

        Performs several basic heuristic simplifications. By default, first removes 
        non-unary functions with only one input, as these have no effect 
        on the output. Secondly, merge any associative gates (e.g. XOR, AND, OR),
        unless the child is a subfunction (does not include leaves).
        Even these basic simplifications can result in a dramatically
        simplified function when applied recursively, and help clean up 
        structures created during function composition or other modifications.
        custom simplifications can be added by overriding the `_merge_redundant()`
        helper function.

        :param in_place: If set to false (by default), the method will return a
            new function, leaving the orignial unmodified (highly recommended usage). 
            If set to true, it will modify the function in place as much as is possible.
            However, because some simplifications change the root node (and thus 
            cannot be done in place), these changes are ommitted. The fully simplified
            function will still be returned, but it is possible for the two to be different
            (i.e. it is possible that `fn = fn.merge_redundant(in_place=True)` is more
            simplified than, and thus not the same as, `fn.merge_redundant(in_place=True)`).
        :type in_place: bool, optional
        :return: A reduced and simplified version of the input function.
        :rtype: BooleanFunction
        """
        # a set, not a list: the hooks only test membership, once per argument,
        # and a list made that quadratic in the number of subfunctions
        subfunctions = set(self.subfunctions())

        new_nodes = {}
        for node in self.postorder():
            if node.is_leaf():
                new_nodes[node] = node
            else:
                new_nodes[node] = node._merge_redundant(
                    new_nodes, subfunctions, in_place = in_place,
                )
        return new_nodes[self]

    # evaluation
    def _eval(self,
        values: dict["BooleanFunction", Any],
        array: IndexableContainer[int, Any]
    ) -> Any:
        """Helper function which determines how to evaluate a single node.

        :param values: A dictionary which maps child nodes to their evaluations
        :type values: dict[BooleanFunction, Any]

        :param array: A container which can be indexed (e.g. array[idx]) to get the value
            of the variable at that index. Often a list (hence the name array), but other 
            types which support indexing (such as dict or numpy.ndarray) also work. 
        :type array: IndexableContainer[int, Any]
        :raises NotImplementedError: If not overriden
        :return: The value of the functions evaluation at this node.
        :rtype: Any
        """
        raise NotImplementedError

    def eval(self,
        array: IndexableContainer[int, Any]
    ) -> Any:
        """Evaluate the function on a given set of inputs.

        :param array: A container which can be indexed (e.g. array[idx]) to get the value
            of the variable at that index. Often a list (hence the name array), but other 
            types which support indexing (such as dict or numpy.ndarray) also work. 
        :type array: IndexableContainer[int, Any]
        :raises NotImplementedError: If not overriden
        :return: The value of the functions evaluation at this node.
        :rtype: Any
        """
        values = {}
        for node in self.postorder():
            values[node] = node._eval(values, array)
        return values[self]

    def _eval_ANF(self,
        values: dict["BooleanFunction", Any],
        array: IndexableContainer[int, Any]
    ) -> Any:
        """Helper function which determines how to evaluate a single node
        using only AND, XOR, and negation.

        :param values: A dictionary which maps child nodes to their evaluations
        :type values: dict[BooleanFunction, Any]

        :param array: A container which can be indexed (e.g. array[idx]) to get the value
            of the variable at that index. Often a list (hence the name array), but other 
            types which support indexing (such as dict or numpy.ndarray) also work. 
        :type array: IndexableContainer[int, Any]
        :raises NotImplementedError: If not overriden
        :return: The value of the functions evaluation at this node.
        :rtype: Any
        """
        raise NotImplementedError

    def eval_ANF(self,
        array: IndexableContainer[int,Any]
    ):
        """Evaluate the function on a given set of inputs using only AND, XOR, and negation.

        :param array: A container which can be indexed (e.g. array[idx]) to get the value
            of the variable at that index. Often a list (hence the name array), but other 
            types which support indexing (such as dict or numpy.ndarray) also work. 
        :type array: IndexableContainer[int, Any]
        :raises NotImplementedError: If not overriden
        :return: The value of the functions evaluation at this node.
        :rtype: Any
        """
        values = {}
        for node in self.postorder():
            values[node] = node._eval_ANF(values, array)
        return values[self]

    def compile(self) -> Any:
        """Just-In-Time compiles a given function, enabling faster evaluation.

        Generates the function as a python function and compiles it with numba's njit
        functionality. The resulting function is then both stored inside the function's
        `_compiled` field, and also returned, allowing for immediate use. Any changes to
        the function will not be reflected until the function is recompiled. Additionally,
        due to the limitations of the numba jit compiler, this only works for evaluations
        of integers in ndarrays, and not the more varied objects supported by `eval` and 
        `eval_ANF`

        :return: The compiled function. The return type is actually a numba CPUDispatcher,
            but it is annotated as `Any` to avoid causing unwarranted type errors in downstream
            applications.
        :rtype: Any
        """
        self._compiled = None
        python_body = "\n    ".join(self.generate_python())

        exec(f"""
@njit(parallel=True)
def _compiled(array):
    {python_body}
    return output
self._compiled = _compiled
""")
        return self._compiled

    # Methods from BooleanANF.py:
    @classmethod
    def from_ANF(cls,
        anf: Any
    ) -> "BooleanFunction":
        """Generates a BooleanFunction from either a BooleanANF or a nested iterable.

        If the passed anf is a BooleanANF, convert it to a BooleanFunction (equivalent
        to `anf.to_BooleanFunction()`). Otherwise, anf is assumed to be a nested iterable,
        which can be converted to BooleanANF via the BooleanANF constructor. Then the
        returned function is equivalent to `BooleanANF(anf).to_BooleanFunction`. This
        enables forming a BooleanFunction from a nested list, or data structure which
        can represent ANF via literals in a convenient way.

        :param anf: The ANF object to be used to create the function. Because there are
            a variety of types available to be parsed by the BooleanANF constructor, the
            type is annotated as `Any` to avoid downstream type errors. 
        :type anf: Any
        :return: The boolean function corresponding to the input ANF.
        :rtype: BooleanFunction
        """
        raise NotImplementedError  # defined in ANF.py

    def translate_ANF(self) -> "BooleanFunction":
        """Convert a BooleanFunction to its ANF representation

        This function translates a BooleanFunction to another BooleanFunction, but in
        ANF form. This means that it is top level XOR gate, the arguments of which are
        either CONST gates or an AND of VAR nodes. Note that this is the same as 
        `BooleanANF.from_BooleanFunction(fn).to_BooleanFunction`, but that the output
        not be confused with BooleanFunction - This is a "round-trip" back to BooleanFunctions
        
        :return: A BooleanFunction which models an ANF equation.
        :rtype: BooleanFunction
        """
        raise NotImplementedError  # defined in ANF.py

    def anf_str(self) -> str:
        """Print a BooleanFunction using the ANF string format 
        
        This is identical to `str(BooleanANF.from_BooleanFunction(fn))`, and is mostly
        included to make the code less verbose. A minor benefit is that it makes it slightly
        easier to reproduce some of the really old experiments

        :return: The string for the BooleanANF of the function
        :rtype: str
        """
        raise NotImplementedError # defined in ANF.py

    def degree(self) -> int:
        """Calculate the algebraic degree of the function

        This is implemented by computing the ANF, and then taking the max of the
        degree of each monomial in the ANF. This is equivalent to
        `BooleanANF.from_BooleanFunction(fn).degree()`. As with all methods which 
        depend on the ANF, this has the potential to become computationally infeasible 
        (but is usually fine).

        :return: The algebraic degree of the function
        :rtype: int
        """
        raise NotImplementedError # defined in ANF.py

    def monomial_count(self) -> int:
        """Calculate the number of monomials in the ANF of the function

        This is implemented by computing the ANF, and then counting the number of the
        monomials in the ANF. This is equivalent to 
        `len(BooleanANF.from_BooleanFunction(fn).degree().terms)`, but is more readible
        and far less verbose. As with all methods which depend on the ANF, this has the 
        potential to become computationally infeasible (but is usually fine).

        :return: The number of monomials in the ANF of the function
        :rtype: int
        """
        raise NotImplementedError # defined in ANF.py

    # def anf_optimize(self, translate=True):
    #     raise NotImplementedError # defined in ANF.py

    # Methods from SAT.py
    def _tseytin_labels(self,
        node_labels: dict["BooleanFunction", list[int]],
        variable_labels: dict[int, int],
        next_idx: int
    ) -> int:
        """Claim the solver variables this node needs for its own wires.

        The default is one wire per intermediate gate of a binary chain, which
        is what every gate in `Gates.py` expands to. A node that encodes itself
        some other way -- a direct n-ary CNF, or a fused subtree -- overrides
        this and claims whatever it needs. The only requirement is that the
        last label claimed is the wire carrying this node's result, since that
        is what a parent wires itself to.

        :param node_labels: Node-to-wires map, extended in place with this node.
        :type node_labels: dict[BooleanFunction, list[int]]
        :param variable_labels: Input-variable-to-wire map; only leaves add to it.
        :type variable_labels: dict[int, int]
        :param next_idx: The next unused solver variable.
        :type next_idx: int
        :return: The next unused solver variable after this node's claim.
        :rtype: int
        """
        num_labels = max(1, len(self.args) - 1)
        node_labels[self] = [next_idx + i for i in range(num_labels)]
        return next_idx + num_labels

    def _tseytin_clauses(self,
        label_map: dict["BooleanFunction", list[int]]
    ) -> list[tuple[int, ...]]:
        """Return the clauses relating this node's wires to its arguments'.

        This is the sole polymorphic entry point for clause generation: the
        walk reaches a node only through here, once, after every argument has
        been labelled. Everything about *how* a class encodes itself is private
        to that class, and this signature is the whole contract -- its own
        wires in, its arguments' result wires in, clauses out.

        The map holds every node labelled so far, so a node takes its own wires
        from `label_map[self]` and an argument's result from `label_map[arg][-1]`.
        Only that last entry is a shared convention -- the layout of the rest of
        a node's list is its own business, which is what lets a gate chain binary
        formulas, a direct n-ary encoding claim one wire, and a fused subtree
        label its whole interior.

        :param label_map: Wires for every node encoded so far, including this
            one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: A list of clauses encoding the wire relationship for the node
        :rtype: list[tuple[int, ...]]
        :raises NotImplementedError: If not implemented for the node class
        """
        raise NotImplementedError

    def tseytin(self,
        prev_clauses: list[tuple[int, ...]] | None = None,
        prev_node_labels: dict['BooleanFunction',list[int]] | None = None,
        prev_variable_labels: dict[int,int] | None = None
    ) -> tuple[
        list[tuple[int, ...]],
        dict["BooleanFunction",list[int]],
        dict[int,int]
    ]:
        """Generate constraints corresponding to the tseytin transformation of the function.

        This function primarily deals with three objects:
        - **Clauses**: A list of constraints, suitable for a sat solver. Each constraint is a tuple
            of integers, which represents an OR clause which must be satisfied (CNF). Negative integers
            represent the negation of a variable.
        - **Node Labels**: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used 
            to convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        - **Variable Labels**: Because the input variables are always supposed to be the same (even
            in different VAR `nodes`), we use this dict to map input variables to the corresponding
            variables in the sat solver. This can be used to convert the sat solution back to a satisfying
            assignment of input variables, or to impose additional conditions based on extra information.

        When called with no inputs, this function will create these three from scratch. However,
        the function can also recieve the output of previous calls as input, and will continue to build
        out the constraints, reusing variables where possible. This means the same function can also be
        used to add onto existing output. For example:
        ```python
        tseytin_data = function_1.tseytin()
        tseytin_data = function_2.tseytin(*tseytin_data)
        tseytin_data = function_3.tseytin(*tseytin_data)
        ...
        ```
        Or equivalently called in a for loop will create a set of constraints which describes all functions
        and which reuses variables where appropriate. This updating lets you create complicated queries,
        if you need to interact with the sat solver more deeply than just calling `fn.sat()`

        
        :param prev_clauses: A list of constraints, suitable for a sat solver. Each constraint is a tuple
            of integers, which represents an OR clause which must be satisfied (CNF). Negative integers
            represent the negation of a variable.
        :type prev_clauses: list[tuple[int, ...]]
        :param prev_node_labels: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used to
            convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        :type prev_node_labels: dict[BooleanFunction,list[int]]
        :param prev_variable_labels: Because the input variables are always supposed to be the same (even
            in different VAR `nodes`), we use this dict to map input variables to the corresponding
            variables in the sat solver. This can be used to convert the sat solution back to a satisfying
            assignment of input variables, or to impose additional conditions based on extra information.
        :type prev_variable_labels: dict[int,int]
        :return: (Clauses, Node Labels, Variable Labels)
        :rtype: tuple[
            list[tuple[int, ...]],
            dict[BooleanFunction,list[int]],
            dict[int,int]
        ]
        """
        raise NotImplementedError  # defined in SAT.py

    def tseytin_labels(self,
        node_labels: dict["BooleanFunction", list[int]] |None = None,
        variable_labels: dict[int,int] | None = None
    ) -> tuple[
        dict["BooleanFunction", list[int]],
        dict[int,int]
    ]:
        """Generate the labels for the tseytin transformation of the function.

        Tseytin primarily deals with three objects: Clauses, Node Labels and variable labels. This
        function is responsible for generating the latter two, which are then used to create clauses.
        In more detail, this  method produces:
        - **Node Labels**: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used 
            to convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        - **Variable Labels**: Because the input variables are always supposed to be the same (even
            in different VAR `nodes`), we use this dict to map input variables to the corresponding
            variables in the sat solver. This can be used to convert the sat solution back to a satisfying
            assignment of input variables, or to impose additional conditions based on extra information.

        When called with no inputs, this function will create these two from scratch. However,
        the function can also recieve the output of previous calls as input, and will continue to build
        out the constraints, reusing variables where possible. This means the same function can also be
        used to add onto existing output. For example:
        ```python
        labels = function_1.tseytin_labels()
        labels = function_2.tseytin_labels(*labels)
        labels = function_3.tseytin_labels(*labels)
        ...
        ```
        Or equivalently called in a for loop will create a set of labels which covers all functions
        and which reuses variables where appropriate. This is an important part of the `tseytin` function.
        
        :param prev_node_labels: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used to
            convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        :type prev_node_labels: dict[BooleanFunction,list[int]]
        :param prev_variable_labels: Because the input variables are always supposed to be the same (even
            in different VAR `nodes`), we use this dict to map input variables to the corresponding
            variables in the sat solver. This can be used to convert the sat solution back to a satisfying
            assignment of input variables, or to impose additional conditions based on extra information.
        :type prev_variable_labels: dict[int,int]
        :return: (Node Labels, Variable Labels)
        :rtype: tuple[
            dict[BooleanFunction,list[int]],
            dict[int,int]
        ]
        """
        raise NotImplementedError  # defined in SAT.py

    def tseytin_clauses(self,
        label_map: dict["BooleanFunction", list[int]]
    ) -> list[tuple[int, ...]]:
        """Generate constraints corresponding to the tseytin transformation of the function.

        `Tseytin` primarily deals with three objects: Clauses, Node Labels and variable labels. This
        function is responsible for generating the clauses, and requires the node labels. The Clauses 
        returned are list of constraints, suitable for a sat solver. Each constraint is a tuple
        of integers, which represents an OR clause which must be satisfied (CNF). Negative integers
        represent the negation of a variable.

        :param node_labels: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used to
            convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        :type node_labels: dict[BooleanFunction,list[int]]
        :return: clauses which encode the given function for a sat solver.
        :rtype: list[tuple[int, ...]],
        """
        raise NotImplementedError  # defined in SAT.py

    def sat(self,
        solver_name: str = "cadical195",
        verbose: bool = False,
    ) -> dict[int,bool] | None:
        """Solve the SAT problem for a given BooleanFunction

        First, the `tseytin` method is used to build a SAT encoding of the function.
        Then, the output is manually asserted to be true, and the resulting clauses are
        used as input for a SAT solver provided by the PySAT package. If more complicated
        sat-based procedures are needed, `tseytin` is exposed, and more complicated instances
        can be created manually. However, because direct SAT solving is the most common
        application, this function provides a much simpler interface, abstracting away
        details of the encoding. 

        :param solver_name: A string giving the name of a sat solver provided by PySAT
            which will be used as the solver, defaults to "cadical195", which we observed
            to work well experimentally.
        :type solver_name: str, optional
        :param verbose: if True, print statistics and timings for debugging, defaults to False
        :type verbose: bool, optional

        :return: if the function is unsatisfiable, return None. otherwise, returns a dictionary
            which maps variables to their boolean values in a satisfying assignment. Any variables 
            which don't appear in the dict are "don't care" .
        :rtype: dict[int,bool] | None
        """
        raise NotImplementedError  # defined in SAT.py

    def enum_models(self,
        solver_name: str = 'cadical195',
        verbose: bool = False
    ) -> Iterator[dict[int,bool]]:
        """Enumerate solutions to the SAT problem for a given BooleanFunction

        First, the `tseytin` method is used to build a SAT encoding of the function.
        Then, the output is manually asserted to be true, and the resulting clauses are
        used as input for a SAT solver provided by the PySAT package. These solvers include
        the ability to enumerate solutions, and this method lifts that to match the PyPR
        interface for sat. 
        
        If more complicated sat-based procedures are needed, `tseytin` is exposed, and 
        more complicated instances can be created manually. However, because direct SAT
        solving and model enumeration are the most common applications, this function provides
        a much simpler interface, abstracting away details of the encoding. 

        :param solver_name: A string giving the name of a sat solver provided by PySAT
            which will be used as the solver, defaults to "cadical195", which we observed
            to work well experimentally.
        :type solver_name: str, optional
        :param verbose: if `True` print statistics and timings for debugging, defaults to False
        :type verbose: bool, optional

        :return: if the function is unsatisfiable, return None. otherwise, on each iteration,
            return a dictionary which maps variables to their boolean values in a satisfying 
            assignment. Any variables which don't appear in the dict are "don't care".
        :rtype: dict[int,bool] | None
        """
        raise NotImplementedError # defined in SAT.py

    def functionally_equivalent(self,
        other: "BooleanFunction"
    ) -> bool:
        """Determines if two functions have the same truth table.

        This function determines if two functions are equivalent
        semantically (as opposed to structurally). This is accomplished
        through the tseytin transform and a SAT solver, as the
        problem is NP-complete. For large functions this may be an expensive
        or infeasible task.

        :param BooleanFunction other: the function to compare to.         
        
        :return equivalent: A boolean representing whether or not the two functions 
            have the same truth table.
        """
        raise NotImplementedError # defined in SAT.py

    # Storage
    def generate_ids(self,
        previous_ids: dict[Any, int] | None = None,
        in_place: bool = True
    ) -> dict[Any, int]:
        """Generate a dictionary which maps each object to a unique ID.

        When called with no inputs, this function will create these IDs from scratch. However,
        the function can also recieve the output of previous calls as input, and will continue 
        to add in new IDs, reusing IDs where possible. This means the same function can also be
        used to add onto existing output. For example:
        ```python
        ids = object_1.generate_ids()
        ids = object_2.generate_ids(ids)
        ids = object_3.generate_ids(ids)
        ...
        ```
        For efficiency reasons, this method both mutates the input and returns the mutated
        list of ids by default. This can be disabled by setting the parameter `in_place=False`

        :param previous_ids: the output of previous calls to the function
        :type previous_ids: dict[Any, int] | None
        :param in_place: whether the input list of ids is mutated in place or not. If false, a copy
            is created and returned. If true (by default) the data is modified in place with no copy.
        :type in_place: bool
        :return: A dict which maps each node to a unique id
        :rtype: dict[Any, int]
        """
        if not previous_ids:
            ids = {}
            next_available_index = 0
        elif in_place:
            ids = previous_ids
            next_available_index = max(previous_ids.values()) + 1
        else:
            # shallow copy to maintain objects, but new id dict
            ids = {k:v for k,v in previous_ids.items()}
            next_available_index = max(previous_ids.values()) + 1

        # nodes given ids by an earlier call keep them, and so does everything beneath them
        for node in self.postorder(stop=ids.__contains__):
            ids[node] = next_available_index
            next_available_index += 1

        return ids

    def _generate_JSON_entry(self,
        node_ids: dict["BooleanFunction", int]
    ) -> dict[str, Any]:
        """Create a JSON entry for a given object

        This should return a JSON encoding with all the data necessary to recreate the object
        to an acceptable degree, with a convention matching the corresponding `parse_JSON_entry`.
        This convention is dependent on the class, and this is the pair of methods to overwrite
        to implement JSON serialization for a custom object.

        For `BooleanFunction` specifically, the `args` attribute is stored as a list of ids,
        each pointing to a previously stored `BooleanFunction`. The `_compiled` function is 
        hard to store (and easy to regenerate with `.compile()`), and so it is simply not stored,
        and is viewed as an acceptable loss. The rest of the attributes in `__dict__` are copied
        with no change. 

        :param node_ids: A dictionary which maps each object to a unique id.
        :type node_ids: dict[Any, int]
        :return: A dictionary which represents the JSON data for one object
        :rtype: dict[str, Any]
        """
        # copy class name and non-nested data
        JSON_data = self.__dict__.copy()
        # use refs for children/nested data:
        if 'args' in JSON_data:
            JSON_data['args'] = [node_ids[arg] for arg in self.args]
        # ignore the compiled version (not serializable)
        if '_compiled' in JSON_data:
            del JSON_data['_compiled']

        return JSON_data

    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str,Any],
        parsed_functions: list["BooleanFunction | None"]
    ) -> Self:
        """Parse a JSON entry back into a given object.

        This method is able to parse JSON generated by `_generate_JSON_entry` for the
        corresponding class, according to some convention. This convention may be different,
        and is decided by the class implementer.

        For `BooleanFunction` specifically, the `args` attribute is stored as a list of ids,
        each pointing to a previously stored `BooleanFunction`. The `_compiled` function is 
        hard to store (and easy to regenerate with `.compile()`), and so it is simply not stored,
        and is viewed as an acceptable loss. The rest of the attributes in `__dict__` are copied
        with no change.

        :param object_data: A dictionary which contains the fields and data of the
            original node object, as generated by `_generate_JSON_entry`.
        :type object_data: dict[str, Any]
        :param parsed_objects: A list which contains the previously parsed objects.
            This can be used to get references to previously stored items, allowing the
            serialization methods to connect the parsed objects together in complex ways.
        :type parsed_objects: list[Any | None]
        :return: The parsed object, with data matching the JSON.
        :rtype: Self
        """
        # intantiate new object:
        new_node = object.__new__(cls)
        for key,value in object_data.items():

            # Use previously parsed functions for args
            if key == 'args':
                new_node.args = tuple([parsed_functions[child_id] for child_id in value])

            # for other fields, just set directly
            else:
                setattr(new_node,key,value)

        return new_node

    def to_JSON(self) -> dict[str,Any]:
        """An alias for `PyPR.JSON_Serialization.generate_JSON(fn)`
        
        This can be used in conjunction with `from_JSON` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to use the full 
        `generate_JSON`/`parse_JSON` methods in PyPR.JSON_Serialization. Check the docstrings on these
        methods for more information on usage and output.

        :return: A JSON object which encodes the input function
        :rtype: dict[str,Any]
        """
        return PyPR.JSON_Serialization.generate_JSON(self)

    @classmethod
    def from_JSON(cls,
        json_object: dict[str,Any]
    ) -> Self:
        """An alias for `PyPR.JSON_Serialization.parse_JSON(json_object)[0]`
        
        This can be used in conjunction with `to_JSON` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to use the full 
        `generate_JSON`/`parse_JSON` methods in PyPR.JSON_Serialization. Check the docstrings on these
        methods for more information on usage and output.

        Although there is no functional difference between `X.from_JSON` and `Y.from_JSON` for two
        classes (`X` and `Y`) which are both serializable, the class you call this method from is used
        to determine type hinting and to clarify the code. Therefore, I choose to throw an error if the
        json encodes a different class than the one you use to decode. This is mostly to enforce 
        readable code and good usage, and to make sure objects are interpreted correctly.

        :param json_object: A dictionary with the expected structure.
        :type json_object: dict[str,Any]
        :return: The BooleanFunction which was used to create the JSON.
        :rtype: BooleanFunction
        """
        return_idx = json_object['return order'][0]
        json_class = json_object['objects'][return_idx]['class']
        subclasses = {
            str(cls)[8:-2] for cls in
            PyPR.JSON_Serialization.all_subclasses(cls)
        }

        if json_class not in subclasses:
            raise ValueError(
                f"JSON encodes {json_class}, which is not " +
                f"a subclass of class {str(cls)[8:-2]}"
            )

        return PyPR.JSON_Serialization.parse_JSON(json_object)[0]

    def to_file(self,
        filename: str
    ) -> None:
        """Writes the output of fn.to_JSON to a file with the given filename.
        
        This can be used in conjunction with `from_file` to reduce verbosity and improve readibility
        when you only want to store/parse one object. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manage I/O manually and use 
        the full `generate_JSON`/`parse_JSON` methods in PyPR.JSON_Serialization. Check the docstrings\
        on these methods for more information on usage and output.

        :param filename: A string which will be used as the name of the generated file 
            (must end with the `.json` file extension)
        :type filename: str
        """
        # json files only:
        if filename[-5:] != ".json":
            raise ValueError("Filename must end with the \".json\" file extension")

        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    @classmethod
    def from_file(cls,
        filename: str
    ) -> Self:
        """Reads a single function from the file with the given filename.
        
        This can be used in conjunction with `to_file` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manage I/O manually and use 
        the full `generate_JSON`/`parse_JSON` methods. Check the docstrings on these methods for more 
        information on usage and output.

        Although there is no functional difference between `X.from_file` and `Y.from_file` for two
        classes (`X` and `Y`) which are both serializable, the class you call this method from is used
        to determine type hinting and to clarify the code. Therefore, I choose to throw an error if the
        json encodes a different class than the one you use to decode. This is mostly to enforce 
        readable code and good usage, and to make sure objects are interpreted correctly.

        :param filename: A string which gives the name of the file to read.
        :type filename: str
        """
        with open(filename, 'r') as f:
            return cls.from_JSON(json.loads(f.read()))


    # statistics and properties:
    def is_leaf(self) -> bool:
        """Returns `True` for leaf nodes (no children), and `False` otherwise.

        A leaf node is one which has no children. In the builtin functions, `VAR` and `CONST`
        always return True, while the gates (e.g. `XOR`, `AND`, etc.) always return False.

        We techinically don't force these gates to have arguments (because we trust the user),
        but it is possible to cause errors by allowing a gate as a leaf, and so it is advised against.

        :return: `True` for leaf nodes (no children), and `False` otherwise.
        :rtype: bool
        """
        return False

    def max_idx(self) -> int:
        """Returns the maximum index used in a variable in the function.

        :return: Returns the maximum index used in a variable in the function.
        :rtype: int
        """
        highest = -1
        for node in self.postorder():
            if node.is_leaf():
                highest = max(highest, node.max_idx())
        return highest

    def idxs_used(self) -> set[int]:
        """Return the set of indices used in variables in the function.

        :return: Return the set of indices used in variables in the function.
        :rtype: set[int]
        """
        indices: set[int] = set()
        for node in self.postorder():
            if node.is_leaf():
                indices |= node.idxs_used()
        return indices

    def num_nodes(self) -> int:
        """Return the number of nodes comprising the input function.

        :return: Return the number of nodes comprising the input function.
        :rtype: int
        """
        return sum(1 for _ in self.postorder())

    def component_count(self) -> dict[str,int]:
        """Return a dict which counts the occurrences of each class in the function DAG

        This is similar to `num_nodes`, but provides a more detailed description. The `__name__`
        attribute of each node to increment the corresponding count in the returned dictionary.

        :return: a dict which counts the occurrences of each class in the function DAG
        :rtype: dict[str,int]
        """
        components = {}
        for node in self.postorder():
            name = type(node).__name__
            if name in components:
                components[name] += 1
            else:
                components[name] = 1
        return components
