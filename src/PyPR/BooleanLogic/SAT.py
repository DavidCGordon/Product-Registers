"""SAT encoding and solving for Boolean functions.

The Tseytin walks (`tseytin`, `tseytin_labels`, `tseytin_clauses`) turn a
BooleanFunction DAG into CNF, one node at a time: every node is given its own
solver variables and emits the clauses that implement it. The solver entry
points (`sat`, `enum_models`, `functionally_equivalent`) are built on those
walks, and all six are attached to `BooleanFunction` as methods at the end of
the encoding section.

Fusion nodes (`TseytinFuse`) come after: a way to control the CNF a subtree
emits without changing the function it computes.

See `docs/architecture/SAT Encoding.md` for the per-node invariant the walks
rely on and fusion nodes exploit.
"""
from collections.abc import Callable, Iterator, Mapping
from typing import Any, Self

from pysat.formula import CNF
from pysat.solvers import Solver

from PyPR.BooleanLogic.BooleanFunction import BooleanFunction
from PyPR.BooleanLogic.FunctionInputs import VAR
from PyPR.BooleanLogic.Gates import AND, OR, XOR
from PyPR.BooleanLogic.Latex import LatexStyle, LatexTerm


def tseytin(self,
    prev_clauses: list[tuple[int, ...]] | None = None,
    prev_node_labels: dict[BooleanFunction,list[int]] | None = None,
    prev_variable_labels: dict[int,int] | None = None,
) -> tuple[
    list[tuple[int, ...]],
    dict[BooleanFunction,list[int]],
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
    clauses: dict[tuple,None]
    if not prev_clauses: clauses = dict.fromkeys([(1,)])
    else: clauses = dict.fromkeys([tuple(x) for x in prev_clauses])

    node_labels, variable_labels = self.tseytin_labels(
        prev_node_labels, prev_variable_labels
    )
    clauses.update(dict.fromkeys(self.tseytin_clauses(node_labels)))
    return list(clauses.keys()), node_labels, variable_labels

def tseytin_labels(self,
    node_labels: dict[BooleanFunction,list[int]] | None = None,
    variable_labels: dict[int,int] | None = None
) -> tuple[
    dict[BooleanFunction,list[int]],
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
    # initialize index maps if needed
    if node_labels == None and variable_labels == None:
        variable_labels = {}
        node_labels = {}

    # if only 1 is passed in, raise an error
    elif node_labels == None:
        raise ValueError("Missing node labels")
    elif variable_labels == None:
        raise ValueError("Missing variable labels")

    # Every node claims its result wire last and highest, so one past the
    # highest wire in the map is always free. Index 1 is reserved: `tseytin`
    # asserts the unit clause (1,), and CONST labels itself [1] or [-1], so a
    # constant is pinned by that clause and needs none of its own.
    # A map holding only CONST(0) has -1 as its highest label, which would put
    # the next wire at 0 -- not a literal at all, since 0 ends a DIMACS clause.
    next_idx = 2 if not node_labels else max(
        2, max(max(labels) for labels in node_labels.values()) + 1
    )

    # a node carried in from a previous call keeps its wires, and so does
    # everything beneath it -- this is what lets labels be threaded across
    # several circuits so shared subexpressions reuse variables
    for node in self.postorder(stop=node_labels.__contains__):
        next_idx = node._tseytin_labels(node_labels, variable_labels, next_idx)

    return node_labels,variable_labels

def tseytin_clauses(self,
    label_map: dict[BooleanFunction, list[int]]
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
    # distinct nodes can emit the same clause (two CONSTs pinning one literal);
    # the dict keeps one copy of each, in the order they were first emitted
    clauses: dict[tuple,None] = {}
    for node in self.postorder():
        clauses.update(dict.fromkeys(node._tseytin_clauses(label_map)))
    return list(clauses.keys())

def satisfiable(self,
    solver_name: str = "cadical195",
    verbose: bool = False
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
    clauses, node_map, var_map = self.tseytin()
    clauses += [(node_map[self][-1],)]
    num_variables = node_map[self][-1] + 1
    num_clauses = len(clauses)

    if verbose:
        print("Tseytin finished")
        print(f'Number of variables: {num_variables}')
        print(f'Number of clauses: {num_clauses}')

    cnf = CNF(from_clauses=clauses)
    with Solver(name = solver_name, bootstrap_with=cnf, use_timer=True) as solver:
        satisfiable = solver.solve()
        assignments: Any = solver.get_model()

    if verbose:
        print(solver.time())

    if satisfiable:
        return {k: (assignments[v-1]>0) for k,v in var_map.items()}
    else:
        return None

def enumerate_models(self,
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
    clauses, node_map, var_map = self.tseytin()
    clauses += [(node_map[self][-1],)]
    num_variables = len(node_map)
    num_clauses = len(clauses)
    cnf = CNF(from_clauses=clauses)

    if verbose:
        print(cnf.nv, len(cnf.clauses))
        print("Tseytin finished")
        print(f'Number of variables: {num_variables}')
        print(f'Number of clauses: {num_clauses}')

    with Solver(name = solver_name, bootstrap_with=cnf, use_timer=True) as solver:
        for assignment in solver.enum_models(): # type: ignore (this is from bad typing in pysat)
            yield {k: (assignment[v-1]>0) for k,v in var_map.items()}

def functionally_equivalent(self,
    other: BooleanFunction
) -> bool:
    """Decide whether two Boolean functions have the same truth table.

    Two functions f and g agree on every input exactly when f XOR g is
    unsatisfiable, so this builds the Tseytin encoding of `XOR(self, other)`
    and asks a SAT solver for a satisfying assignment: none means equivalent,
    and any one is an input on which they differ. The comparison is semantic,
    not structural. Deciding it is coNP-complete in general, so for large
    functions it may be expensive.

    :param other: The function to compare against.
    :type other: BooleanFunction

    :return: `True` if the two functions agree on every input.
    :rtype: bool
    """
    return ((satisfiable(XOR(self,other))) == None)

# add functions to BooleanFunction class
BooleanFunction.tseytin = tseytin
BooleanFunction.tseytin_labels = tseytin_labels
BooleanFunction.tseytin_clauses = tseytin_clauses
BooleanFunction.sat = satisfiable
BooleanFunction.enum_models = enumerate_models
BooleanFunction.functionally_equivalent = functionally_equivalent


# ── fusion nodes ─────────────────────────────────────────────────────────────
# A `_TseytinFuse` wraps a template function and the arguments it is applied to,
# and presents itself to the encoder as a single node. Because it owns the whole
# domain it fuses over, it can emit any CNF it likes for that domain -- including
# one with fewer wires and fewer clauses than expanding the subtree gate by gate
# would produce. It changes the encoding, never the function.
#
# Fusion is never automatic. A caller builds these nodes deliberately, which
# keeps the choice of fusion domain explicit and keeps every other node encoding
# exactly as it did before.
#
# See `BooleanFunction._tseytin_labels` and `._tseytin_clauses` for the contract
# these nodes implement.

class _TseytinFuse(BooleanFunction):
    """A subtree presented to the SAT encoder as one node.

    The stored `template` is over `VAR(0) .. VAR(n-1)`, and the node's
    arguments supply those inputs positionally, so the node denotes exactly
    `template.compose({i: args[i]})`. That equivalence is the whole semantics, and
    `expand()` produces it.

    What changes is the encoding. An ordinary gate hands the walk one wire per
    intermediate and relates them pairwise; a fusion node claims whatever wires
    its own encoding needs and emits the clauses for the entire domain at once.
    Wires interior to the domain therefore stop being nodes, but they are still
    published in this node's label list, so a caller can still read or constrain
    every wire in the solution -- it just has no node object to key them by.

    :ivar template: The template this node encodes, over `VAR(0) .. VAR(n-1)`.
    :vartype template: BooleanFunction
    :ivar args: The arguments supplying the template's inputs, in order.
    :vartype args: tuple[BooleanFunction, ...]
    """

    template: BooleanFunction

    # the template is a BooleanFunction of its own, stored by reference
    _JSON_references = ("template",)

    def __init__(self, template: BooleanFunction, *args: BooleanFunction) -> None:
        """Wrap `template` applied to `args`.

        :param template: A template over `VAR(0) .. VAR(n-1)`.
        :type template: BooleanFunction
        :param args: One argument per template input, in index order.
        :type args: BooleanFunction
        :raises ValueError: If the template reads an input with no argument.
        """
        # An input the arguments don't supply would otherwise be encoded as a
        # free interior wire -- a different function from the one `eval` sees,
        # which fails on the missing index instead.
        unsupplied = sorted(i for i in template.idxs_used() if not 0 <= i < len(args))
        if unsupplied:
            raise ValueError(
                f"the template reads inputs {unsupplied}, but only "
                f"{len(args)} argument(s) were given"
            )

        # A copy, so the node cannot be changed out from under itself by a
        # caller still holding the template, and so two nodes built from one
        # factory do not share interior objects -- the label map is keyed by
        # node identity, so shared interiors would collide. Sharing is meant to
        # happen at this node: reuse the _TseytinFuse itself and the walk gives
        # it one set of wires and one set of clauses, like any other node.
        self.template = template.__copy__()
        self.args = args
        self.arg_limit = None
        self._build_skeleton()

    def _generate_JSON_entry(self, node_ids: dict[Any, int]) -> dict[str, Any]:
        """Store the template (by reference) and the arguments, but not the skeleton.

        The skeleton fields are derived from the template at construction, and
        `_arg_wires` would not survive JSON anyway (its keys are ints), so they
        are rebuilt on parsing instead of stored.

        :param node_ids: A map from each object to its id.
        :type node_ids: dict[Any, int]
        :return: The node's data.
        :rtype: dict[str, Any]
        """
        data = super()._generate_JSON_entry(node_ids)
        for derived in ("_arg_wires", "_own_wires", "_skeleton"):
            data.pop(derived, None)
        return data

    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str, Any],
        parsed_functions: list[Any]
    ) -> Self:
        """Rebuild the node through its constructor, which rebuilds the skeleton.

        :param object_data: The data written for this node.
        :type object_data: dict[str, Any]
        :param parsed_functions: The objects rebuilt so far, indexed by id.
        :type parsed_functions: list[Any]
        :return: The rebuilt fusion node.
        :rtype: _TseytinFuse
        """
        template = parsed_functions[object_data["__refs__"]["template"]]
        args = [parsed_functions[arg_id] for arg_id in object_data["args"]]
        return cls(template, *args)

    def _copy(
        self,
        child_copies: dict[BooleanFunction, BooleanFunction]
    ) -> Self:
        """Rebuild this node over already-copied arguments.

        The template is not a child: it is this node's own data, so the default
        (which passes only the children to the constructor) would put the first
        child where the template goes. The constructor takes a copy of the
        template it is given, so the result shares nothing with this node and
        the copy is deep.

        :param child_copies: Copies of the node's arguments, keyed by original.
        :type child_copies: dict[BooleanFunction, BooleanFunction]
        :return: A fusion node over a copy of the template and the copied arguments.
        :rtype: _TseytinFuse
        """
        return type(self)(self.template, *(child_copies[arg] for arg in self.args))

    def expand(self) -> BooleanFunction:
        """The ordinary subtree this node stands for.

        Substituting the arguments into the template gives a function built
        from normal gates, which denotes the same thing and encodes the
        ordinary way. This is what a fused encoding must agree with, and what
        `verify` checks against.

        :return: The template with its inputs replaced by this node's arguments.
        :rtype: BooleanFunction
        """
        return self.template.compose({i: arg for i, arg in enumerate(self.args)})

    def verify(self) -> bool:
        """Whether this node's CNF encodes the same function as `expand()`.

        A fused encoding is written by hand, so nothing but a check like this
        establishes that it agrees with the subtree it replaces. The check is a
        SAT call on the XOR of the two, and it is deliberately *not* run during
        clause emission: emission is on the hot path, the cost is unbounded, and
        the check itself encodes this node, which would recurse.

        :return: True if the fused and expanded forms are equivalent.
        :rtype: bool
        """
        return self.functionally_equivalent(self.expand())

    def _eval(self, values: dict[BooleanFunction, Any], array: Any) -> Any:
        """Evaluate the template on the already-evaluated arguments.

        :param values: Evaluations of this node's children.
        :type values: dict[BooleanFunction, Any]
        :param array: The input container the walk was given.
        :type array: Any
        :return: The template's value on those arguments.
        :rtype: Any
        """
        return self.template.eval([values[arg] for arg in self.args])

    def _eval_ANF(self, values: dict[BooleanFunction, Any], array: Any) -> Any:
        """Evaluate the template over ANF-valued arguments.

        :param values: ANF evaluations of this node's children.
        :type values: dict[BooleanFunction, Any]
        :param array: The input container the walk was given.
        :type array: Any
        :return: The template's ANF value on those arguments.
        :rtype: Any
        """
        return self.template.eval_ANF([values[arg] for arg in self.args])

    def _template_code(
        self,
        render: Callable[[BooleanFunction, Mapping[BooleanFunction, Any]], Any],
        arg_values: Mapping[BooleanFunction, Any]
    ) -> Any:
        """Generate the template's expression over its arguments' expressions.

        Fusion only changes the SAT encoding; generated code is the template's
        own gates, with each input replaced by the code already generated for
        the argument that supplies it. The template is walked in post-order and
        every node other than an input is rendered the way the caller renders
        its own nodes.

        :param render: Renders one template node given its arguments' renderings,
            e.g. `lambda node, strings: node._generate_c(strings, array_name)`.
        :type render: Callable[[BooleanFunction, Mapping[BooleanFunction, Any]], Any]
        :param arg_values: The renderings already made for this node's arguments.
        :type arg_values: Mapping[BooleanFunction, Any]
        :return: One expression computing this node.
        :rtype: Any
        """
        values: dict[BooleanFunction, Any] = {}
        for node in self.template.postorder():
            if isinstance(node, VAR):
                values[node] = arg_values[self.args[node.index]]
            else:
                values[node] = render(node, values)
        return values[self.template]

    def _generate_c(self, c_strings: Mapping[BooleanFunction, str], array_name: str) -> str:
        return self._template_code(lambda node, strings: node._generate_c(strings, array_name), c_strings)

    def _generate_VHDL(self, vhdl_strings: Mapping[BooleanFunction, str], array_name: str) -> str:
        return self._template_code(lambda node, strings: node._generate_VHDL(strings, array_name), vhdl_strings)

    def _generate_python(self, python_strings: Mapping[BooleanFunction, str], array_name: str) -> str:
        return self._template_code(lambda node, strings: node._generate_python(strings, array_name), python_strings)

    def _generate_latex(self, terms: Mapping[BooleanFunction, LatexTerm | str], style: LatexStyle) -> LatexTerm:
        return self._template_code(lambda node, inner: node._generate_latex(inner, style), terms)

    def _binarize(self, cache: dict[BooleanFunction, BooleanFunction]) -> "_TseytinFuse":
        """Binarize below the fusion boundary, but not across it.

        Decomposing the domain into two-input gates is the thing this node
        exists to avoid, so the template is left whole and only the arguments
        are replaced by their binarized forms.

        :param cache: Already-binarized nodes, keyed by original.
        :type cache: dict[BooleanFunction, BooleanFunction]
        :return: A fusion node over binarized arguments.
        :rtype: _TseytinFuse
        """
        return type(self)(self.template, *[cache[arg] for arg in self.args])

    def add_arguments(self, *new_args: BooleanFunction) -> None:
        """Not supported: the arguments are fixed by the template's arity.

        :param new_args: Unused.
        :type new_args: BooleanFunction
        :raises ValueError: Always.
        """
        raise ValueError(
            "a fusion node's arguments are fixed by the template it wraps; "
            "build a new node over the template you want instead"
        )

    def remove_arguments(self, *old_args: BooleanFunction) -> None:
        """Not supported: the arguments are fixed by the template's arity.

        :param old_args: Unused.
        :type old_args: BooleanFunction
        :raises ValueError: Always.
        """
        raise ValueError(
            "a fusion node's arguments are fixed by the template it wraps; "
            "build a new node over the template you want instead"
        )

    # ── encoding ─────────────────────────────────────────────────────────

    def _build_skeleton(self) -> None:
        """Encode the template once, in wire numbers relative to this node.

        The CNF a fusion node emits never changes shape -- only which solver
        variables it is written over -- so it is built here, at construction,
        and reduced to a substitution at encoding time. That leaves the walk
        doing no work per node beyond renumbering.

        The skeleton is stored over *placeholder* wires. `_arg_wires` maps an
        argument position to the placeholder standing for that argument's
        result, and `_own_wires` lists the placeholders this node will need
        solver variables for, its own result last. Constant literals are not
        placeholders and survive renumbering untouched.
        """
        template_labels, input_labels = self.template.tseytin_labels()
        self._arg_wires = {
            i: input_labels[i] for i in range(len(self.args)) if i in input_labels
        }
        argument_wires = set(self._arg_wires.values())
        result = template_labels[self.template][-1]

        direct = self._direct_skeleton(result)
        if direct is not None:
            self._own_wires, self._skeleton = direct
            return

        clauses = self.template.tseytin_clauses(template_labels)
        interior = {
            abs(literal)
            for clause in clauses for literal in clause
            if abs(literal) != 1 and abs(literal) not in argument_wires
        }

        # A template that is a bare input or constant carries its result on a
        # wire this node does not own, so it takes one of its own and asserts
        # the equivalence; every node has to publish a result wire it owns.
        if result in argument_wires or abs(result) == 1:
            owned_result = max(interior | argument_wires | {1}) + 1
            clauses = clauses + [(-owned_result, result), (owned_result, -result)]
            interior.add(owned_result)
            result = owned_result

        self._own_wires = sorted(interior - {result}) + [result]
        self._skeleton = clauses

    def _direct_skeleton(
        self,
        result: int
    ) -> tuple[list[int], list[tuple[int, ...]]] | None:
        """A fused skeleton for the templates that have an easy one, else None.

        An n-input AND or OR is definable in one wire and n+1 clauses, where
        chaining two-input gates spends n-1 wires and 3(n-1) clauses. The
        conjunction case is `o -> a_i` for each argument together with
        `(and a_i) -> o`; the disjunction case is its dual.

        This is the placeholder pair, not a general fusion engine: anything else
        falls back to the ordinary expansion, which is always correct.

        :param result: The placeholder the ordinary encoding gave the template's
            result, used only to pick a fresh placeholder clear of it.
        :type result: int
        :return: `(own wires, clauses)` over placeholders, or None.
        :rtype: tuple[list[int], list[tuple[int, ...]]] | None
        """
        template = self.template
        if not isinstance(template, (AND, OR)):
            return None
        # the template must be exactly the gate applied to its inputs in order
        if not template.args or len(template.args) != len(self.args):
            return None
        if not all(isinstance(a, VAR) and a.index == i for i, a in enumerate(template.args)):
            return None
        if len(self._arg_wires) != len(self.args):
            return None

        args = [self._arg_wires[i] for i in range(len(self.args))]
        output = max(args + [result]) + 1

        if isinstance(template, AND):
            clauses = (
                [(-output, arg) for arg in args]
                + [tuple([output] + [-arg for arg in args])]
            )
        else:
            clauses = (
                [(-arg, output) for arg in args]
                + [tuple([-output] + list(args))]
            )
        return [output], clauses

    def _tseytin_labels(
        self,
        node_labels: dict[BooleanFunction, list[int]],
        variable_labels: dict[int, int],
        next_idx: int
    ) -> int:
        """Claim one solver variable per wire the skeleton uses.

        The count was fixed at construction, so this is a block of consecutive
        variables and nothing is walked. They are this node's wires and no
        other's: the interior of a fused domain has no nodes to be keyed by, so
        it is published here and is opaque to the caller by design, with the
        node's result last as every node's is.

        :param node_labels: Node-to-wires map, extended in place.
        :type node_labels: dict[BooleanFunction, list[int]]
        :param variable_labels: Input-variable-to-wire map; untouched here.
        :type variable_labels: dict[int, int]
        :param next_idx: The next unused solver variable.
        :type next_idx: int
        :return: The next unused solver variable after this node's claim.
        :rtype: int
        """
        node_labels[self] = [next_idx + i for i in range(len(self._own_wires))]
        return next_idx + len(self._own_wires)

    def _tseytin_clauses(
        self,
        label_map: dict[BooleanFunction, list[int]]
    ) -> list[tuple[int, ...]]:
        """Renumber the skeleton onto the wires this encoding actually uses.

        Placeholders map to the arguments' result wires and to the block claimed
        in `_tseytin_labels`; a literal's sign is carried through, and constants
        are left alone.

        :param label_map: Wires for every node encoded so far, including this
            one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the template over those arguments.
        :rtype: list[tuple[int, ...]]
        """
        renumber = {
            placeholder: label_map[self.args[i]][-1]
            for i, placeholder in self._arg_wires.items()
        }
        renumber.update(zip(self._own_wires, label_map[self]))

        return [
            tuple(
                renumber[abs(literal)] * (1 if literal > 0 else -1)
                if abs(literal) != 1 else literal
                for literal in clause
            )
            for clause in self._skeleton
        ]


def TseytinFuse(template: BooleanFunction):
    """Turn a template into a node constructor that fuses it.

    `TseytinFuse(template)` behaves like a gate class: calling it with
    arguments builds a node denoting `template` applied to them, which the SAT encoder treats
    as a single unit::

        Maj = TseytinFuse(OR(AND(VAR(0), VAR(1)), AND(VAR(1), VAR(2))))
        fn = XOR(Maj(VAR(3), VAR(4), VAR(5)), VAR(6))

    :param template: A template over `VAR(0) .. VAR(n-1)`.
    :type template: BooleanFunction
    :return: A constructor taking the arguments to apply the template to.
    :rtype: Callable[..., _TseytinFuse]
    """
    def build(*args: BooleanFunction) -> _TseytinFuse:
        return _TseytinFuse(template, *args)
    return build
