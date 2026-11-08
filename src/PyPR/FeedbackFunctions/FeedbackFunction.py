
from collections.abc import Callable
from copy import deepcopy
from functools import cached_property
from typing import Any, Self

import numba

# for compiling to python
import numpy as np  # noqa: F401 -- used by the source compile() execs
from numba.core import types as nb_types

from PyPR.JSON_Serialization import Serializable

from PyPR.BooleanLogic import VAR, BooleanANF, BooleanFunction
from PyPR.BooleanLogic.Latex import LatexStyle, fill_template, partial_name
from PyPR.BooleanLogic.Optimization import operation_count
from PyPR.BooleanLogic.Optimization import optimize as optimize_functions

u8 = nb_types.u8
void = nb_types.void



class FeedbackFunction(Serializable):
    fn_list: list[BooleanFunction]
    size: int
    _compiled: Any
    _compiled_inplace: Any

    def __init__(self, fn_list):
        #convert update to a list of ANF<int> objects:
        self.fn_list = fn_list
        self.size = len(fn_list)
        self._compiled = None
        self._compiled_inplace = None

    def __deepcopy__(self, memo: dict[int, Any]) -> Self:
        """Create a copy which shares no mutable state with the original.

        Every attribute is deep-copied through the shared `memo`, so anything
        reachable twice -- a node referenced by several bits, or an attribute
        pointing into `fn_list` -- is copied once and stays shared in the copy.
        The compiled functions are the exception: they are immutable numba
        dispatchers generated from `fn_list`, which the copy reproduces
        exactly, so the copy shares them and runs compiled without recompiling.

        :param memo: The `copy.deepcopy` memo, mapping `id` of an original to its copy.
        :type memo: dict[int, Any]
        :return: A copy of the feedback function.
        :rtype: Self
        """
        new_obj = object.__new__(type(self))
        memo[id(self)] = new_obj
        for name, value in self.__dict__.items():
            if name in ('_compiled', '_compiled_inplace'):
                new_obj.__dict__[name] = value
            else:
                new_obj.__dict__[name] = deepcopy(value, memo)
        return new_obj

    # copies are always deep: a shallow copy would share fn_list with the original
    def __copy__(self) -> Self:
        return deepcopy(self)

    def copy(self) -> Self:
        return deepcopy(self)

    def __getitem__(self, idx): return self.fn_list[idx]

    def __setitem__(self, idx, val): self.fn_list[idx] = val

    def __len__(self): return self.size

    # Strings:
    def __str__(self):
        outstr = ""
        for i in range(self.size-1,-1,-1):
            outstr += str(i) + "="
            outstr += str(self.fn_list[i]) + ";\n"
        return outstr[:-1]

    def pretty_str(self):
        outstr = ""
        for i in range(self.size-1,-1,-1):
            outstr += f"Bit {i} updates according to:\n"
            outstr += self.fn_list[i].pretty_str() + "\n\n\n"
        return outstr[:-3]

    def dense_str(self):
        outstr = ""
        for i in range(self.size-1,-1,-1):
            outstr += f"{i} = "
            outstr += self.fn_list[i].dense_str() + ";\n"
        return outstr[:-1]

    def anf_str(self):
        outstr = ""
        for i in range(self.size-1,-1,-1):
            outstr += str(i) + "="
            outstr += self.fn_list[i].anf_str() + ";\n"
        return outstr[:-1]


    # Convenient Manipulations
    def flip(self):
        """Reverse the bit labelling: bit i becomes bit size-1-i.

        This relabels the register rather than changing what it computes, but it
        does replace `fn_list`, so anything derived from that list is discarded
        along with it. Note it does **not** touch `is_inverted`: flipping is not
        inversion, and a subclass whose `invert` reads that flag to decide which
        end the shift discards (see `Fibonacci._inverse_feedback`) will read it
        as describing the pre-flip layout.
        """
        new_indices = {i: self.size-1-i for i in range(self.size)}
        self.fn_list = [f.remap_indices(new_indices) for f in self.fn_list][::-1]

        # derived caches are computed from fn_list, which was just replaced
        for name in list(self.__dict__):
            if isinstance(getattr(type(self), name, None), cached_property):
                del self.__dict__[name]


    # Storage
    def generate_ids(self,
        previous_ids: dict[Any, int] | None = None,
        in_place = True
    ) -> dict[Any, int]:
        """Generate a dictionary which maps each object to a unique ID.

        When called with no inputs, this function will create these IDs from scratch. However,
        the function can also recieve the output of previous calls as input, and will continue 
        to add in new IDs, reusing IDs where possible. This means the same function can also be
        used to add onto existing output. For example:

        .. code-block:: python

            ids = object_1.generate_ids()
            ids = object_2.generate_ids(ids)
            ids = object_3.generate_ids(ids)
            ...

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
        elif in_place:
            ids = previous_ids
        else:
            # shallow copy to maintain objects, but new id dict
            ids = {k:v for k,v in previous_ids.items()}

        # create ids for all functions:
        for fn in self.fn_list:
            ids = fn.generate_ids(ids)

        # lastly, add id for self -- unless an earlier call already gave it one:
        # renumbering would leave its old id empty and dangle references to it
        if self not in ids:
            ids[self] = max(ids.values()) + 1
        return ids

    def _generate_JSON_entry(self,
        ids: dict[Any, int]
    ) -> dict[str, Any]:
        """Create a JSON entry for a given object

        This should return a JSON encoding with all the data necessary to recreate the object
        to an acceptable degree, with a convention matching the corresponding `parse_JSON_entry`.
        This convention is dependent on the class, and this is the pair of methods to overwrite
        to implement JSON serialization for a custom object.

        For `FeedbackFunction` specifically, the `fn_list` attribute is stored as a list of ids,
        each pointing to a previously stored `BooleanFunction`. The `_compiled` and `_compiled_inplace`
        functions are hard to store (and easy to regenerate with `.compile()`), and so they are
        simply not stored, and are viewed as an acceptable loss. The rest of the attributes in 
        `__dict__` are copied with no change. 

        :param node_ids: A dictionary which maps each object to a unique id.
        :type node_ids: dict[Any, int]
        :return: A dictionary which represents the JSON data for one object
        :rtype: dict[str, Any]
        """
        # copy class name and non-nested data
        JSON_data = self.__dict__.copy()
        # fn_list:
        if 'fn_list' in JSON_data:
            JSON_data['fn_list'] = [ids[fn] for fn in self.fn_list]
                # ignore the compiled version (not serializable)
        if '_compiled' in JSON_data:
            del JSON_data['_compiled']
        if '_compiled_inplace' in JSON_data:
            del JSON_data['_compiled_inplace']

        # cached_property results are stored in __dict__ alongside real
        # attributes, and several are ndarrays or other non-JSON types
        # (update_matrix, CMPR's matrix caches). They are all derived from
        # fn_list, so dropping them costs nothing but a recomputation.
        for name in list(JSON_data):
            if isinstance(getattr(type(self), name, None), cached_property):
                del JSON_data[name]

        return JSON_data

    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str,Any],
        parsed_objects: list[Any | None]
    ) -> Self:
        """Parse a JSON entry back into a given object.

        This method is able to parse JSON generated by `_generate_JSON_entry` for the
        corresponding class, according to some convention. This convention may be different,
        and is decided by the class implementer.

        For `FeedbackFunction` specifically, the `fn_list` attribute is stored as a list of ids,
        each pointing to a previously stored `BooleanFunction`. The `_compiled` and `_compiled_inplace`
        functions are hard to store (and easy to regenerate with `.compile()`), and so they are
        simply not stored, and are viewed as an acceptable loss. The rest of the attributes in 
        `__dict__` are copied with no change. 

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
        new_obj = object.__new__(cls)
        for key,value in object_data.items():
            # Use previously parsed functions for args
            if key == 'fn_list':
                new_obj.fn_list = [parsed_objects[fn_id] for fn_id in value]

            # for other fields, just set directly
            else:
                setattr(new_obj,key,value)

        return new_obj

    # text generation
    # TODO: rename to match BF naming
    def write_VHDL(self, filename):
        #writes a VHDL file
        #Credit: Anna Hemingway

        # generate the assignments first: each subfunction split out of a bit
        # becomes a signal, and every one of them has to be declared up front
        overrides = {}
        assignments = []
        for i in range(self.size - 1, -1 , -1):
            assignments += self.fn_list[i].generate_VHDL(
                output_name = f"next_state({i})",
                array_name = "curr_state",
                subfunction_prefix = f"fn_{i}",
                overrides = overrides
            )

            # later bits refer to this bit's subfunctions by their signal
            for j, node in enumerate(self.fn_list[i].subfunctions()):
                if node not in overrides:
                    overrides[node] = f'fn_{i}_{j+1}'

        subfunction_signals = [
            line.split(" <= ")[0] for line in assignments
            if not line.startswith("next_state(")
        ]
        declarations = ""
        if subfunction_signals:
            declarations = f"    signal {', '.join(subfunction_signals)}: std_logic;\n"

        vhdl_str = f"""
library ieee;
use ieee.std_logic_1164.all;

entity fpr is
    port (
    i_clk :in std_logic;
    i_rst : in std_logic;
    i_seed_data: in std_logic_vector( {self.size - 1} downto 0);
    output: out std_logic_vector({self.size - 1} downto 0)
    );
end entity fpr;

architecture run of fpr is

    signal curr_state, next_state:std_logic_vector({self.size - 1} downto 0);
{declarations}

begin

    statereg: process(i_clk, i_rst)
    begin
        if (i_rst = '1') then
            curr_state <= i_seed_data;
        elsif (i_clk = '1' and i_clk'event) then
            curr_state <= next_state;
        end if;
    end process;\n"""

        vhdl_str += "\n    " + "".join(line + "\n    " for line in assignments)
        vhdl_str += """
    output <= curr_state;

end run;

"""
        with open(filename, "w") as f:
            f.write(vhdl_str)

    def generate_latex(self,
        output_name: str | Callable[..., str] = "x_{$bit}[t+1]",
        subfunction_name: str | Callable[..., str] = "g_{$bit,$index}",
        style: LatexStyle | None = None,
        environment: str | None = "align*",
        inline_subfunctions: bool = False
    ) -> str:
        """Write the register's update as a system of LaTeX equations, one per bit.

        Bits are written from the highest index down, matching `str` and
        `write_VHDL`. Each bit's function is rendered by
        `BooleanFunction.generate_latex`; a subfunction shared with a bit already
        written is referred to by the name it was given there rather than
        defined again.

        :param output_name: Template for the left-hand side of bit `$bit`'s
            equation, or a callable taking `bit`. Defaults to 'x_{$bit}[t+1]'.
        :type output_name: str | Callable[..., str]
        :param subfunction_name: Template naming the `$index`-th subfunction
            (1-based) split out of bit `$bit`'s function, or a callable taking
            `bit` and `index`. Defaults to 'g_{$bit,$index}'.
        :type subfunction_name: str | Callable[..., str]
        :param style: The notation to write in. Defaults to `LatexStyle` with
            state variables written `x_{$index}[t]`.
        :type style: LatexStyle | None
        :param environment: The environment to wrap the equations in, or None to
            return the lines joined by `\\\\` with no environment. Defaults to 'align*'.
        :type environment: str | None
        :param inline_subfunctions: If True, write every bit as a single
            expression with no subfunctions split out. Defaults to False.
        :type inline_subfunctions: bool
        :return: The LaTeX for the whole update.
        :rtype: str
        """
        if style is None:
            style = LatexStyle(variable="x_{$index}[t]")

        overrides: dict[BooleanFunction, str] = {}
        lines = []
        for i in range(self.size - 1, -1, -1):
            fn = self.fn_list[i]
            names = partial_name(subfunction_name, bit=i)
            lines += fn.generate_latex(
                output_name = fill_template(output_name, bit=i),
                subfunction_name = names,
                style = style,
                overrides = overrides,
                inline_subfunctions = inline_subfunctions,
            )

            # later bits refer to this bit's subfunctions by name
            if not inline_subfunctions:
                for j, node in enumerate(fn.subfunctions()):
                    if node not in overrides:
                        overrides[node] = fill_template(names, index=j+1)

        body = " \\\\\n".join(lines)
        if environment is None:
            return body
        return f"\\begin{{{environment}}}\n{body}\n\\end{{{environment}}}"

    # Compilation
    def optimize(self) -> "FeedbackFunction":
        """A feedback function with the same update, rewritten to take fewer operations.

        The bits are optimized together (`PyPR.BooleanLogic.Optimization.optimize`),
        so logic several bits have in common is one node shared between their
        functions. The update is the same on every state of bits; only the
        structure of `fn_list` differs. `generate_python(optimize=True)` and
        `compile` do this themselves, so this is for inspecting or exporting
        the optimized form.

        :return: A plain `FeedbackFunction` over the optimized functions. It
            does not keep the subclass or its attributes (the blocks of a CMPR, say).
        :rtype: FeedbackFunction
        """
        return FeedbackFunction(optimize_functions(self.fn_list))

    def _straight_line_program(self,
        functions: list[BooleanFunction],
        output_name: str,
        array_name: str
    ) -> list[str]:
        """Python statements computing `functions` from `array_name` into `output_name`.

        Every gate object is computed once, however many of the functions use
        it: a gate with several users gets a local, as does any expression long
        enough that nesting it further would only deepen the source.

        The layout otherwise follows the unoptimized code, because measurement
        favoured it. Bits are written from the highest index down, each as soon
        as its value is ready, and a state bit is read from the array at every
        use, not held in a local. Reading every state bit into a local first
        measured 5-20% slower on registers that are mostly shifts, which then
        hold more values at once than there are machine registers to keep them
        in. Like the unoptimized code, this writes before it has finished
        reading, so the two arrays must be distinct.

        :param functions: One function per bit, sharing nodes where they share logic.
        :type functions: list[BooleanFunction]
        :param output_name: The name of the array the new state is written to.
        :type output_name: str
        :param array_name: The name of the array holding the current state.
        :type array_name: str
        :return: The statements, in order, unindented.
        :rtype: list[str]
        """
        uses: dict[BooleanFunction, int] = {}
        order = []
        for function in functions:
            for node in function.postorder(stop=uses.__contains__):
                uses[node] = 0
                order.append(node)
        for node in order:
            for arg in node.args:
                uses[arg] += 1
        for function in functions:
            uses[function] += 1

        lines = []
        values: dict[BooleanFunction, str] = {}
        temporaries = 0
        for i in range(len(functions) - 1, -1, -1):
            for node in functions[i].postorder(stop=values.__contains__):
                text = node._generate_python(values, array_name)
                if not node.is_leaf() and (uses[node] > 1 or len(text) > 120):
                    values[node] = f"t{temporaries}"
                    temporaries += 1
                    lines.append(f"{values[node]} = {text}")
                else:
                    values[node] = text
            lines.append(f"{output_name}[{i}] = {values[functions[i]]}")
        return lines

    def generate_python(self,
        output_name: str = "next_state",
        array_name: str = "curr_state",
        optimize: bool = False
    ) -> list[str]:
        """Generates Python statements computing the whole update.

        The statements read the current state from the array `array_name` and
        write bit i of the next state to `output_name[i]`, highest bit first.
        They use only `^`, `&`, `|` and `1 - x`, so they are valid numba code
        over integer bits; `compile` is built on this. The two arrays must be
        distinct, since bits are written before all of them have been read.

        :param output_name: The name of the array the next state is written to.
            Defaults to 'next_state'.
        :type output_name: str
        :param array_name: The name of the array holding the current state.
            Defaults to 'curr_state'.
        :type array_name: str
        :param optimize: If False (the default), each bit's function is written
            as it stands, by `BooleanFunction.generate_python`: every
            subfunction of a bit gets a variable, and one already named by a
            higher bit is reused. If True, the bits are first rewritten
            together to share common logic and use fewer operations (see
            `optimize`), and every gate shared between any of them gets a
            variable. The update computed is the same on every state of bits,
            and `fn_list` is not changed. When the rewrite saves no operations
            the result is the same as with False.
        :type optimize: bool
        :return: One Python statement per line, dependencies first.
        :rtype: list[str]
        """
        if optimize:
            functions = optimize_functions(self.fn_list)
            # nothing saved: the functions as written are already as good, and
            # writing those keeps the code in step with fn_list
            if operation_count(functions) < operation_count(self.fn_list):
                return self._straight_line_program(functions, output_name, array_name)

        lines = []
        overrides: dict[BooleanFunction, str] = {}
        for i in range(self.size - 1, -1, -1):
            lines += self.fn_list[i].generate_python(
                output_name = f"{output_name}[{i}]",
                array_name = array_name,
                subfunction_prefix = f"fn_{i}",
                overrides = overrides
            )
            # a subfunction named here is reused by the bits written after it
            for j, node in enumerate(self.fn_list[i].subfunctions()):
                if node not in overrides:
                    overrides[node] = f"fn_{i}_{j+1}"
        return lines

    def compile(self, optimize: bool = True):
        """Compile the update with numba, and return the compiled function.

        Two functions are built and stored: `_compiled(state)` returns the next
        state as a new array, and `_compiled_inplace(state, buffer)` writes it
        into `buffer`, which is what `FeedbackRegister` clocks with. Both are
        the statements of `generate_python`, wrapped in a function.

        :param optimize: Passed to `generate_python`. True (the default)
            compiles the optimized form of the update when it saves operations:
            updates with shared or redundant logic then clock faster, and
            compile faster too, there being less code to compile, while an
            update the optimizer cannot improve is compiled exactly as written.
            False always compiles each bit's function as written.
        :type optimize: bool
        :return: The compiled `_compiled` function.
        :rtype: Callable
        """
        exec_locals = {}
        for name, signature, setup, output_name, result in (
            ("_compiled", "curr_state", "next_state = np.zeros_like(curr_state)", "next_state", "return next_state"),
            ("_compiled_inplace", "curr_state,output_buffer", "", "output_buffer", "return"),
        ):
            lines = self.generate_python(output_name, "curr_state", optimize=optimize)
            lines = ([setup] if setup else []) + lines + [result]
            exec_str = f"def {name}({signature}):\n" + "".join(f"    {line}\n" for line in lines)
            exec(exec_str, globals(), exec_locals)

        # u8[:](u8[:])
        self._compiled = numba.njit()(
            exec_locals["_compiled"]
        )

        # void(u8[:],u8[:])
        self._compiled_inplace = numba.njit()(
            exec_locals["_compiled_inplace"]
        )

        return self._compiled

    # Function unrolling (possibly remove)
    def iterator(self, n):
        fns = [VAR(i) for i in range(self.size)]
        yield fns

        for i in range(1,n+1):
            fns = [self.fn_list[b].compose(fns) for b in range(self.size)]
            yield fns


    # Statistics:
    def gate_summary(self):
        # get and merge counts from all fns
        dicts = [f.component_count() for f in self.fn_list]
        unified_keys = set.union(*(set(d.keys()) for d in dicts))
        output = {}
        for key in unified_keys:
            output[key] = 0
            for d in dicts:
                # 0 as a default value (if key not in d)
                output[key] += d.get(key, 0)
        return output

    def is_linear(self, allow_affine: bool = False) -> bool:
        """Whether every bit's update is a linear (or affine) function of the state.

        Decided on each bit's ANF, not on the gate types: ANF-derived DAGs wrap
        every monomial in an AND, even single-variable ones, so a linear register
        can contain AND nodes. An update is affine when every monomial has degree
        at most 1, and linear when it is affine and has no constant monomial
        (the empty term), so that it maps the zero state to zero.

        :param allow_affine: Whether to accept a constant term, i.e. test for affine
            rather than linear updates.
        :type allow_affine: bool
        :return: True if every bit's update is linear (affine, if `allow_affine`).
        :rtype: bool
        """
        for f in self.fn_list:
            anf = BooleanANF.from_BooleanFunction(f)
            if anf.degree() > 1:
                return False
            if not allow_affine and frozenset() in anf.terms:
                return False
        return True
