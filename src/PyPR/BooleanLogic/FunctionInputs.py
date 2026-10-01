from collections.abc import Mapping
from typing import Any, Self

from PyPR.BooleanLogic.BooleanFunction import BooleanFunction, IndexableContainer
from PyPR.BooleanLogic.Latex import LatexStyle, LatexTerm


class CONST(BooleanFunction):
    # a constant may be a plain bit or an object such as a BooleanANF
    # (remap_constants puts those in); the latter is stored by reference
    _JSON_references = ("value",)

    def __init__(self, value):
        self.args = ()
        self.arg_limit = 0
        self.value = value

    def is_leaf(self) -> bool:
        return True
    def max_idx(self) -> int:
        return -1
    def idxs_used(self) -> set[int]:
        return set()


    def _eval(self,
        values: dict[BooleanFunction,Any],
        array: IndexableContainer[int, Any]
    ) -> Any:
        return self.value
    def _eval_ANF(self,
        values: dict[BooleanFunction,Any],
        array: IndexableContainer[int, Any]
    ) -> Any:
        return self.value


    def _generate_c(self,
        c_strings: Mapping[BooleanFunction, str],
        array_name: str
    ) -> str:
        return f"{self.value}"
    def _generate_VHDL(self,
        vhdl_strings: Mapping[BooleanFunction, str],
        array_name: str
    ) -> str:
        return f" '{self.value}' "
    def _generate_python(self,
        python_strings: Mapping[BooleanFunction, str],
        array_name: str
    ) -> str:
        return f"{self.value}"
    def _generate_latex(self,
        terms: Mapping[BooleanFunction, LatexTerm | str],
        style: LatexStyle
    ) -> LatexTerm:
        return style.render_constant(self.value)
    def dense_str(self) -> str:
        return f"CONST({self.value})"

    # this is using this input to allow unhashable constants
    # TODO: rework to be more elegant
    def _remap_constants(self,
        const_map: list[tuple[Any,Any]]
    ):
        for key,new_value in const_map:
            try:
                if self.value == key:
                    self.value = new_value
                    break
            except:
                pass
    def _remap_indices(self,
        index_map: IndexableContainer[int,int]
    ):
        pass
    def _shift_indices(self,
        shift_amount: int
    ):
        pass

    # overwriting BooleanFunction
    def __copy__(self) -> Self:
        return type(self)(self.value)

    # overwriting BooleanFunction
    def _compose(self,
        input_map: IndexableContainer[int,BooleanFunction],
        in_place:bool = False
    ) -> "CONST":
        if in_place:
            return self
        else:
            return CONST(self.value)

    # overwriting BooleanFunction
    def component_count(self) -> dict[str,int]:
        return {"CONST":1}

    #overwriting BooleanFunction
    def inputs(self) -> list[BooleanFunction]:
        return [self]

    def _binarize(self,
        cache: dict["BooleanFunction", "BooleanFunction"]
        ) -> Self:
        return self

    def _tseytin_labels(self,
        node_labels: dict['BooleanFunction',list[int]],
        variable_labels: dict[int,int],
        next_idx: int
    ) -> int:
        """A helper function to generate the labels for this leaf

        Generates the labels, and returns the next available index, which
        can be used for future calls to the same function. `node_labels` and
        `variable_labels` are updated inside the call

        :param node_labels: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used to
            convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        :type node_labels: dict[BooleanFunction,list[int]]
        :param prev_variable_labels: Because the input variables are always supposed to be the same (even
            in different VAR `nodes`), we use this dict to map input variables to the corresponding
            variables in the sat solver. This can be used to convert the sat solution back to a satisfying
            assignment of input variables, or to impose additional conditions based on extra information.
        :type prev_variable_labels: dict[int,int]
        :param next_idx: the next available index which can be used for this node
        :type next_idx: int
        :return: the next available index for future calls to use.
        :rtype: int
        """
        if self.value == 0:
            node_labels[self] = [-1]
        elif self.value == 1:
            node_labels[self] = [1]
        else:
            raise ValueError("Bad Constant")
        return next_idx

    def _tseytin_clauses(self,
        label_map: dict[BooleanFunction, list[int]]
    ) -> list[tuple[int, ...]]:
        """Return no clauses: a leaf's wire carries its value by definition.

        `_tseytin_labels` has already bound this leaf to a solver variable, and
        there is no relationship to assert about it -- a variable is free and a
        constant is pinned where it was labelled. This is the base case of the
        walk, reached like every other node rather than special-cased by it.

        :param label_map: Wires for every node encoded so far; unused here.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: An (empty) list of clauses
        :rtype: list[tuple[int, ...]]
        """
        return []

    def num_nodes(self) -> int:
        return 1




class VAR(BooleanFunction):
    def __init__(self, index):
        self.args = ()
        self.arg_limit = 0
        self.index = index

    def is_leaf(self):
        return True
    def max_idx(self):
        return self.index
    def idxs_used(self):
        return {self.index}


    def _eval(self,
        values: dict[BooleanFunction,Any],
        array: IndexableContainer[int, Any]
    ):
        return array[self.index]
    def _eval_ANF(self,
        values: dict[BooleanFunction,Any],
        array: IndexableContainer[int, Any]
    ):
        return array[self.index]


    def _generate_c(self,
        c_strings: Mapping[BooleanFunction, str],
        array_name: str
    ) -> str:
        return f"{array_name}[{self.index}]"
    def _generate_VHDL(self,
        vhdl_strings: Mapping[BooleanFunction, str],
        array_name: str
    ) -> str:
        return f"{array_name}({self.index})"
    def _generate_python(self,
        python_strings: Mapping[BooleanFunction, str],
        array_name: str
    ) -> str:
        return f"{array_name}[{self.index}]"
    def _generate_latex(self,
        terms: Mapping[BooleanFunction, LatexTerm | str],
        style: LatexStyle
    ) -> LatexTerm:
        return style.render_variable(self.index)
    # overwriting BooleanFunction str methods
    def pretty_lines(self, depth:int = 0) -> list[str]:
        return [f"VAR({self.index})"]
    def dense_str(self) -> str:
        return f"VAR({self.index})"


    # overwriting BooleanFunction
    def _remap_constants(self,
        const_map: list[tuple[Any,Any]]
    ):
        pass
    def _remap_indices(self,
        index_map: IndexableContainer[int,int]
    ):
        if self.index in index_map:
            self.index = index_map[self.index]
    def _shift_indices(self,
        shift_amount:int
    ):
        self.index = self.index + shift_amount

    # overwriting BooleanFunction
    def __copy__(self) -> "VAR":
        return VAR(self.index)

    # overwriting BooleanFunction
    def _compose(self,
        input_map: IndexableContainer[int,BooleanFunction],
        in_place: bool = False
    ) -> BooleanFunction:
        try:
            return input_map[self.index]
        except:
            if in_place:
                return self
            else:
                return VAR(self.index)


    # overwriting BooleanFunction
    def component_count(self) -> dict[str,int]:
        return {"VAR":1}

    # overwriting BooleanFunction
    def inputs(self) -> list[BooleanFunction]:
        return [self]

    def _binarize(self,
        cache: dict["BooleanFunction", "BooleanFunction"]
        ) -> "VAR":
        return self

    def _tseytin_labels(self,
        node_labels: dict[BooleanFunction,list[int]],
        variable_labels: dict[int,int],
        next_idx: int
    ) -> int:
        """A helper function to generate the labels for this leaf

        Generates the labels, and returns the next available index, which
        can be used for future calls to the same function. `node_labels` and
        `variable_labels` are updated inside the call

        :param node_labels: The Tseytin transform introduces new variables for each gate in the function.
            This dict maps every node to the list of variables which are used to encode it, and can be used to
            convert the sat solution back to values in the wires of the circuit, or to impose additional
            conditions based on extra information.
        :type node_labels: dict[BooleanFunction,list[int]]
        :param prev_variable_labels: Because the input variables are always supposed to be the same (even
            in different VAR `nodes`), we use this dict to map input variables to the corresponding
            variables in the sat solver. This can be used to convert the sat solution back to a satisfying
            assignment of input variables, or to impose additional conditions based on extra information.
        :type prev_variable_labels: dict[int,int]
        :param next_idx: the next available index which can be used for this node
        :type next_idx: int
        :return: the next available index for future calls to use.
        :rtype: int
        """
        if self.index in variable_labels:
            node_labels[self] = [variable_labels[self.index]]
            return next_idx

        variable_labels[self.index] = next_idx
        node_labels[self] = [next_idx]
        return next_idx + 1

    def _tseytin_clauses(self,
        label_map: dict[BooleanFunction, list[int]]
    ) -> list[tuple[int, ...]]:
        """Return no clauses: a leaf's wire carries its value by definition.

        `_tseytin_labels` has already bound this leaf to a solver variable, and
        there is no relationship to assert about it -- a variable is free and a
        constant is pinned where it was labelled. This is the base case of the
        walk, reached like every other node rather than special-cased by it.

        :param label_map: Wires for every node encoded so far; unused here.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: An (empty) list of clauses
        :rtype: list[tuple[int, ...]]
        """
        return []

    def num_nodes(self) -> int:
        return 1
