# TYPE ANNOTATIONS: TRUE
# DOCSTRINGS: FALSE

from typing import Self, Optional, Any, Protocol
from collections.abc import Iterator

# njit actually used (in compile) even if greyed out
from numba import njit

import json


class IndexableContainer[K,V](Protocol):
    def __getitem__(self, key: K, /) -> V: ...

class BooleanFunction:
    args: tuple["BooleanFunction", ...]
    arg_limit: Optional[int]

    def __init__(self,
        *args: "BooleanFunction",
        arg_limit: Optional[int] = None
    ):
        self.args = args
        self.arg_limit = arg_limit
    
    @classmethod
    def _copy(cls, 
        fn: "BooleanFunction", 
        child_copies: dict["BooleanFunction","BooleanFunction"] 
    ) -> Self:
        """Helper function which specifies how to copy a specific node in a 
        boolean function.

        Given a dictionary of child copies, form a new node and return it.
        This is similar to recursion but is used with iteration and memoization
        in `.copy()` and `.__copy__()`. This is also the function to override
        if you want custom behavior (such as copying custom attributes) for
        a child class which extends boolean function.

        :param BooleanFunction fn: the node you wish to return a copy of.
        :param dict child_copies: A dict containing copies of necessary nodes, 
            which are needed to make the new copy.                
        
        :return output: A copy of the input fn.
        """
        return cls(
            *(child_copies[arg] for arg in fn.args),
            arg_limit = fn.arg_limit
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
        is the same subclass as the input.       
        
        :return copy: a copy of the input fn.
        """
        copies: dict[Any,Any] = {}
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in copies:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                copies[curr_node] = type(curr_node)._copy(curr_node,copies)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                copies[curr_node] = curr_node.__copy__()
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return copies[self]
    
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
            self.args = tuple()
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
        subfuncs = []
        visited = set()
        stack: list[Any] = [self]
        last = None

        # used to sort subfuncs by when they appear:
        order = {}
        idx = 0
        
        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting an already visited gate:
            elif curr_node in visited and not curr_node.is_leaf():
                if curr_node not in subfuncs:
                    subfuncs.append(curr_node)
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
                visited.add(curr_node)
                order[curr_node] = idx
                idx += 1

                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                visited.add(curr_node)
                order[curr_node] = idx
                idx += 1
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return sorted(subfuncs, key = lambda x: order[x])

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
        leaves = []
        visited = set()
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in visited:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
                visited.add(curr_node)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                leaves.append(curr_node)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return leaves


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
        subfuncs = self.subfunctions() + [self]
        fn_strings = {root:f"(subfunction {i+1})" for i,root in enumerate(subfuncs)}
        fn_strings[self] = f"(Main Function)"
        pretty_strings = []

        for root in subfuncs:
            # no need for visited set (subfuncs is the same info):
            stack: list[Any] = [root]
            last = None

            pstr = f"{fn_strings[root]} = (\n"
            indent_lvl = 0

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                if curr_node in subfuncs and curr_node != root:
                    pstr += (
                        "   |" * indent_lvl + "   " + 
                        fn_strings[curr_node] + '\n'
                    )

                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    pstr += (
                        "   |" * indent_lvl + "   " +
                        curr_node.dense_str() + '\n'
                    )

                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    pstr += ("   |" * (indent_lvl-1) + "   )\n")
                    indent_lvl -= 1
                    last = stack.pop()
                    continue
                
                # before moving down to children:
                else:
                    # print prefix for functions with args:
                    pstr += (
                        "   |" * (indent_lvl) + "   "
                        f"{type(curr_node).__name__} (\n"
                    )

                    indent_lvl += 1
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            # finish string and add it in:
            pstr += ")\n"
            pretty_strings.append(pstr)

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
        fn_strings = {
            root:f"(subfunction {i+1})" 
            for i,root in enumerate(subfuncs)
        }
        
        # no need for visited set (subfuncs is the same info)
        stack: list[Any] = [self]
        last = None
        out_str = ""

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a subfunc:
            if curr_node in subfuncs:
                out_str += (fn_strings[curr_node] + ",")
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                out_str += (curr_node.dense_str() + ",")
                last = stack.pop()
                continue
                
            # moving up the tree after finishing children:
            elif last == False:
                # strip trailing comma
                if out_str[-1] == ',':
                    out_str = out_str[:-1]
                
                out_str += "),"
                last = stack.pop()
                continue
            
            # before moving down to children:
            else:
                # print prefix for functions with args:
                out_str += f"{type(curr_node).__name__}("

                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                continue

        # strip trailing comma
        if out_str[-1] == ',':
            out_str = out_str[:-1]
        return out_str

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

    def _generate_c(self,
        c_strings: dict["BooleanFunction", str],  
        array_name: str
    ) -> str:
        """Helper function which specifies how to generate VHDL for a BooleanFunction subclass

        Given a dictionary of child VHDL strings and the name of the evaluation array,
        form the VHDL for node and return it. This acts a hook, which is called in `generate_VHDL`
        and allows custom subclasses to generate valid VHDL as long as they have an implementation

        :param dict[BooleanFunction, str] vhdl_strings: a dictionary mapping child nodes
            to their corresponding VHDL strings.
        :param str array_name: What to use as the array name, if the evaluation array is referenced
            this is necessary so that `VAR` nodes can return the correct string.

        :return output: A string representing the computation for this node in VHDL.
        """
        raise NotImplementedError  # implemented for each subclass
            
    def generate_c(self,
        output_name: str = 'output',
        subfunction_prefix: str = 'fn', 
        array_name: str = 'array',
        overrides: dict["BooleanFunction", str]  = {}
    ) -> list[str]:
        """Generates valid python for the given function

        This will generate a list of strings which is valid python for a given function. there are
        various parameters that can be tweaked in order to allow functions to interact with
        the same or different variables in different ways. Although this is exposed to the 
        user, its a little difficult to use manually, and is mostly used elsewhere in the library
        (such as in `compile` method in the FeedbackFunction Class).

        :param str output_name: The string to be used in the generate python for 
            the output variable. The default value is 'output'. 
        :param str subfunction_prefix: In the generated python, subfunction variables
            will have the form {subfunction_prefix}_{index}. This string allows you
            to set the prefix. The default is 'fn'.
        :param str array_name: The string to be used for the name of the array holding
            the current state in the generated python. The default is 'array'.
        :param dict[BooleanFunction, str] overrides: A dictionary containing string overrides.
            if a node is in this dictionary, its corresponding override string will be used in
            place of the generated VHDL. Modify carefully, as this can cause the generated python
            to be invalid.

        :returns python_lines: A list containing several lines, each of which is a string 
        of valid python. These lines describe the computational DAG of the circuit.
        """
        subfuncs = self.subfunctions() + [self]
        fn_strings = {
            root:f"{subfunction_prefix}_{i+1}" 
            for i,root in enumerate(subfuncs)
        }

        fn_strings[self] = f"{output_name}"
        subfunction_lines = []
        c_strings = {}

        for root in subfuncs:
            # dont rederive an expression we already have:
            if root in overrides:
                continue


            # no need for visited set (subfuncs is the same info):
            stack: list[Any] = [root]
            last = None

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # override node already has a string
                elif curr_node in overrides:
                    c_strings[curr_node] = overrides[curr_node]
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                elif curr_node in subfuncs and curr_node != root:
                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    c_strings[curr_node] = curr_node._generate_c(c_strings,array_name)
                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    c_strings[curr_node] = curr_node._generate_c(c_strings,array_name)
                    last = stack.pop()
                    continue
            
                else:
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            subfunction_lines.append(f"{fn_strings[root]} = {c_strings[root]};")
            c_strings[root] = fn_strings[root]
        return subfunction_lines

    # def generate_tex(self):
        pass

    def _generate_VHDL(self,
        vhdl_strings: dict["BooleanFunction", str],  
        array_name: str
    ) -> str:
        """Helper function which specifies how to generate VHDL for a BooleanFunction subclass

        Given a dictionary of child VHDL strings and the name of the evaluation array,
        form the VHDL for node and return it. This acts a hook, which is called in `generate_VHDL`
        and allows custom subclasses to generate valid VHDL as long as they have an implementation

        :param dict[BooleanFunction, str] vhdl_strings: a dictionary mapping child nodes
            to their corresponding VHDL strings.
        :param str array_name: What to use as the array name, if the evaluation array is referenced
            this is necessary so that `VAR` nodes can return the correct string.

        :return output: A string representing the computation for this node in VHDL.
        """
        raise NotImplementedError  # implemented for each subclass

    def generate_VHDL(self,
        output_name: str = 'output',
        subfunction_prefix: str = 'fn', 
        array_name: str = 'array',
        overrides: dict["BooleanFunction", str]  = {}
    ) -> list[str]:
        """Generates valid VHDL for the given function

        This will generate a list of string which is valid VHDL for a given function. there are
        various parameters that can be tweaked in order to allow functions to interact with
        the same or different variables in different ways. Although this is exposed to the 
        user, its a little difficult to use manually, and is mostly used elsewhere in the library
        (such as in `write_VHDL` method in the FeedbackFunction Class).

        :param str output_name: The string to be used in the generate VHDL for 
            the output variable. The default value is 'output'. 
        :param str subfunction_prefix: In the generated VHDL, subfunction variables
            will have the form {subfunction_prefix}_{index}. This string allows you
            to set the prefix. The default is 'fn'.
        :param str array_name: The string to be used for the name of the array holding
            the current state in the generated VHDL. The default is 'array'.
        :param dict[BooleanFunction, str] overrides: A dictionary containing string overrides.
            if a node is in this dictionary, its corresponding override string will be used in
            place of the generated VHDL. Modify carefully, as this can cause the generated VHDL
            to be invalid.

        :returns VHDL_lines: A list containing several lines, each of which is a string 
        of valid VHDL. These lines describe the computational DAG of the circuit
        """
        subfuncs = self.subfunctions() + [self]
        fn_strings = {
            root:f"{subfunction_prefix}_{i+1}" 
            for i,root in enumerate(subfuncs)
        }

        fn_strings[self] = f"{output_name}"
        subfunction_lines = []
        vhdl_strings = {}

        for root in subfuncs:
            # dont rederive an expression we already have:
            if root in overrides:
                continue

            # no need for visited set (subfuncs is the same info):
            stack: list[Any] = [root]
            last = None

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # override node already has a string
                elif curr_node in overrides:
                    vhdl_strings[curr_node] = overrides[curr_node]
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                elif curr_node in subfuncs and curr_node != root:
                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    vhdl_strings[curr_node] = curr_node._generate_VHDL(vhdl_strings,array_name)
                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    vhdl_strings[curr_node] = curr_node._generate_VHDL(vhdl_strings,array_name)
                    last = stack.pop()
                    continue
            
                else:
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            subfunction_lines.append(f"{fn_strings[root]} <= {vhdl_strings[root]};")
            vhdl_strings[root] = fn_strings[root]
        return subfunction_lines

    def _generate_python(self,
        python_strings: dict["BooleanFunction", str],  
        array_name: str
    ) -> str:
        """Helper function which specifies how to generate python for a BooleanFunction subclass

        Given a dictionary of child python strings and the name of the evaluation array,
        form the python for node and return it. This acts a hook, which is called in `generate_python`
        and allows custom subclasses to generate valid python as long as they have an implementation

        :param python_strings: a dictionary mapping child nodes to their corresponding python strings.
        :type python_strings: dict[BooleanFunction, str] 
        :param array_name: What to use as the array name, if the evaluation array is referenced
            this is necessary so that `VAR` nodes can return the correct string.
        :type array_name:

        :return: A string representing the computation for this node in python.
        :rtype: str
        """
        raise NotImplementedError  # implemented for each subclass
        
    def generate_python(self,
        output_name: str = 'output',
        subfunction_prefix: str = 'fn', 
        array_name: str = 'array',
        overrides: dict["BooleanFunction", str]  = {}
    ) -> list[str]:
        """Generates valid python for the given function

        This will generate a list of strings which is valid python for a given function. there are
        various parameters that can be tweaked in order to allow functions to interact with
        the same or different variables in different ways. Although this is exposed to the 
        user, its a little difficult to use manually, and is mostly used elsewhere in the library
        (such as in `compile` method in the FeedbackFunction Class).

        :param str output_name: The string to be used in the generate python for 
            the output variable. The default value is 'output'. 
        :param str subfunction_prefix: In the generated python, subfunction variables
            will have the form {subfunction_prefix}_{index}. This string allows you
            to set the prefix. The default is 'fn'.
        :param str array_name: The string to be used for the name of the array holding
            the current state in the generated python. The default is 'array'.
        :param dict[BooleanFunction, str] overrides: A dictionary containing string overrides.
            if a node is in this dictionary, its corresponding override string will be used in
            place of the generated VHDL. Modify carefully, as this can cause the generated python
            to be invalid.

        :returns python_lines: A list containing several lines, each of which is a string 
        of valid python. These lines describe the computational DAG of the circuit.
        """
        subfuncs = self.subfunctions() + [self]
        fn_strings = {
            root:f"{subfunction_prefix}_{i+1}" 
            for i,root in enumerate(subfuncs)
        }

        fn_strings[self] = f"{output_name}"
        subfunction_lines = []
        py_strings = {}

        for root in subfuncs:
            # dont rederive an expression we already have:
            if root in overrides:
                continue


            # no need for visited set (subfuncs is the same info):
            stack: list[Any] = [root]
            last = None

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # override node already has a string
                elif curr_node in overrides:
                    py_strings[curr_node] = overrides[curr_node]
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                elif curr_node in subfuncs and curr_node != root:
                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    py_strings[curr_node] = curr_node._generate_python(py_strings,array_name)
                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    py_strings[curr_node] = curr_node._generate_python(py_strings,array_name)
                    last = stack.pop()
                    continue
            
                else:
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            subfunction_lines.append(f"{fn_strings[root]} = {py_strings[root]}")
            py_strings[root] = fn_strings[root]
        return subfunction_lines

    # def generate_tex(self):
        pass

    # def _generate_tex(self,
    #     tex_strings: dict["BooleanFunction", str],
    #     array_name: str
    # ) -> str:
    #     """Helper function which specifies how to generate tex for a BooleanFunction subclass

    #     Given a dictionary of child tex strings and the name of the evaluation array,
    #     form the VHDL for node and return it. This acts a hook, which is called in `generate_tex`
    #     and allows custom subclasses to generate valid tex as long as they have an implementation

    #     :param tex_strings: a dictionary mapping child nodes to their corresponding python strings.
    #     :type tex_strings: dict[BooleanFunction, str]
    #     :param array_name: What to use as the array name, if the evaluation array is referenced
    #         this is necessary so that `VAR` nodes can return the correct string.
    #     :type array_name: str

    #     :return output: A string representing the computation for this node in python.
    #     """
    #     raise NotImplementedError  # implemented for each subclass
        
    # def generate_tex(self,
    #     output_name: str = 'output',
    #     subfunction_prefix: str = 'fn', 
    #     array_name: str = 'array',
    #     overrides: dict["BooleanFunction", str]  = {}
    # ) -> list[str]:
    #     """Generates valid python for the given function

    #     This will generate a list of strings to help format a given function into tex. There are
    #     various parameters that can be tweaked to help make the output more flexible, but they are
    #     limited; for anything beyond relatively basic usage it's encouraged you write your own
    #     function based on the BooleanFunction structure (as opposed to trying to coerce this method
    #     into doing something it wasn't built for). 

    #     :param str output_name: The string to be used in the generate python for 
    #         the output variable. The default value is 'output'. 
    #     :param str subfunction_prefix: In the generated python, subfunction variables
    #         will have the form {subfunction_prefix}_{index}. This string allows you
    #         to set the prefix. The default is 'fn'.
    #     :param str array_name: The string to be used for the name of the array holding
    #         the current state in the generated python. The default is 'array'.
    #     :param dict[BooleanFunction, str] overrides: A dictionary containing string overrides.
    #         if a node is in this dictionary, its corresponding override string will be used in
    #         place of the generated VHDL. Modify carefully, as this can cause the generated python
    #         to be invalid.

    #     :returns python_lines: A list containing several lines, each of which is a string 
    #     of valid python. These lines describe the computational DAG of the circuit.
    #     """
    #     subfuncs = self.subfunctions() + [self]
    #     fn_strings = {
    #         root:f"{subfunction_prefix}_{i+1}" 
    #         for i,root in enumerate(subfuncs)
    #     }

    #     fn_strings[self] = f"{output_name}"
    #     subfunction_lines = []
    #     py_strings = {}

    #     for root in subfuncs:
    #         # dont rederive an expression we already have:
    #         if root in overrides:
    #             continue


    #         # no need for visited set (subfuncs is the same info):
    #         stack: list[Any] = [root]
    #         last = None

    #         while stack:
    #             curr_node = stack[-1]

    #             # dont interact with sentinel values
    #             if curr_node == False:
    #                 last = stack.pop()
    #                 continue

    #             # override node already has a string
    #             elif curr_node in overrides:
    #                 py_strings[curr_node] = overrides[curr_node]
    #                 last = stack.pop()
    #                 continue

    #             # hitting a subfunc contained in the current one:
    #             elif curr_node in subfuncs and curr_node != root:
    #                 last = stack.pop()
    #                 continue

    #             # hitting a leaf:
    #             elif curr_node.is_leaf():
    #                 py_strings[curr_node] = curr_node._generate_tex(py_strings,array_name)
    #                 last = stack.pop()
    #                 continue
                    
    #             # moving up the tree after finishing children:
    #             elif last == False:
    #                 py_strings[curr_node] = curr_node._generate_tex(py_strings,array_name)
    #                 last = stack.pop()
    #                 continue
            
    #             else:
    #                 # set up children to process:
    #                 stack.append(False) # sentinel value
    #                 for child in reversed(curr_node.args):
    #                     stack.append(child)
    #                 continue

    #         subfunction_lines.append(f"{fn_strings[root]} = {py_strings[root]}")
    #         py_strings[root] = fn_strings[root]
    #     return subfunction_lines


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
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in new_nodes:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                new_nodes[curr_node] = curr_node._binarize(new_nodes)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                new_nodes[curr_node] = curr_node.__copy__()
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

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
        const_map: Any
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
        constant_map: Any, 
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
        input_map: Any, 
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
        input_map: Any, 
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
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            if curr_node in new_nodes:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
            
                # new node:
                if in_place:
                    curr_node.args = tuple([new_nodes[arg] for arg in curr_node.args])
                    new_nodes[curr_node] = curr_node
                else:
                    new_nodes[curr_node] = type(curr_node)._copy(curr_node,new_nodes)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                new_nodes[curr_node] = curr_node._compose(input_map,in_place)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return new_nodes[self]
   
    def _merge_redundant(self,
        cache: dict["BooleanFunction","BooleanFunction"], 
        subfunctions: list["BooleanFunction"], 
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
        :param subfunctions: A list of the subfunctions of the root node on which 
            `merge_redundant` was called.
        :type subfunctions: list[BooleanFunction]
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
            return type(self)._copy(self, cache)
    
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
        subfunctions = self.subfunctions()

        new_nodes = {}
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            if curr_node in new_nodes:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
            
                # new node:
                new_nodes[curr_node] = curr_node._merge_redundant(
                    new_nodes, subfunctions, in_place = in_place,
                )

            # hitting a leaf:
            elif curr_node.is_leaf():
                new_nodes[curr_node] = curr_node
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

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
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in values:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                values[curr_node] = curr_node._eval(values, array)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                values[curr_node] = curr_node._eval(values, array)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

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
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in values:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                values[curr_node] = curr_node._eval_ANF(values, array)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                values[curr_node] = curr_node._eval_ANF(values, array)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

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
    @classmethod
    def tseytin_formula(cls, 
        output: int,
        *args: int
    ) -> list[tuple[int]]:
        """Given the variables, return the clauses associated with this node 

        Output is the int label for the output wire, which is given first as a
        convention. The rest of the args are also ints, but are variable args to
        allow classes to define their tseyting representation more flexibly.
        
        :raises NotImplementedError: If not implemented for the node class
        :return: A list of clauses encoding the wire relationship for the node
        :rtype: list[tuple[int]]
        """
        raise NotImplementedError
    
    @classmethod
    def tseytin_unroll(cls,
        gate_labels: list[int],
        arg_labels: list[int]
    ) -> list[tuple[int]]:
        """Given labels for the gate wires and nodes, unroll a binary tseytin formula
        to describe a multi-argument gate.

        This is used to expand gates with multiple arguments (such as large `ANDs` or `XORs`),
        while using the binary formula in `tseytin_formula`. 

        :param gate_labels: A list with labels for the gate wires (length `n-1`)
        :type gate_labels: list[int]
        :param arg_labels: A list with labels for the argument wires (length `n`)
        :type arg_labels: list[int]
        :return: A list of clauses encoding the wire relationship for the node
        :rtype: list[tuple[int]]
        :raises NotImplementedError: 
        """
        raise NotImplementedError

    def tseytin(self, 
        prev_clauses: list[tuple[int]] | None = None, 
        prev_node_labels: dict['BooleanFunction',list[int]] | None = None, 
        prev_variable_labels: dict[int,int] | None = None
    ) -> tuple[
        list[tuple[int]],
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
        :type prev_clauses: list[tuple[int]]
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
            list[tuple[int]],
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
    ) -> list[tuple[int]]:
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
        :rtype: list[tuple[int]],
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
        previous_ids: dict["BooleanFunction", int] | None = None
    ) -> dict["BooleanFunction", int]:
        """Generate a dictionary which maps each node to a unique ID.

        When called with no inputs, this function will create these IDs from scratch. However,
        the function can also recieve the output of previous calls as input, and will continue 
        to add in new IDs, reusing IDs where possible. This means the same function can also be
        used to add onto existing output. For example:
        ```python
        node_ids = function_1.generate_ids()
        node_ids = function_2.generate_ids(node_ids)
        node_ids = function_3.generate_ids(node_ids)
        ...
        ```

        :param previous_ids: the output of previous calls to the function
        :type previous_ids: dict[BooleanFunction, int]
        :return: A dict which maps each node to a unique id
        :rtype: dict[BooleanFunction, int]
        """
        if previous_ids:
            node_labels = previous_ids
            next_available_index = max(previous_ids.values()) + 1
        else:
            node_labels = {}
            next_available_index = 0
        
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in node_labels:
                last=stack.pop()
                continue

            # moving up the tree after finishing children or hitting a leaf:
            elif last == False or curr_node.is_leaf():
                if curr_node not in node_labels:
                    node_labels[curr_node] = next_available_index
                    next_available_index += 1
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return node_labels
    
    def _generate_JSON_entry(self,
        node_ids: dict["BooleanFunction", int]
    ) -> dict[str, Any]:
        """Create a JSON entry for a given node

        The default behavior is that data for the node is that most fields are simply stored
        in the "data" field of the dict as a nested dict which matches the object fields.
        The fields which dont map nicely to JSON values (i.e. `args` and potentially `_compiled`)
        are the only ones which are handled differently. args uses the node ids instead, and
        the `_compiled` field is simply dropped, as it is impractical and counterintuitive to store.
        Other node types can override this to create custom JSON, if they have other fields or need
        to be handled in special ways (e.g. `VAR`,`CONST`, or custom nodes).

        :param node_ids: A dictionary which maps each BooleanFunction to a unique id.
        :type node_ids: dict[BooleanFunction, int]
        :return: A dictionary which represents the JSON for one node
        :rtype: dict[str, Any]
        """
        # copy class name and non-nested data
        JSON_object = {
            'class': type(self).__name__,
            'data': self.__dict__.copy()
        }

        # recurse on any children/nested data:
        if 'args' in JSON_object['data']:
            JSON_object['data']['args'] = [node_ids[arg] for arg in self.args]

        # ignore the compiled version (not serializable)
        if '_compiled' in JSON_object['data']:
            del JSON_object['data']['_compiled']

        return JSON_object
    
    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str,Any],
        parsed_functions: list["BooleanFunction | None"]
    ) -> Self:
        """Parse a JSON entry back into a BooleanFunction.

        This method is able to parse JSON generated by `_generate_JSON_entry` back into
        a BooleanFunction, and is used in `parse_JSON` to parse nodes back to their
        original type.

        :param object_data: A dictionary which contains the fields and data of the
            original node object, as generated by `_generate_JSON_entry`.
        :type object_data: dict[str, Any]
        :param parsed_functions: A list which contains the previously parsed functions.
            this can be used to make sure child functions are properly linked. Functions
            which have not yet been parsed will have None at their index instead.
        :type parsed_functions: list[BooleanFunction | None]
        :return: The parsed node, with data matching the JSON.
        :rtype: BooleanFunction
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

    @classmethod
    def generate_JSON(cls, 
        *fns: "BooleanFunction"
    ) -> dict[str, Any]:
        """Given a set of functions generate a JSON file which stores their information.

        This method is paired with `parse_JSON`. The parse method can reverse the JSON generated
        by this method back into an equivalent structure, and if you override this method, you
        should also change the parse method to correspond. By default, the generated JSON is a dict,
        structured as follows:

        - **"Return IDs"**: A list of integers, each specifying the ID of a function to return, 
            corresponding to the functions passed as input (usually the terminal nodes of the function DAG)
        - **"Node Data"**: A list, where each entry is a dict containing the information needed to reconstruct
            the node whose id matches that index in the list. The node data contains two fields,

            - **"class"**: The class name. This is used to identify and construct the appropriate
                BooleanFunction subclass and call the appropriate implementations of `_parse_JSON_entry`.
            - **"data"**: The additional data which is used to recreate the fields of the node 
                (in key:value pairs). For example:
                - **"args"**: a list containing the id's of the child nodes, which can be used to reconstruct the nodes arguments.
                - **"arg_limit"**: preserved as is from the node.
                - any additional fields (e.g. value/index for `CONST`/`VAR` nodes)

        One benefit of using this structure is that we can parse all nodes from their data,
        and then return the ones we want using the stored return IDs. When we have multiple
        functions which are enmeshed in a single DAG, this allows faithful reconstruction,
        which is not possible if we generate JSON for each function by itself. Note that the
        caller can control which functions are stored together, and which common nodes get copied
        by which functions they group into a single call (separate calls will never be part of
        a single DAG).

        For single functions, where this generality is not needed, there are some aliases
        which simplify the process and inputs. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manually manage file_IO 
        and using the full `generate_JSON`/`parse_JSON` methods.
        - `fn.to_JSON()`: an alias for `BooleanFunction.generate_JSON(fn)`
        - `BooleanFunction.from_JSON(json)`: an alias for `BooleanFunction.parse_JSON(json)[0]` 
        - `fn.to_file(filename)`: writes the output of `fn.to_JSON()` to a `.json` file
        - `BooleanFunction.from_file(filename)`: parses a `.json` file using `from_JSON()`

        :param fns: A variable number of booleanfunctions to store into the JSON file
        :type: BooleanFunction
        :return: The JSON object encoding the data of the input functions.
        :rtype: dict[str, Any]
        """
        node_ids = {}
        for fn in fns:
            node_ids = fn.generate_ids(node_ids)

        num_nodes = max(node_ids.values())+1
        json_node_list: list[Any] = [None for i in range(num_nodes)]

        for node,id in node_ids.items():
            json_node_list[id] = node._generate_JSON_entry(node_ids)
        
        return {
            "Return IDs": [node_ids[fn] for fn in fns],
            "Node Data": json_node_list
        }

    @classmethod
    def parse_JSON(cls,
        json_object: dict[str,Any]
    ) -> tuple["BooleanFunction"]:
        """Parses a JSON object back into the BooleanFunctions which generated it.

        This method is paired with `generate_JSON`. This method must be able to reverse the Generated
        JSON back into an equivalent structure, and if you override this method, you should ensure the
        generate ist still in correspondence. By default, the JSON is expected to have the following structure:
        
        - **"Return IDs"**: A list of integers, each specifying the ID of a function to return
            corresponding to the functions passed as input (usually the terminal nodes of the function DAG)
        - **"Node Data"**: A list, where each entry is a dict containing the information needed to reconstruct
            the node whose id matches that index in the list. The node data contains two fields,
            
            - **"class"**: The class name. This is used to identify and construct the appropriate
                BooleanFunction subclass and call the appropriate implementations of `_parse_JSON_entry`.
            - **"data"**: The additional data which is used to recreate the fields of the node 
                (in key:value pairs). For example:
                - **"args"**: a list containing the id's of the child nodes, which can be used to reconstruct the nodes arguments.
                - **"arg_limit"**: preserved as is from the node.
                - any additional fields (e.g. value/index for `CONST`/`VAR` nodes)

        This function will return a tuple of functions, in the same template as they were passed to
        `generate_functions`. Length checks may be needed if you are loading from a file that you don't
        know the structure of. For single functions, where this generality is not needed, there are some aliases
        which simplify the process and inputs. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manually manage file_IO 
        and using the full `generate_JSON`/`parse_JSON` methods.
        - `fn.to_JSON()`: an alias for `BooleanFunction.generate_JSON(fn)`
        - `BooleanFunction.from_JSON(json)`: an alias for `BooleanFunction.parse_JSON(json)[0]` 
        - `fn.to_file(filename)`: writes the output of `fn.to_JSON()` to a `.json` file
        - `BooleanFunction.from_file(filename)`: parses a `.json` file using `from_JSON()`

        :json_object: A dictionary with the expected structure.
        :type: dict[str,Any]
        :return: A tuple of booleanFunctions, parsed from the JSON.
        :rtype: tuple[BooleanFunction]
        """
        # parse object class and data
        return_ids = json_object["Return IDs"]
        json_node_list = json_object["Node Data"]
        num_nodes = len(json_node_list)
        parsed_functions: list[Any] = [None for i in range(num_nodes)]
        for node_id in range(num_nodes):
            node_data = json_node_list[node_id]
            
            # create information for the python object for this node
            object_class = None
            object_data = node_data['data']

            # find the appropriate subclass of BooleanFunction for the node
            for subcls in cls.__subclasses__():
                if subcls.__name__ == node_data['class']:
                    object_class = subcls

            # throw a better error if no class found
            if object_class == None:
                raise TypeError(f"Type \'{node_data['class']}\' is not a valid BooleanFunction")

            # put data into new object and add it to the parsed functions
            parsed_functions[node_id] = object_class._parse_JSON_entry(
                object_data, parsed_functions
            )
        
        # the root node is the last one in the list:
        return tuple([parsed_functions[node_id] for node_id in return_ids])

    def to_JSON(self) -> dict[str,Any]:
        """An alias for `BooleanFunction.generate_JSON(fn)`
        
        This can be used in conjunction with `from_JSON` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to use the full 
        `generate_JSON`/`parse_JSON` methods. Check the docstrings on these methods for more 
        information on usage and output.

        :return: A JSON object which encodes the input function
        :rtype: dict[str,Any]
        """
        return BooleanFunction.generate_JSON(self)

    @classmethod
    def from_JSON(cls, 
        json_object: dict[str,Any]
    ) -> "BooleanFunction":
        """An alias for `BooleanFunction.parse_JSON(json_object)[0]`
        
        This can be used in conjunction with `to_JSON` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to use the full 
        `generate_JSON`/`parse_JSON` methods. Check the docstrings on these methods for more 
        information on usage and output.

        :json_object: A dictionary with the expected structure.
        :type: dict[str,Any]
        :return: The BooleanFunction which was used to create the JSON.
        :rtype: BooleanFunction
        """
        return BooleanFunction.parse_JSON(json_object)[0]
    
    def to_file(self,
        filename: str
    ) -> None:
        """Writes the output of fn.to_JSON to a file with the given filename.
        
        This can be used in conjunction with `from_file` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manage I/O manually and use 
        the full `generate_JSON`/`parse_JSON` methods. Check the docstrings on these methods for more 
        information on usage and output.

        :param filename: A string which will be used as the name of the generated file 
            (a `.json` suffix is highly recommended)
        :type filename: str
        """
        with open(filename, 'w') as f:
            # also equivalent: f.write(json.dumps(self.to_JSON(), indent = 2))
            f.write(json.dumps(BooleanFunction.generate_JSON(self), indent = 2))

    @classmethod
    def from_file(cls, 
        filename: str
    ) -> "BooleanFunction":
        """Reads a single BooleanFunction from the file with the given filename.
        
        This can be used in conjunction with `to_file` to reduce verbosity and improve readibility
        when you only want to store/parse one function. These are useful shortcuts for a lot of cases,
        but once the use case becomes complex enough, its preferred to manage I/O manually and use 
        the full `generate_JSON`/`parse_JSON` methods. Check the docstrings on these methods for more 
        information on usage and output.

        :param filename: A string which givens the name of the file to read.
        :type filename: str
        """
        with open(filename, 'r') as f:
            # also equivalent: return BooleanFunction.from_json(json.loads(f.read()))
            return BooleanFunction.parse_JSON(json.loads(f.read()))[0]


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
        return max((arg.max_idx() for arg in self.args), default=-1)
    
    def idxs_used(self) -> set[int]:
        """Return the set of indices used in variables in the function.

        :return: Return the set of indices used in variables in the function.
        :rtype: set[int]
        """
        return set().union(*(arg.idxs_used() for arg in self.args))

    def num_nodes(self) -> int:
        """Return the number of nodes comprising the input function.

        :return: Return the number of nodes comprising the input function.
        :rtype: int
        """
        visited = set()
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in visited:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False or curr_node.is_leaf():
                visited.add(curr_node)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return len(visited)
   
    def component_count(self) -> dict[str,int]:
        """Return a dict which counts the occurrences of each class in the function DAG

        This is similar to `num_nodes`, but provides a more detailed description. The `__name__`
        attribute of each node to increment the corresponding count in the returned dictionary.

        :return: a dict which counts the occurrences of each class in the function DAG
        :rtype: dict[str,int]
        """
        components = {}
        visited = set()
        stack: list[Any] = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in visited:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False or curr_node.is_leaf():
                name = type(curr_node).__name__
                if name in components:
                    components[name] += 1
                else:
                    components[name] = 1
                    
                visited.add(curr_node)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return components
