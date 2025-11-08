from typing import Any, Self

from PyPR.BooleanLogic import BooleanFunction, VAR
from PyPR.BooleanLogic.BooleanANF import BooleanANF
import PyPR.JSON_Serialization

# for compiling to c to iterate faster
import tempfile
import subprocess
from shutil import rmtree
import contextlib

# for compiling to python
import numpy as np
from numba import njit

# For Storing and loading as JSON files.
import json

class FeedbackFunction:
    fn_list: list[BooleanFunction]
    size: int
    _compiled: Any
    _compiled_inplace: Any

    def __init__(self, fn_list):
        #convert update to a list of ANF<int> objects:
        self.fn_list = fn_list
        self.size = len(fn_list)
    
    def __copy__(self):
        new_obj = object.__new__(type(self))
        new_obj.__dict__ = self.__dict__
        new_obj.fn_list = [f.__copy__() for f in self.fn_list]
        return new_obj

    #TODO: expand on this
    def copy(self):
        return self.__copy__()
    
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
        new_indices = {i: self.size-1-i for i in range(self.size)}
        self.fn_list = [f.remap_indices(new_indices) for f in self.fn_list][::-1]


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
        elif in_place: 
            ids = previous_ids
        else:
            # shallow copy to maintain objects, but new id dict
            ids = {k:v for k,v in previous_ids.items()}

        # create ids for all functions:
        for fn in self.fn_list:
            ids = fn.generate_ids(ids)
        
        # lastly, add id for self:
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
        
    def to_JSON(self):
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
        :return: The FeedbackFunction which was used to create the JSON.
        :rtype: FeedbackFunction
        """
        return_idx = json_object['return order'][0]
        json_class = json_object['objects'][return_idx]['class']
        subclasses = set((
            str(cls)[8:-2] for cls in 
            PyPR.JSON_Serialization.all_subclasses(cls)
        ))

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


    # text generation
    # TODO: rename to match BF naming
    def write_VHDL(self, filename):
        #writes a VHDL file
        #Credit: Anna Hemingway
        overrides = {}
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


begin

    statereg: process(i_clk, i_rst)
    begin
        if (i_rst = '1') then
            curr_state <= i_seed_data;
        elsif (i_clk = '1' and i_clk'event) then
            curr_state <= next_state;
        end if;
    end process;\n"""
        
        vhdl_str += "\n    "
        for i in range(self.size - 1, -1 , -1):
            vhdl_str += ("\n    ".join(self.fn_list[i].generate_VHDL(
                output_name = f"next_state({i})",
                array_name = "curr_state",
                subfunction_prefix = f"fn_{i}",
                overrides = overrides
            )) + "\n    ")

            for j, node in enumerate(self.fn_list[i].subfunctions()):
                if node not in overrides:
                    overrides[node] = f'fn_{i}_{j+1}'

        vhdl_str += """
    output <= currstate;

end run;

"""
        with open(filename, "w") as f:
            f.write(vhdl_str)

    # def write_tex(self, filename):
    #     with open(filename, "w") as f:
    #         for i in range(self.size - 1, -1 , -1):
    #             f.write(f"c_{{{str(i)}}}[t+1] &= {self.fn_list[i].generate_tex()}\\\\\n")

    # Compilation
    def compile(self):
        self._compiled = None
        self._compiled_inplace = None

        # return a new answer
        overrides = {}
        exec_str = """
@njit(parallel=True)
def _compiled(curr_state):
    next_state = np.zeros_like(curr_state)
"""
        exec_str += "\n    "
        for i in range(self.size - 1, -1 , -1):
            exec_str += ("\n    ".join(self.fn_list[i].generate_python(
                output_name = f"next_state[{i}]",
                array_name = "curr_state",
                subfunction_prefix = f"fn_{i}",
                overrides = overrides
            )) + "\n    ")

            for j, node in enumerate(self.fn_list[i].subfunctions()):
                if node not in overrides:
                    overrides[node] = f'fn_{i}_{j+1}'
            
        exec_str += "return next_state\n\n"
        exec_str += "self._compiled = _compiled"
        exec(exec_str)

        # write to an existing buffer
        overrides = {}
        exec_str = """
@njit(parallel=True)
def _compiled_inplace(curr_state,output_buffer):
"""
        exec_str += ("    ")
        for i in range(self.size - 1, -1 , -1):
            exec_str += ("\n    ".join(self.fn_list[i].generate_python(
                output_name = f"output_buffer[{i}]",
                array_name = "curr_state",
                subfunction_prefix = f"fn_{i}",
                overrides = overrides
            )) + "\n    ")

            for j, node in enumerate(self.fn_list[i].subfunctions()):
                if node not in overrides:
                    overrides[node] = f'fn_{i}_{j+1}'
        exec_str += "return\n\n"
        exec_str += "self._compiled_inplace = _compiled_inplace"
        exec(exec_str)

        return self._compiled

    # Function unrolling (possibly remove)
    def iterator(self, n):
        fns = [VAR(i) for i in range(self.size)]
        yield fns

        for i in range(1,n+1): 
            fns = [self.fn_list[b].compose(fns) for b in range(self.size)]
            yield fns


    # Statistics:
    def gateSummary(self):
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

    def isLinear(self, allowAfine = False):
        for component in self.gateSummary().keys():
            if component not in ['XOR','CONST','VAR']:
                return False
            else:
                return True