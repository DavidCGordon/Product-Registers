from typing import Any

from PyPR.BooleanLogic import BooleanFunction, VAR
from PyPR.BooleanLogic.BooleanANF import BooleanANF

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
    def __init__(self, fn_list):
        #convert update to a list of ANF<int> objects:
        self.fn_list = fn_list
        self.size = len(fn_list)
    
    def __copy__(self):
        new_obj = object.__new__(type(self))
        new_obj.__dict__ = self.__dict__
        new_obj.fn_list = [f.__copy__() for f in self.fn_list]
        return new_obj

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


    # Storage:
    def to_JSON(self):
        # copy class name and non-nested data
        JSON_object = {
            'class': type(self).__name__,
            'data': self.__dict__.copy()
        }

        # convert fn_list/nested data:
        if 'fn_list' in JSON_object['data']:
            # JSON_object['data']['fn_list'] = [f.to_JSON() for f in self.fn_list]
            JSON_object['data']['fn_list'] = BooleanFunction.generate_JSON(*self.fn_list)

        # ignore the compiled version (not serializable)
        if '_compiled' in JSON_object['data']:
            del JSON_object['data']['_compiled']

        return JSON_object
    
    @classmethod
    def from_JSON(cls, JSON_object):
        # parse object class and data
        object_data = JSON_object['data']
        object_class = None
        for subcls in cls.__subclasses__():
            if subcls.__name__ == JSON_object['class']:
                object_class = subcls

        # throw a better error if no class found
        if object_class == None:
            raise TypeError(f"Type \'{JSON_object['class']}\' is not a valid FeedbackFunction")

        # put data into new object
        output = object.__new__(object_class)
        for key,value in object_data.items():
            if key == "fn_list":
                #output.fn_list = [BooleanFunction.from_JSON(f) for f in value]
                output.fn_list = list(BooleanFunction.parse_JSON(value))
            else:
                setattr(output,key,value)
    
        return output
    
    def to_file(self, filename):
        # json files only:
        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    @classmethod
    def from_file(cls, filename):
        # json files only:
        with open(filename, 'r') as f:
            return FeedbackFunction.from_JSON(json.loads(f.read()))


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

    def write_tex(self, filename):
        with open(filename, "w") as f:
            for i in range(self.size - 1, -1 , -1):
                f.write(f"c_{{{str(i)}}}[t+1] &= {self.fn_list[i].generate_tex()}\\\\\n")

    # Compilation
    @contextlib.contextmanager
    def compiled_to_c(self):
        self._data_store = tempfile.mkdtemp(prefix="ProductRegisters_")

        with open(self._data_store + "function_source.c","w") as f:
            f.write(f"""
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {{

    //parse number of cycles from cmd
    unsigned long long limit = strtoll(argv[1],NULL,10);


    //parse and init arrays:
    unsigned short arr1[{self.size}];
    unsigned short arr2[{self.size}];

    for(int i = 0; argv[2][i] != 0; i++) {{
        arr1[i] = (unsigned short)(argv[2][i] - '0');
    }}
    
    unsigned short (*currstate)[{self.size}] = &arr1;
    unsigned short (*nextstate)[{self.size}] = &arr2;
    unsigned short (*temporary)[{self.size}];

    for (int i = 0; i<limit; i++) {{\n""")

            for i in range(self.size - 1, -1 , -1):
                f.write(f"        (*nextstate)[{str(i)}] = {self.fn_list[i].generate_c()};\n")

            f.write(f"""
        for (int j = 0; j < {self.size}; j++) {{
            printf("%hu", (*currstate)[j]);
        }}
        printf("\\n");

        temporary = currstate;
        currstate = nextstate;
        nextstate = temporary;
    }}
    return 0;
}}""")

        subprocess.run(
            ["gcc", self._data_store +"function_source.c", "-o", self._data_store + "function_iteration.exe"],
            creationflags=subprocess.CREATE_NO_WINDOW
        )

        yield
        print(f"deleting {self._data_store}")
        rmtree(self._data_store)
        del self._data_store

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

    # Probably remove
    def anf_iterator_1(
        self,
        rounds,
        bits=None,
        initialization = None,
    ):
        if bits == None: 
            bits = list(range(self.size))

        # if initialization == None:
            
        # else:
        #     fns = initialization

        fns: list[Any] = [VAR(b) if b in bits else None for b in range(self.size)]
        out: list[Any] = [None]*self.size
        for b in bits:
            out[b] = fns[b].compose(initialization).translate_ANF()

        yield out
        for i in range(1,rounds+1):
            #fns = [self.fn_list[b].compose(fns).translate_ANF() for b in range(self.size)]
            
            for b in bits:
                fns[b] = fns[b].compose(self.fn_list).translate_ANF()
                out[b] = fns[b].compose(initialization).translate_ANF()
            # fns = [
            #     fns[b].compose(self.fn_list).translate_ANF() if b in bits else None
            #     for b in range(self.size)
            # ]
            
            yield out
            
    def anf_iterator_2(
        self,
        rounds,
        initialization = None,
    ):
        #optimized_list = [self.fn_list[b].anf_optimize() for b in range(self.size)]
        #optimized_list = self.fn_list

        eval_list = [
            f.anf_optimize().remap_constants([
                (0, BooleanANF()),
                (1, BooleanANF([True]))
            ]) for f in self.fn_list
        ]

        if initialization == None:
            fns = [BooleanANF([b]) for b in range(self.size)]
        else:
            fns = [BooleanANF.from_BooleanFunction(f) for f in initialization]

        yield [f.to_BooleanFunction() for f in fns]
        for i in range(1,rounds+1):
            fns = [eval_list[b].eval_ANF(fns) for b in range(self.size)]
            yield [f.to_BooleanFunction() for f in fns]
    
        # return the number of gates before any optimization (VERY rough estimate of size)
    
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