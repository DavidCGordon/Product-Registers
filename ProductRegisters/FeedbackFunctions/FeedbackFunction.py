from ProductRegisters.BooleanLogic import BooleanFunction, VAR

# for compiling to c to iterate faster
import tempfile
import subprocess
from shutil import rmtree
import contextlib

# for compiling to python
import numpy as np
from numba import njit, jit
import types

# For Storing and loading as JSON files.
from copy import deepcopy
import json

class FeedbackFunction:
    def __init__(self, fn_list):
        #convert update to a list of ANF<int> objects:
        self.fn_list = fn_list
        self.size = len(fn_list)

    def __getitem__(self, idx): return self.fn_list[idx]

    def __setitem__(self, idx, val): self.fn_list[idx] = val

    def __len__(self): return self.size

    def __eq__(self, other): return self.fn_list == other.fn_list

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
            outstr += self.fn_list[i].pretty_str() + ";\n\n\n"
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

    #FUNCTION MANIPULATION:

    #new function that generates the same states with bit order reversed
    def flip(self):
        new_indices = {i: self.size-1-i for i in range(self.size)}
        self.fn_list = [f.remap_indices(new_indices) for f in self.fn_list][::-1]


    # return the number of gates before any optimization (VERY rough estimate of size)
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

    #return true if the function is linear
    def isLinear(self, allowAfine = False):
        for component in self.gateSummary().keys():
            if component not in ['XOR','CONST','VAR']:
                return False
            else:
                return True






    def to_JSON(self):
        # copy class name and non-nested data
        JSON_object = {
            'class': type(self).__name__,
            'data': self.__dict__.copy()
        }

        # convert fn_list/nested data:
        if 'fn_list' in JSON_object['data']:
            JSON_object['data']['fn_list'] = [f.to_JSON() for f in self.fn_list]

        # ignore the compiled version (not serializable)
        if '_compiled' in JSON_object['data']:
            del JSON_object['data']['_compiled']

        return JSON_object
    
    @classmethod
    def from_JSON(self, JSON_object):
        # parse object class and data
        object_data = JSON_object['data']
        object_class = None
        for subcls in self.__subclasses__():
            if subcls.__name__ == JSON_object['class']:
                object_class = subcls

        # throw a better error if no class found
        if object_class == None:
            raise TypeError(f"Type \'{JSON_object['class']}\' is not a valid FeedbackFunction")

        # put data into new object
        output = object.__new__(object_class)
        for key,value in object_data.items():
            if key == "fn_list":
                output.fn_list = [BooleanFunction.from_JSON(f) for f in value]
            else:
                setattr(output,key,value)
    
        return output

    # json files only:
    def to_file(self, filename):
        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    # json files only:
    @classmethod
    def from_file(self, filename):
        with open(filename, 'r') as f:
            return FeedbackFunction.from_JSON(json.loads(f.read()))







    @contextlib.contextmanager
    def compiled_to_c(self):
        self._data_store = tempfile.mkdtemp(prefix="ProductRegisters_")

        with open(self._data_store + "function_source.c","w") as f:
            f.write(f"""
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {{

    //parse number of cycles
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





    #writes a VHDL file
    #Credit: Anna Hemingway
    def write_VHDL(self, filename):
        with open(filename, "w") as f:
            f.write(f"""
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

    signal currstate, nextstate:std_logic_vector({self.size - 1} downto 0);


begin

    statereg: process(i_clk, i_rst)
    begin
        if (i_rst = '1') then
            currstate <= i_seed_data;
        elsif (i_clk = '1' and i_clk'event) then
            currstate <= nextstate;
        end if;
    end process;\n""")
        
            for i in range(self.size - 1, -1 , -1):
                f.write(f"    nextstate({str(i)}) <= {self.fn_list[i].generate_VHDL()};\n")
            f.write("""

    output <= currstate;

end run;

""")


    def __copy__(self):
        new_obj = object.__new__(type(self))
        new_obj.__dict__ = self.__dict__
        new_obj.fn_list = [f.__copy__() for f in self.fn_list]
        return new_obj

    # compile and run in python
    def compile(self):
        self._compiled = None

        exec_str = """
@njit(parallel=True)
def _compiled(currstate):
    nextstate = np.zeros_like(currstate)
"""
        for i in range(self.size - 1, -1 , -1):
                exec_str += f"    nextstate[{str(i)}] = {self.fn_list[i].generate_python()}\n"
        exec_str += "    return nextstate\n\n"

        exec_str += "self._compiled = _compiled"
        exec(exec_str)
        return self._compiled


    # # returns the nth, but not anything before it
    # def repeat(self, n, bits=0):
    #     # input handling
    #     if not bits:
    #         bits = list(range(self.size))

    #     # base case 
    #     if n == 0:
    #         res = [None for i in range(self.size)]
    #         for b in bits:
    #             res[b] = VAR(b)
    #         return res
        
    #     # recursive case
    #     prev = self.repeat(n-1, *bits)
    #     res = [None for i in range(self.size)]
    #     for bit in bits:
    #         res[bit] = prev[bit].compose(self.fn_list).to_ANF()
    #     return res



    # def iterator(self,n):
    #     # base case 
    #     if n == 0:
    #         yield [VAR(b) for b in range(self.size)]
            
    #     else:
    #         # yield all previous answers
    #         for prev in self.iterator(n-1):
    #             yield prev
            
    #         # build the newest answer using last layer
    #         yield [self.fn_list[b].compose(prev) for b in range(self.size)]
    def iterator(self, n):
        fns = [VAR(i) for i in range(self.size)]
        yield fns

        for i in range(1,n+1): 
            fns = [self.fn_list[b].compose(fns) for b in range(self.size)]
            yield fns

    # def anf_iterator(self,n,bits=None):
    #     """
    #     The composition order is reversed, because this leads to significantly smaller ANF operations
    #     because ANF is not deeply nested, we dont care about reference-reuse, so the order of the previous
    #     iterator makes less sense
    #     """

    #     # input handling
    #     if not bits:
    #         bits = list(range(self.size))

    #     # base case 
    #     if n == 0:
    #         res = [None for i in range(self.size)]
    #         for b in bits:
    #             res[b] = VAR(b)
    #         yield res

    #     else:
    #         # yield all previous answers
    #         for prev in self.anf_iterator(n-1,bits):
    #             yield prev
            
    #         # build the newest answer using last layer
    #         res = [None for i in range(self.size)]
    #         for bit in bits:
    #             res[bit] = prev[bit].compose(self.fn_list).to_ANF()
    #         yield res

    def anf_iterator(self,n,bits=None):
        #input handling:
        if not bits:
             bits = list(range(self.size))

        fns = [VAR(b) if b in bits else None for b in range(self.size)]
        yield fns

        for i in range(1,n+1):
            
            fns = [fns[b].compose(self.fn_list).translate_ANF() for b in bits]
            yield fns