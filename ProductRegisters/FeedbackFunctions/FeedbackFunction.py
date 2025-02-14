from ProductRegisters.BooleanLogic import BooleanFunction, VAR

# for compiling to c to iterate faster
import tempfile
import subprocess
from shutil import rmtree
import contextlib

# for compiling to python
import numpy as np
from numba import njit
import types

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

    # def __eq__(self, other): return self.fn_list == other.fn_list

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
    
    def to_file(self, filename):
        # json files only:
        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    @classmethod
    def from_file(self, filename):
        # json files only:
        with open(filename, 'r') as f:
            return FeedbackFunction.from_JSON(json.loads(f.read()))

    #writes a VHDL file
    #Credit: Anna Hemingway
    def write_VHDL(self, filename):
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

    # generate C implementation for a CMPR
    # by default, the C implementation uses an initial state of all 1's and clocks for 25 cycles (including the initial state)
    def write_C(self, filename, limit=100):
        with open(filename, "w") as f:
            f.write(f"""
#include <stdio.h>
#include <stdlib.h>
int main() {{
    unsigned long long limit = {limit};
    unsigned short arr1[{self.size}]; // initial state: index 0 is the LSB
    for (int i = 0; i < {self.size}; i++) {{
        arr1[i] = 1;
    }}
    unsigned short arr2[{self.size}]; // next state: always initialize to all 0's
    for (int i = 0; i < {self.size}; i++) {{
        arr2[i] = 0;
    }}
    
    unsigned short (*currstate)[{self.size}] = &arr1;
    unsigned short (*nextstate)[{self.size}] = &arr2;
    unsigned short (*temporary)[{self.size}];
    for (unsigned long long cycle = 0; cycle < limit; cycle++) {{\n""")

            for i in range(self.size - 1, -1 , -1):
                f.write(f"        (*nextstate)[{str(i)}] = {self.fn_list[i].generate_c()};\n")

            f.write(f"""
        for (int j = 0; j < {self.size}; j++) {{
            // currstate[0] is the LSB of the state
            printf("%hu", (*currstate)[j]);
        }}
        printf("\\n");
        temporary = currstate;
        currstate = nextstate;
        nextstate = temporary;
    }}
    return 0;
}}""")

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



    def iterator(self, n):
        fns = [VAR(i) for i in range(self.size)]
        yield fns

        for i in range(1,n+1): 
            fns = [self.fn_list[b].compose(fns) for b in range(self.size)]
            yield fns

    def anf_iterator(self,n):
        fns = [VAR(b) for b in range(self.size)]
        yield fns

        for i in range(1,n+1):
            fns = [self.fn_list[b].compose(fns).translate_ANF() for b in range(self.size)]
            yield fns
            
            