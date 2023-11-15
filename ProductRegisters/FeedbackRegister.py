from ProductRegisters.FeedbackFunctions import FeedbackFunction


import numpy as np

import subprocess

from numba import jit

import json
#import system


class FeedbackRegister:
    #INITIALIATION/DATA:
    def __init__(self, seed, fn):
        
        # attributes
        self.fn = fn
        self.size = len(fn)

        # set seed
        self.seed(seed)
        # set state to seed
        self._state = self._seed.copy()

        # set next state appropriately
        self._prev_state = np.zeros_like(self._seed)

    def __len__(self): return self.size


    #REGISTER SEED:
    def seed(self, seed):
        # For a given seed
        if type(seed) == int:
            self._seed = np.asarray([int(x) for x in format(seed, f'0{self.size}b')[::-1]], dtype='uint8')
        if type(seed) == list:
            self._seed = np.asarray(seed, dtype='uint8')

        # Allowing random.random to be used as a seed:
        if type(seed) == float and 0 < seed  and seed < 1:
            closest_int = round(seed * 2**self.size)
            self.seed(closest_int)
        
    def reset(self):
        self._state = self._seed.copy()

    #TYPE CONVERSIONS / CASTING:
    def __str__(self): return "".join(str(x) for x in self._state[::-1])
    def __int__(self): return int(str(self),2)
    def __list__(self): return self._state[:]

    #STATE MANIPULATION:
    def __getitem__(self, key): return self._state[key].copy()
    def __setitem__(self, key, val): self._state[key] = val
    
    #ITERATION THROUGH REGISTER BITS:
    def __iter__(self): return iter(self._state)
    def __reversed__(self): return reversed(self._state)


    #CLOCKING AND RUNNING THE REGISTER:
    def clock(self):
        # swap pointers to move _state to _prev_state
        self._state, self._prev_state = self._prev_state, self._state

        # overwrite the new nextstate
        for i in range(len(self.fn.fn_list)):
            self._state[i] = self.fn.fn_list[i].eval(self._prev_state)

    def clock_compiled(self):
        if not hasattr(self.fn, "_compiled"):
            print("Compile first!")
            return
        self._state, self._prev_state = self._prev_state, self._state
        self._state = self.fn._compiled(self._prev_state)

    #generate a sequence of states, in order
    def run(self, limit = None):
        #number of iterations to run
        if type(limit) == int:
            for _ in range(limit):
                yield self
                self.clock()
                
        #no limit
        elif limit == None:
            while True:
                yield self
                self.clock()

    def run_compiled(self, arg = None):
        if not hasattr(self.fn, "_compiled"):
            print("Compile first!")
            return

        update_fn = self.fn._compiled
        # number of iterations to run
        if type(arg) == int:
            for _ in range(arg):
                yield self
                # swap pointers to move _state to _prev_state
                self._state, self._prev_state = self._prev_state, self._state
                # overwrite _state with the new state
                self._state = update_fn(self._prev_state)
                
        #no limit
        elif arg == None:
            while True:
                yield self
                # swap pointers to move next_state to curr_state
                self._state, self._prev_state = self._prev_state, self._state
                # overwrite the new nextstate
                self._state = update_fn(self._prev_state)


                

    #DIAGNOSTIC AND EXTRA INFO     
    #return the period of the register:
    def period(self, lim = 2**18):
        first_state = self._state.copy()
        self.clock()
        count = 1
        while not(all(self._state == first_state)):
            self.clock()
            count += 1
            if count > lim:
                self._state = first_state.copy()
                return None

        self._state = first_state.copy()
        return count

    def period_compiled(self, iter_lim = 2**22, bit_lim = None):
        if not hasattr(self.fn, "_compiled"):
            print("Compile first!")
            return
        
        first_state = self._state.copy()
        self.clock()
        count = 1

        for state in self.run_compiled(iter_lim):
            if all(state._state[bit_lim:] == first_state[bit_lim:]):
                break
            count += 1

        self._state = first_state.copy()
        
        if count > iter_lim:
            return None
        else:
            return count


    # must have finite limit
    # WORSE THAN RUN_COMPILED()
    def runC(self,lim):
        if not hasattr(self.fn, "_data_store"):
            print("Compile first!")
            return 

        process = subprocess.Popen(
            [self.fn._data_store + "function_iteration.exe", f"{lim}", "".join(str(i) for i in self)],
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW
        )

        try:
            while True:
                line = process.stdout.readline()
                if not line:
                    break
                current = [int(x) for x in line.decode('utf8').strip()]
                #current = [int(x) for x in reversed(line.decode('utf8').strip())]
                yield(current)

        finally:
            process.kill()
            self._state = np.asarray(current, dtype='uint8')
            self.clock()
    


    def to_JSON(self):
        # copy class name and non-nested data
        JSON_object = {
            'class': type(self).__name__,
            'data': self.__dict__.copy()
        }

        # add data:
        JSON_object['data']['fn'] = self.fn.to_JSON()
        JSON_object['data']['_seed'] = self._seed.tolist()
        JSON_object['data']['_state'] = self._state.tolist()
        JSON_object['data']['_prev_state'] = self._prev_state.tolist()
        return JSON_object
    
    @classmethod
    def from_JSON(self, JSON_object):
        # parse object class and data
        object_data = JSON_object['data']
        if JSON_object['class'] != 'FeedbackRegister':
            raise TypeError(f"Expected type \'FeedbackRegister\', but got \'{JSON_object['class']}\'")

        # put data into new object
        output = object.__new__(self)
        for key,value in object_data.items():
            if key == "fn":
                output.fn = FeedbackFunction.from_JSON(value)
            elif key == "_seed":
                output._seed = np.asarray(value,dtype='uint8')
            elif key == "_state":
                output._state = np.asarray(value,dtype='uint8')
            elif key == "_prev_state":
                output._prev_state = np.asarray(value,dtype='uint8')
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
            return FeedbackRegister.from_JSON(json.loads(f.read()))