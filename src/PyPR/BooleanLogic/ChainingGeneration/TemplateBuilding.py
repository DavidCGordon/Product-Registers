from typing import Any

from PyPR.BooleanLogic import BooleanFunction
from PyPR.BooleanLogic.FunctionInputs import VAR,CONST
from itertools import product

import numpy as np
import random

class TEMPLATE:
    def __init__(self):
        pass

    def _sample(self, allow_empty_return=False):
        raise NotImplementedError
    def _clean_up(self):
        raise NotImplementedError
    
    def sample(self):
        output = self._sample()
        self._clean_up()    
        return output


class VALUE(TEMPLATE):
    def __init__(self, value):
        self.value = value
    
    def _sample(self, allow_empty_return = False):
        return [self.value]
    
    def _clean_up(self):
        pass

                 

class GATE(TEMPLATE):
    def __init__(self,
        parameters = {},
        sources = []
    ):
        if "gate_class" in parameters:
            self.gate_class = parameters["gate_class"]
        else:
            raise ValueError('GATE template object requires the parameter "gate_class", but this was not supplied.')
        
        self.sources = sources
        if not self.sources:
            raise ValueError("No Sources Defined")

    def _sample(self, allow_empty_return = False):
        outputs = []
        for source in self.sources:
            outputs += source._sample()

        return [self.gate_class(*outputs)]
    
    def _clean_up(self):
        for s in self.sources:
            s._clean_up()


class FUNCTION(TEMPLATE):
    def __init__(self,
        parameters = {},
        sources = []
    ):
        if "fn" in parameters:
            self.fn = parameters["fn"]
        else:
            raise ValueError('FUNCTION template object requires the parameter "fn", but this was not supplied.')
        
        self.sources = sources
        if not self.sources:
            raise ValueError("No Sources Defined")
        
    def _sample(self, allow_empty_return = False):
        outputs = []
        for source in self.sources:
            outputs += source._sample()
            
        return [self.fn.compose(outputs)]

    def _clean_up(self):
        for s in self.sources:
            s._clean_up()

        

class SAMPLE(TEMPLATE): # pick one child according to distribution
    def __init__(self,
        parameters = {
            "num_sampled": 1
        },
        source_distribution = []
    ):
        self.num_sampled = parameters.get("num_sampled",1)
        self.source_distribution = source_distribution

        if not self.source_distribution:
            raise ValueError("No Sources Defined")

    def _sample(self, allow_empty_return = False):
        probs = np.array([x[0] for x in self.source_distribution])
        sources = np.array([x[1] for x in self.source_distribution])

        selected_sources = np.random.choice(
            a=sources,
            size=self.num_sampled,
            replace=False,
            p=probs
        )

        outputs = []
        for source in selected_sources:
            for output in source._sample():
                outputs.append(output)

        return outputs
    
    def _clean_up(self):
        for p,source in self.source_distribution:
            source._clean_up()



# Filters
class NONCONSTANT(TEMPLATE):
    source: TEMPLATE
    def __init__(self,
        parameters = {
            "attempt_limit": 1000,
            "disable_on_failure": False
        },
        source: Any = None
    ):
        self.attempt_limit = parameters.get("attempt_limit",1000)
        self.disable_on_failure = parameters.get("disable_on_failure", False)
        self.source = source
        if self.source == None:
            raise ValueError("No Source Defined")

    def _sample(self, allow_empty_return = False):
        
        for i in range(self.attempt_limit):
            output = self.source._sample()
            if len(output) == 0 and not allow_empty_return:
                raise ValueError("NONCONSTANT template object expects an output from its source, but none was returned.")
            if len(output) > 1:
                raise ValueError("NONCONSTANT template object expects a single output from its source, but several were returned.")

            is_constant = False
            for constant in [CONST(0),CONST(1)]:
                if output[0].functionally_equivalent(constant):
                    is_constant = True
                    break

            if not is_constant:
                return output
            if allow_empty_return:
                return []
            
        if self.disable_on_failure:
            return output
        else:
            raise ValueError(f"NONCONSTANT template object was unable to find an output after {self.attempt_limit} attempts.")
    
    def _clean_up(self):    
        self.source._clean_up()



class DISTINCT(TEMPLATE):
    _function_cache = {}
    source: TEMPLATE
    def __init__(self,
        parameters = {
            "group_id": 0,
            "attempt_limit": 1000,
            "disable_on_failure": False
        },
        source: Any = None
    ):
        self.group_id = parameters.get("group_id", 0)
        if not (self.group_id in DISTINCT._function_cache):
            DISTINCT._function_cache[self.group_id] = []

        self.attempt_limit = parameters.get("attempt_limit",1000)
        self.disable_on_failure = parameters.get("disable_on_failure", False)
        self.source = source
        if self.source == None:
            raise ValueError("No Source Defined")
        

    def _sample(self, allow_empty_return = False):
        for i in range(self.attempt_limit):
            output = self.source._sample()
            if len(output) == 0 and not allow_empty_return:
                raise ValueError("DISTINCT template object expects an output from its source, but none was returned.")
            if len(output) > 1:
                raise ValueError("DISTINCT template object expects a single output from its source, but several were returned.")

            is_distinct = True
            for fn in DISTINCT._function_cache[self.group_id]:
                if output[0].functionally_equivalent(fn):
                    is_distinct = False
                    break

            if is_distinct:
                DISTINCT._function_cache[self.group_id].append(output[0])
                return output
            if allow_empty_return:
                return []
            
        if self.disable_on_failure:
            return output
        else:
            raise ValueError(f"DISTINCT template object was unable to find an output after {self.attempt_limit} attempts.")
    

    def _clean_up(self):
        for group in DISTINCT._function_cache:
            if DISTINCT._function_cache[group] != []: 
                DISTINCT._function_cache[group]  = [] 
        self.source._clean_up()

# the difference between unique and distinct is that:
# distinct clears cache after each call,
# while unique does not.
class UNIQUE(TEMPLATE):
    _function_cache = {}
    source: TEMPLATE

    def __init__(self,
        parameters = {
            "group_id": 0,
            "attempt_limit": 1000,
            "disable_on_failure": False
        },
        source: Any = None
    ):
        self.group_id = parameters.get("group_id", 0)
        if not (self.group_id in UNIQUE._function_cache):
            UNIQUE._function_cache[self.group_id] = []

        self.attempt_limit = parameters.get("attempt_limit",1000)
        self.disable_on_failure = parameters.get("disable_on_failure", False)

        self.source = source
        if self.source == None:
            raise ValueError("No Source Defined")
        

    def _sample(self, allow_empty_return = False):
        for i in range(self.attempt_limit):
            output = self.source._sample()
            if len(output) == 0 and not allow_empty_return:
                raise ValueError("UNIQUE template object expects an output from its source, but none was returned.")
            if len(output) > 1:
                raise ValueError("UNIQUE template object expects a single output from its source, but several were returned.")

            is_distinct = True
            for fn in UNIQUE._function_cache[self.group_id]:
                if output[0].functionally_equivalent(fn):
                    is_distinct = False
                    break

            if is_distinct:
                UNIQUE._function_cache[self.group_id].append(output[0])
                return output
            if allow_empty_return:
                return []

        if self.disable_on_failure:
            return output
        else:  
            raise ValueError(f"UNIQUE template object was unable to find an output after {self.attempt_limit} attempts.")
    
    def _clean_up(self):    
        self.source._clean_up()

class FIXED(TEMPLATE):
    source: TEMPLATE

    def __init__(self,
        parameters = {},
        source: Any = None
    ):
        self.returned = None
        self.source = source
        if self.source == None:
            raise ValueError("No Source Defined")

    def _sample(self, allow_empty_return=False):
        if self.returned:
            return self.returned

        else:
            output = self.source._sample()
            self.returned = output
            return output
        
    def _clean_up(self):    
        self.source._clean_up()

class REPEAT(TEMPLATE):
    source: TEMPLATE
    def __init__(self,
        parameters = {},
        source: Any = None
    ):
        if "iterations" in parameters:
            self.iterations = parameters["iterations"]
        else:
            raise ValueError('REPEAT template object requires the parameter "iterations", but this was not supplied.')
        
        self.source = source
        if self.source == None:
            raise ValueError("No Source Defined")


    def _sample(self, allow_empty_return = False):
        outputs = []
        for i in range(self.iterations):
            outputs += self.source._sample()
        
        return outputs
    
    def _clean_up(self):    
        self.source._clean_up()
    

class OPTIONAL(TEMPLATE):
    source: TEMPLATE

    def __init__(self,
        parameters = {"drop_chance": 0},
        source: Any = None
    ):
        self.drop_chance = parameters.get("drop_chance",0)
        self.source = source
        if self.source == None:
            raise ValueError("No Source Defined")

    def _sample(self):
        if random.random() < self.drop_chance:
            return []
        else:
            return self.source._sample(allow_empty_return = True)
    
    def _clean_up(self):
        self.source._clean_up()
