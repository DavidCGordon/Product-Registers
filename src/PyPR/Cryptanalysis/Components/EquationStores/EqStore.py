from typing import Any
from PyPR.BooleanLogic import BooleanFunction, CONST

import numpy as np

class EqStore:
    def __init__(self, comb_to_idx, integrated_constants=False):        
        self.comb_to_idx = {**comb_to_idx}
        self.idx_to_comb = {v:k for k,v in comb_to_idx.items()}
        self.equation_ids = {} # mapping of equation to identifier (usually a clock cycle)
        if integrated_constants and tuple() not in comb_to_idx:
            self.comb_to_idx[tuple()] = len(comb_to_idx)
            self.idx_to_comb[len(comb_to_idx)] = tuple()

        self.num_vars = len(comb_to_idx)
        self.num_eqs = 0

        self.equations = np.zeros([256,self.num_vars], dtype = 'uint8')
        self.constants = np.zeros([256], dtype = 'uint8')


    def insert_equation(self, 
        equation: BooleanFunction | np.ndarray, 
        extra_const: int = 0, 
        identifier: Any = None, 
        translate_ANF = True
    ):
        const_val = extra_const

        # shortcut to allow accepting an array directly
        if type(equation) == np.ndarray and equation.shape == (self.num_vars,):
            coef_vector = equation
            
        elif isinstance(equation, BooleanFunction):
            if translate_ANF:
                equation_anf = equation.translate_ANF()
            else:
                equation_anf = equation

            coef_vector = np.zeros([self.num_vars], dtype=np.uint8)
            const_val = extra_const
            for term in equation_anf.args:

                # handle constant values
                if type(term) == CONST:
                    if tuple() in self.comb_to_idx:
                        coef_vector[self.comb_to_idx[tuple()]] = term.value
                    else:
                        const_val ^= term.value
                    continue

                # do not need to expand for new variables:
                comb = tuple(sorted([var.index for var in term.args])) #type: ignore
                coef_vector[self.comb_to_idx[comb]] = 1


        # expand number of equations as necessary:
        if self.num_eqs == self.equations.shape[0]:
            # initialize
            new_equations = np.zeros([self.num_eqs* 2, self.num_vars],  dtype=np.uint8)
            new_constants = np.zeros([self.num_eqs* 2],  dtype=np.uint8)

            # copy
            new_equations[:self.num_eqs] = self.equations
            new_constants[:self.num_eqs] = self.constants

            # replace
            self.equations = new_equations
            self.constants = new_constants
        
        self.equation_ids[self.num_eqs] = identifier
        self.equations[self.num_eqs] = coef_vector
        self.constants[self.num_eqs] = const_val
        self.num_eqs += 1
        return True