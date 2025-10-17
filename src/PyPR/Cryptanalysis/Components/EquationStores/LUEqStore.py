from PyPR.BooleanLogic import BooleanFunction
from PyPR.BooleanLogic import CONST

import numpy as np
from numpy.typing import NDArray


import numba
u8 = numba.types.uint8
u64 = numba.types.uint64
b1 = numba.types.b1

from typing import Any

# general utilities:
@numba.njit(numba.types.Tuple((b1,u64))(
        u8[:,:],u8[:,:],u8[:],u8[:],u64,
        u8[:],u8,b1)
)
def LU_reduction(
    upper,lower,constants,solved_for,num_vars,
    coef_vector,const_val,consistent
):
    linearly_independent = False
    coef_vector = coef_vector.copy()
    modification_vector = np.zeros_like(coef_vector)
    for idx in range(len(coef_vector)):
        if coef_vector[idx] == 1:
            modification_vector[idx] = 1
            if solved_for[idx]:
                # this changes constants to actually be L^{-1} C
                const_val ^= constants[idx]
                for i in range(idx,num_vars):
                    coef_vector[i] ^= upper[idx,i]
                    
            else:
                linearly_independent = True

                solved_for[idx] = 1
                upper[idx] = coef_vector
                lower[idx] = modification_vector        
                constants[idx] = const_val
                break

    if (
        consistent
        and not linearly_independent
        and const_val == 1
    ):
        raise ValueError("Inconsistent!")

    return linearly_independent, idx


class LUEqStore:
    def __init__(
        self, 
        comb_to_idx: dict[tuple[int], int],
        consistent: bool = False
    ):        
        self.consistent = consistent
        self.comb_to_idx = {k:v for k,v in comb_to_idx.items()}
        self.idx_to_comb = {v:k for k,v in comb_to_idx.items()}
        self.equation_ids = {} # mapping of equation to identifier (usually a clock cycle)

        self.num_vars = len(comb_to_idx)
        self.num_eqs = 0

        self.lower_matrix = np.eye(self.num_vars, dtype = 'uint8')
        self.upper_matrix = np.eye(self.num_vars, dtype = 'uint8')
        self.constants = np.zeros([self.num_vars], dtype = 'uint8')
        self.solved_for = np.zeros([self.num_vars], dtype = 'uint8')


    def insert_equation(
        self, 
        equation: BooleanFunction | list[int] | NDArray, 
        extra_const: int = 0, 
        identifier: Any = None, 
        translate_ANF: bool = True
    ):
        const_val = extra_const

        # shortcut to allow accepting an array directly
        if type(equation) == np.ndarray and equation.shape == (self.num_vars,):
            coef_vector = equation
            
        elif type(equation) == BooleanFunction:
            if translate_ANF:
                equation_anf = equation.translate_ANF()
            else:
                equation_anf = equation

            coef_vector = np.zeros([self.num_vars], dtype=np.uint8)
            for term in equation_anf.args:

                # handle constants
                if type(term) == CONST:
                    const_val ^= term.value
                    continue

                # do not need to expand for new variables:
                comb = tuple(sorted([var.index for var in term.args])) #type: ignore
                coef_vector[self.comb_to_idx[comb]] = 1

        linearly_independent, insertion_idx =  LU_reduction(
            self.upper_matrix, self.lower_matrix, self.constants, self.solved_for, self.num_vars,
            coef_vector, const_val, self.consistent
        )

        # update fields which can't be modified in numba function
        if linearly_independent:
            self.equation_ids[insertion_idx] = identifier
            self.num_eqs += 1

        return linearly_independent

    @property
    def rank(self):
        return self.num_eqs
