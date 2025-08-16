from PyPR.BooleanLogic import CONST

import numpy as np
import numba
u8 = numba.types.uint8
u64 = numba.types.uint64
b1 = numba.types.b1

# general utilities:
@numba.njit(numba.types.Tuple((b1,u64))(
        u8[:,:],u8[:,:],u8[:],u8[:],u64,
        u8[:],u8)
)
def LU_reduction(
    upper,lower,constants,solved_for,num_vars,
    coef_vector,const_val
):
    linearly_independent = False
    modification_vector  = np.zeros_like(coef_vector)
    for idx in range(len(coef_vector)):
        if coef_vector[idx] == 1:
            modification_vector[idx] = 1
            if solved_for[idx]:
                for i in range(idx,num_vars):
                    coef_vector[i] ^= upper[idx,i]
            else:
                linearly_independent = True

                solved_for[idx] = 1
                upper[idx] = coef_vector
                lower[idx] = modification_vector        
                constants[idx] = const_val
                break
            
    return linearly_independent, idx


class LUDynamicEqStore:
    def __init__(self, integrated_constants=False):
        self.comb_to_idx = {}  # mapping of monomial -> index
        self.idx_to_comb = {}  # mapping of index -> 
        self.equation_ids = {} # mapping of equation to identifier (usually a clock cycle)
        if integrated_constants and tuple() not in self.comb_to_idx:
            self.comb_to_idx[tuple()] = len(self.comb_to_idx)
            self.idx_to_comb[len(self.comb_to_idx)] = tuple()

        self.num_vars = len(self.comb_to_idx)
        self.num_eqs = 0

        self.linked_stores = set([self])

        self.lower_matrix = np.eye(256, dtype = 'uint8')
        self.upper_matrix = np.eye(256, dtype = 'uint8')
        self.constants = np.zeros([256], dtype = 'uint8')
        self.solved_for = np.zeros([256], dtype = 'uint8')


    def insert_equation(self, equation, extra_const = 0, identifier=None, translate_ANF = True):
        if translate_ANF:
            equation_anf = equation.translate_ANF()
        else:
            equation_anf = equation

        coef_vector = np.zeros([self.lower_matrix.shape[1]], dtype=np.uint8)
        const_val = extra_const
        for term in equation_anf.args:

            # handle constant values
            if type(term) == CONST:
                if tuple() in self.comb_to_idx:
                    coef_vector[self.comb_to_idx[tuple()]] = term.value
                else:
                    const_val ^= term.value
                continue

            comb = tuple(sorted([var.index for var in term.args]))

            if comb not in self.comb_to_idx:
                # expand matrices for new variables as necessary:
                # because LU only adds independent vectors
                # num_eqs <= num_vars always, and so only
                # the one resize (new variables) is needed.
                if self.num_vars == self.lower_matrix.shape[1]:
                    # initialize
                    new_upper_matrix = np.eye(self.num_vars * 2,  dtype=np.uint8)
                    new_lower_matrix = np.eye(self.num_vars * 2,  dtype=np.uint8)
                    new_constants = np.zeros([self.num_vars * 2], dtype=np.uint8)
                    new_solved_for = np.zeros([self.num_vars * 2], dtype=np.uint8)
                    new_coef_vector = np.zeros([self.num_vars * 2], dtype=np.uint8)

                    # copy
                    new_upper_matrix[:self.num_vars, :self.num_vars] = self.upper_matrix
                    new_lower_matrix[:self.num_vars, :self.num_vars] = self.lower_matrix
                    new_constants[:self.num_vars] = self.constants
                    new_solved_for[:self.num_vars] = self.solved_for
                    new_coef_vector[:self.num_vars] = coef_vector

                    # replace
                    self.upper_matrix = new_upper_matrix
                    self.lower_matrix = new_lower_matrix
                    self.constants = new_constants
                    self.solved_for = new_solved_for
                    coef_vector = new_coef_vector

                # insert variable into variable maps
                self.idx_to_comb[self.num_vars] = comb
                self.comb_to_idx[comb] = self.num_vars
                self.num_vars += 1

                # also send the variable to any linked stores
                for store in self.linked_stores:
                    if store != self:
                        store._update_from_linked(comb)

            # after expanding as necessary, still set the appropriate var:
            coef_vector[self.comb_to_idx[comb]] = 1

        linearly_independent, insertion_idx =  LU_reduction(
            self.upper_matrix,self.lower_matrix,self.constants,self.solved_for, self.num_vars,
            coef_vector,const_val,
        )

        # have to insert identifier outside of numba optimized loop
        # this is because equation_ids is an untyped dictionary
        if linearly_independent:
            self.equation_ids[insertion_idx] = identifier
            self.num_eqs += 1
        return linearly_independent
    
    # attest this class is a dynamic store
    @classmethod
    def _is_dynamic_store(cls):
        return True
    
    # link to another dynamic store
    def link(self, other_store):
        if not type(other_store)._is_dynamic_store():
            raise ValueError("Can only link dynamic stores")

        self.linked_stores |= other_store.linked_stores
        other_store.linked_stores |= self.linked_stores

    # recieve a new variable from a linked store:
    def _update_from_linked(self, comb):
        if comb not in self.comb_to_idx:
            if self.num_vars == self.lower_matrix.shape[1]:
                # initialize
                new_upper_matrix = np.eye(self.num_vars * 2,  dtype=np.uint8)
                new_lower_matrix = np.eye(self.num_vars * 2,  dtype=np.uint8)
                new_constants = np.zeros([self.num_vars * 2], dtype=np.uint8)
                new_solved_for = np.zeros([self.num_vars * 2], dtype=np.uint8)

                # copy
                new_upper_matrix[:self.num_vars, :self.num_vars] = self.upper_matrix
                new_lower_matrix[:self.num_vars, :self.num_vars] = self.lower_matrix
                new_constants[:self.num_vars] = self.constants
                new_solved_for[:self.num_vars] = self.solved_for

                # replace
                self.upper_matrix = new_upper_matrix
                self.lower_matrix = new_lower_matrix
                self.constants = new_constants
                self.solved_for = new_solved_for

            # insert variable into variable maps
            self.idx_to_comb[self.num_vars] = comb
            self.comb_to_idx[comb] = self.num_vars
            self.num_vars += 1

    @property
    def rank(self):
        return self.num_eqs
