from PyPR.BooleanLogic import BooleanFunction
from PyPR.BooleanLogic import CONST

import numpy as np

import numba
u8 = numba.types.uint8
u64 = numba.types.uint64
b1 = numba.types.b1

from typing import Any

@numba.njit(numba.types.Tuple((b1,u64))(u8[:,:],u8[:,:],u8[:],u8[:],u8[:],u64,u64))
def _LU_reduction_consistent(
    upper: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    lower: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    extra_constants: np.ndarray[tuple[int],np.dtype[np.uint8]],
    solved_for: np.ndarray[tuple[int],np.dtype[np.uint8]],
    coef_vector: np.ndarray[tuple[int],np.dtype[np.uint8]],
    num_vars: int,
    const_idx: int
) -> tuple[bool, int]:
    linearly_independent = False
    coef_vector = coef_vector.copy()
    modification_vector = np.zeros_like(coef_vector)
    
    initial_constant = coef_vector[const_idx]
    for idx in range(const_idx):
        # skip coefficients which don't need to be cancelled:
        if coef_vector[idx] == 0:
            continue

        modification_vector[idx] = 1

        if solved_for[idx]:
            # update coefficient vectors
            for i in range(idx,num_vars):
                coef_vector[i] ^= upper[idx,i]  
        else:
            linearly_independent = True
            solved_for[idx] = 1
            upper[idx] = coef_vector
            lower[idx] = modification_vector

            # top of extra constants stores the augmented constants
            # accumulated constants are naturally computed in the upper matrix
            extra_constants[idx] = initial_constant
            break

    # cancel constants as if they were at the end, to enable checking consistency
    if linearly_independent:
        return linearly_independent, idx

    derived_constant = coef_vector[const_idx]
    for idx in range(const_idx+1,len(coef_vector)):
        # skip coefficients which don't need to be cancelled:
        if coef_vector[idx] == 0:
            continue

        modification_vector[idx] = 1

        if solved_for[idx]:
            coef_vector[const_idx] ^= extra_constants[idx]
            for i in range(idx,num_vars):
                coef_vector[i] ^= upper[idx,i]  
            
        else:
            linearly_independent = True
            solved_for[idx] = 1
            upper[idx] = coef_vector
            lower[idx] = modification_vector
            
            # bottom of extra constants stores the augmented constants
            # the derived constants are stored in the lower matrix to allow proper reconstruction.
            extra_constants[idx] = coef_vector[const_idx]
            lower[idx, const_idx] = derived_constant
            upper[idx, const_idx] = 0
            break

    if (
        not linearly_independent and 
        coef_vector[const_idx] == 1   # accumulated constant == 1
    ):
        raise ValueError("Inconsistent!")

    return linearly_independent, idx

@numba.njit(numba.types.Tuple((b1,u64))(u8[:,:],u8[:,:],u8[:],u8[:],u64))
def _LU_reduction(
    upper: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    lower: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    solved_for: np.ndarray[tuple[int],np.dtype[np.uint8]],
    coef_vector: np.ndarray[tuple[int],np.dtype[np.uint8]],
    num_vars: int,
) -> tuple[bool, int]:
    """This function performs LU reduction, and inserts the results into the corresponding matrix.

    This variant is for the inconsistent version. CONST(1) is treated the same as any other monomial,
    and is simply decomposed into L and U (extra constants are not used, and remain all zero).

    :param upper: The Upper matrix of the Equation store. The reduced coefficients will be
        inserted here
    :type upper: 
    :param lower: The Lower matrix of the Equation store. The reduction coefficients will be
        inserted here
    :type lower: _type_
    :param solved_for: The vector describing which variables have an associated equation.
    :type solved_for: _type_
    :param num_vars: The number of variables the equation store is currently tracking (which
        might differ from the length of the coefficient vector). This allows the loop to be
        shorter for dynamic equation stores.
    :type num_vars: int
    :param coef_vector: The vector of coefficients for the equation being inserted
    :type coef_vector: _type_
    :return: Whether or not the equation was inserted, and if so, at what index 
        (if not inserted the index will returned but has no meaning and shouldn't be used)
    :rtype: tuple[bool, int]
    """
    linearly_independent = False
    coef_vector = coef_vector.copy()
    modification_vector = np.zeros_like(coef_vector)
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
                break

    return linearly_independent, idx

class LUEqStore:
    def _enforce_consistency(self):
        if (self.dynamic == False and tuple() not in self.comb_to_idx):
            raise ValueError(
                "Consistency checks are only enabled when CONST(1) has an associated index. " +
                "Please update your variable map if you wish to enable this feature, or " +
                "set the flag `consistent = False` to ignore it."
            )
        
        # set consistency flag:
        self.consistent=True

        # add constant column in index maps if necessary:
        if not tuple() in self.comb_to_idx:
            self.comb_to_idx[tuple()] = self.num_vars
            self.idx_to_comb[self.num_vars] = tuple()
            self.num_vars += 1

        # add equation by marking solved if necessary:
        if not self.solved_for[self.comb_to_idx[tuple()]]:
            self.solved_for[self.comb_to_idx[tuple()]] = 1
            # self.equation_ids[self.num_eqs] = "CONSISTENCY"
            self.num_eqs += 1
        else:
            raise ValueError("Haven't implemented the case where const already has a definition.")

    def __init__(
        self, 
        comb_to_idx: dict[tuple[int,...], int] | None = None,
        consistent: bool = False
    ):        
        self.consistent = consistent
        self.linked_stores = set([self])

        # Dynamic Stores
        if comb_to_idx == None:
            self.dynamic = True
            self.comb_to_idx = {}  # mapping of monomial -> index
            self.idx_to_comb = {}  # mapping of index -> monomial

            self.equation_ids = {}
            self.num_vars = len(self.comb_to_idx)
            self.num_eqs = 0

            self.lower_matrix = np.eye(256, dtype = 'uint8')
            self.upper_matrix = np.eye(256, dtype = 'uint8')
            self.solved_for = np.zeros([256], dtype = 'uint8')
            self._extra_constants = np.zeros([256], dtype = 'uint8')

            if consistent:
                self._enforce_consistency()
        
        # Static Stores
        else:
            self.dynamic = False
            self.comb_to_idx = {k:v for k,v in comb_to_idx.items()}
            self.idx_to_comb = {v:k for k,v in comb_to_idx.items()}
        
            self.equation_ids = {}
            self.num_vars = len(self.comb_to_idx)
            self.num_eqs = 0

            self.lower_matrix = np.eye(self.num_vars, dtype = 'uint8')
            self.upper_matrix = np.eye(self.num_vars, dtype = 'uint8')
            self.solved_for = np.zeros([self.num_vars], dtype = 'uint8')
            self._extra_constants = np.zeros([self.num_vars], dtype = 'uint8')

            if consistent:
                self._enforce_consistency()

    def insert_equation(
        self,
        equation: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]], 
        identifier: Any = None, 
        translate_ANF: bool = True
    ): 
        # Allow accepting an array directly when the indices are known
        if type(equation) == np.ndarray:
            if self.dynamic:
                raise ValueError(
                    "Cannot insert an ndarray into a dynamic equation store. " +
                    "An index mapping must be provided in order to enable this feature."
                )
            if equation.shape != (self.num_vars,):
                raise ValueError(
                    f"Equation store has {self.num_vars} variables, " +
                    f"but the provided ndarray has shape {equation.shape}."
                )
            else:
                coef_vector = equation
            
        elif isinstance(equation, BooleanFunction):
            # translate equations to ANF if necessary:
            if translate_ANF:
                equation_anf = equation.translate_ANF()
            else:
                equation_anf = equation

            coef_vector = np.zeros([self.lower_matrix.shape[1]], dtype=np.uint8)
            for term in equation_anf.args:

                # handle consts if necessary:
                if type(term) != CONST:
                    comb = tuple(sorted([var.index for var in term.args])) #type: ignore
                elif term.value == 1:
                    comb = tuple()
                elif term.value == 0:
                    continue
                
                # update maps / expand if necessary:
                if comb not in self.comb_to_idx:
                    self._update_known_monomials(comb)
                
                    # update coef vector size if needed:
                    if len(coef_vector) != self.lower_matrix.shape[1]:
                        new_coefs = np.zeros([self.lower_matrix.shape[1]], dtype=np.uint8)
                        new_coefs[:len(coef_vector)] = coef_vector
                        coef_vector = new_coefs

                # after expanding as necessary, still set the appropriate var:
                coef_vector[self.comb_to_idx[comb]] = 1
        else:
            raise TypeError(
                f"Equation must be a BooleanFunction or ndarray, not {type(equation)}."
            )
        
        # insert / LU reduce
        if (
            len(coef_vector) > min(self.upper_matrix.shape[0], self.lower_matrix.shape[0]) or
            len(coef_vector) > min(self.upper_matrix.shape[1], self.lower_matrix.shape[1]) or
            len(coef_vector) > len(self._extra_constants) or
            len(coef_vector) > len(self.solved_for)
        ):
            print(self.upper_matrix.shape, self.lower_matrix.shape)
            print(self._extra_constants.shape, self.solved_for.shape)
            print(coef_vector.shape)
            raise ValueError("COEF TOO LONG!!!")
        if self.consistent:
            linearly_independent, insertion_idx = _LU_reduction_consistent(
                self.upper_matrix,
                self.lower_matrix,
                self._extra_constants,
                self.solved_for, 
                coef_vector,
                self.num_vars,
                self.comb_to_idx[tuple()]
            )
        else:
            linearly_independent, insertion_idx =  _LU_reduction(
                self.upper_matrix,
                self.lower_matrix,
                self.solved_for, 
                coef_vector,
                self.num_vars,  
            )
        
        # update fields which can't be modified in numba function
        if linearly_independent:
            self.equation_ids[insertion_idx] = identifier
            self.num_eqs += 1

        return linearly_independent

    # link stores:
    def link(self, other_store):
        self.linked_stores.add(other_store)

        # immediately after linking, send over all known monomials:
        for monomial in self.comb_to_idx:
            other_store._update_known_monomials(monomial)
        
    # recieve a new variable from a linked store:
    def _update_known_monomials(self, comb):
        if comb not in self.comb_to_idx:
            if not self.dynamic:
                raise ValueError(f"Monomial {comb} not in index maps for static {type(self)}")
           
            if self.num_vars == self.lower_matrix.shape[1]:
                # initialize
                new_upper_matrix = np.eye(self.num_vars * 2,  dtype=np.uint8)
                new_lower_matrix = np.eye(self.num_vars * 2,  dtype=np.uint8)
                new_solved_for = np.zeros([self.num_vars * 2], dtype=np.uint8)
                new_extra_consts = np.zeros([self.num_vars * 2], dtype=np.uint8)

                # copy
                new_upper_matrix[:self.num_vars, :self.num_vars] = self.upper_matrix
                new_lower_matrix[:self.num_vars, :self.num_vars] = self.lower_matrix
                new_solved_for[:self.num_vars] = self.solved_for
                new_extra_consts[:self.num_vars] = self._extra_constants

                # replace
                self.upper_matrix = new_upper_matrix
                self.lower_matrix = new_lower_matrix
                self.solved_for = new_solved_for
                self._extra_constants = new_extra_consts

            # insert variable into variable maps
            self.idx_to_comb[self.num_vars] = comb
            self.comb_to_idx[comb] = self.num_vars
            self.num_vars += 1

            # recursively pass monomials along to all other linked stores
            # doing this only when "surprised" prevents infinite recursion
            for store in self.linked_stores:
                if not (store is self):
                    store._update_known_monomials(comb)

    @property
    def rank(self):
        return self.num_eqs
