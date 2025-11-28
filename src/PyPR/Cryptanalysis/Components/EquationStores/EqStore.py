from typing import Any
from PyPR.BooleanLogic import BooleanFunction, CONST

import numpy as np
from numpy.typing import NDArray

import numpy as np


class EqStore:
    def __init__(self, 
        comb_to_idx = None,
        consistent = False
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

            self.equations = np.zeros([256,256], dtype = 'uint8')
     
        # Static Stores
        else:
            self.dynamic = False
            self.comb_to_idx = {k:v for k,v in comb_to_idx.items()}
            self.idx_to_comb = {v:k for k,v in comb_to_idx.items()}
        
            self.equation_ids = {}
            self.num_vars = len(self.comb_to_idx)
            self.num_eqs = 0

            self.equations = np.zeros([256,self.num_vars], dtype = 'uint8')

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

            coef_vector = np.zeros([self.equations.shape[1]], dtype=np.uint8)
            for term in equation_anf.args:

                # handle consts if necessary:
                if type(term) != CONST:
                    comb = tuple(sorted([var.index for var in term.args])) #type: ignore
                elif term.value == 1:
                    comb = tuple()
                elif term.value == 0:
                    continue
                
                # update maps / expand if necessary:
                # additional check prevents function call most times, at the cost of
                # an extra check when triggered - this should be beneficial on avg.
                if comb not in self.comb_to_idx:
                    self._update_known_monomials(comb)
                
                    # expand coef vector size as necessary:
                    if len(coef_vector) != self.equations.shape[1]:
                        new_coefs = np.zeros([self.equations.shape[1]], dtype=np.uint8)
                        new_coefs[:len(coef_vector)] = coef_vector
                        coef_vector = new_coefs

                # after expanding variables as necessary, still set the appropriate var:
                coef_vector[self.comb_to_idx[comb]] = 1
         
        else:
            raise TypeError(
                f"Equation must be a BooleanFunction or ndarray, not {type(equation)}."
            )
        
        # expand number of equations as necessary:
        if self.num_eqs == self.equations.shape[0]:
            new_equations = np.zeros([
                self.equations.shape[0] * 2, 
                self.equations.shape[1]
            ],  dtype=np.uint8)
            new_equations[:self.num_eqs,:self.num_vars] = self.equations[:self.num_eqs,:self.num_vars]
            self.equations = new_equations

        self.equation_ids[self.num_eqs] = identifier
        self.equations[self.num_eqs] = coef_vector
        self.num_eqs += 1
        return True
    
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
            
            # expand variables if needed:
            if self.num_vars == self.equations.shape[1]:
                new_equations = np.zeros([
                    self.equations.shape[0],
                    self.equations.shape[1] * 2
                ], dtype=np.uint8)

                new_equations[:self.num_eqs,:self.num_vars] = self.equations[:self.num_eqs,:self.num_vars]
                self.equations = new_equations

            # insert variable into variable maps
            self.idx_to_comb[self.num_vars] = comb
            self.comb_to_idx[comb] = self.num_vars
            self.num_vars += 1

            # recursively pass monomials along to all other linked stores
            # doing this only when "surprised" prevents infinite recursion
            for store in self.linked_stores:
                if not store is self:
                    store._update_known_monomials(comb)
