from PyPR.BooleanLogic import CONST

import numpy as np

class DynamicEqStore:
    def __init__(self, integrated_constants=False):
        self.comb_to_idx = {} # mapping of monomial -> index
        self.idx_to_comb = {} # mapping of index -> monomial
        self.equation_ids = {} # mapping of equation to identifier (usually a clock cycle)
        if integrated_constants and tuple() not in self.comb_to_idx:
            self.comb_to_idx[tuple()] = len(self.comb_to_idx)
            self.idx_to_comb[len(self.comb_to_idx)] = tuple()

        self.num_vars = len(self.comb_to_idx)
        self.num_eqs = 0
        
        self.linked_stores = set([self])

        self.equations = np.zeros([256,256], dtype = 'uint8')
        self.constants = np.zeros([256], dtype = 'uint8')


    def insert_equation(self, equation, extra_const = 0, identifier=None, translate_ANF = True):
        if translate_ANF:
            equation_anf = equation.translate_ANF()
        else:
            equation_anf = equation

        coef_vector = np.zeros([self.equations.shape[1]], dtype=np.uint8)
        const_val = extra_const
        for term in equation_anf.args:

            # handle constant values
            if type(term) == CONST:
                if tuple() in self.comb_to_idx:
                    coef_vector[self.comb_to_idx[tuple()]] = term.value
                else:
                    const_val ^= term.value
                continue
            try:
                comb = tuple(sorted([var.index for var in term.args]))
            except:
                print(term)
                raise KeyboardInterrupt

            # insert variable into self data
            if comb not in self.comb_to_idx:
                # expand matrices for new variables as necessary:
                if self.num_vars == self.equations.shape[1]:
                    # initialize
                    new_equations = np.zeros([self.equations.shape[0], self.num_vars * 2],  dtype=np.uint8)
                    new_coef_vector = np.zeros([self.num_vars * 2], dtype=np.uint8)

                    # copy
                    new_equations[:,:self.num_vars] = self.equations
                    new_coef_vector[:self.num_vars] = coef_vector

                    # replace
                    self.equations = new_equations
                    coef_vector = new_coef_vector

                # insert variable into variable maps
                self.idx_to_comb[self.num_vars] = comb
                self.comb_to_idx[comb] = self.num_vars
                self.num_vars += 1

                # also send the variable to any linked stores
                for store in self.linked_stores:
                    if store != self:
                        store._update_from_linked(comb)

            # after expanding variables as necessary, still set the appropriate var:
            coef_vector[self.comb_to_idx[comb]] = 1


        # expand equations as necessary:
        if self.num_eqs == self.equations.shape[0]:
            # initialize
            new_equations = np.zeros([self.num_eqs* 2, self.equations.shape[1]],  dtype=np.uint8)
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
            # expand matrices for new variables as necessary:
            if self.num_vars == self.equations.shape[1]:
                # initialize
                new_equations = np.zeros([self.equations.shape[0], self.num_vars * 2],  dtype=np.uint8)

                # copy
                new_equations[:,:self.num_vars] = self.equations

                # replace
                self.equations = new_equations

            # insert variable into variable maps
            self.idx_to_comb[self.num_vars] = comb
            self.comb_to_idx[comb] = self.num_vars
            self.num_vars += 1