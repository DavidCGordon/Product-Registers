from typing import Any
from PyPR.BooleanLogic import BooleanFunction

import numpy as np
from numpy.typing import NDArray

from PyPR.Cryptanalysis.Components.EquationStores.IndexedEqStore import IndexedEqStore
from PyPR.Cryptanalysis.Components.Adapters.equation_repr import (
    extract_monomials,
    boolean_function_to_coef_vector,
)


class EqStore(IndexedEqStore):
    """Bag of coefficient vectors — no reduction or solving on insertion.

    Equations are stored as rows of a dense coefficient matrix.
    Each column corresponds to a monomial tracked by the index maps
    inherited from :class:`IndexedEqStore`.
    """

    def __init__(self,
        comb_to_idx = None,
        consistent = False
    ):
        super().__init__(comb_to_idx, consistent)

        if self.dynamic:
            self.equations = np.zeros([256, 256], dtype='uint8')
        else:
            self.equations = np.zeros([256, self.num_vars], dtype='uint8')

    def _expand_storage(self):
        """Double the column capacity of the equations matrix."""
        if self.num_vars < self.equations.shape[1]:
            return

        new_equations = np.zeros([
            self.equations.shape[0],
            self.equations.shape[1] * 2
        ], dtype=np.uint8)
        new_equations[:self.num_eqs, :self.num_vars] = (
            self.equations[:self.num_eqs, :self.num_vars]
        )
        self.equations = new_equations

    def insert_equation(
        self,
        equation: BooleanFunction | np.ndarray[tuple[int], np.dtype[np.uint8]],
        identifier: Any = None,
        translate_ANF: bool = True
    ):
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
            coef_vector = equation

        elif isinstance(equation, BooleanFunction):
            for comb in extract_monomials(equation, translate_ANF):
                if comb not in self.comb_to_idx:
                    self._update_known_monomials(comb)

            coef_vector = np.zeros(self.equations.shape[1], dtype=np.uint8)
            coef_vector[:self.num_vars] = boolean_function_to_coef_vector(
                equation, self.comb_to_idx, self.num_vars, translate_ANF
            )
        else:
            raise TypeError(
                f"Equation must be a BooleanFunction or ndarray, not {type(equation)}."
            )

        if self.num_eqs == self.equations.shape[0]:
            new_equations = np.zeros([
                self.equations.shape[0] * 2,
                self.equations.shape[1]
            ], dtype=np.uint8)
            new_equations[:self.num_eqs, :self.num_vars] = (
                self.equations[:self.num_eqs, :self.num_vars]
            )
            self.equations = new_equations

        self.equation_ids[self.num_eqs] = identifier
        self.equations[self.num_eqs] = coef_vector
        self.num_eqs += 1
        return True
