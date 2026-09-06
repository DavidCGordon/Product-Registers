from PyPR.BooleanLogic import BooleanFunction, BooleanANF

from PyPR.Cryptanalysis.Components.EquationStores.IndexedEqStore import IndexedEqStore
from PyPR.Cryptanalysis.Components.Adapters.equation_repr import (
    extract_monomials,
    coef_vector_to_anf,
)
import numpy as np

from typing import Any


class SymbolicEqStore(IndexedEqStore):
    """Bag of BooleanANF equations — the symbolic analog of EqStore.

    Equations are stored as a list of :class:`BooleanANF` with no
    reduction or solving on insertion. To solve, feed the equations
    to a batch solver (Grob_Solver) or convert to a matrix via
    :func:`~PyPR.Cryptanalysis.Components.Adapters.store_repr.to_coef_matrix`.

    Maintains ``comb_to_idx`` / ``idx_to_comb`` index maps inherited
    from :class:`IndexedEqStore`, enabling linking with matrix-based
    stores and conversion via the converters module.
    """

    equations: list[BooleanANF]

    def __init__(
        self,
        comb_to_idx: dict[tuple[int, ...], int] | None = None,
        consistent: bool = False,
    ):
        super().__init__(comb_to_idx, consistent)
        self.equations: list[BooleanANF] = []

    def _expand_storage(self):
        """No-op — symbolic storage (list) grows naturally."""
        pass

    def insert_equation(
        self,
        equation: BooleanFunction | np.ndarray[tuple[int], np.dtype[np.uint8]],
        identifier: Any = None,
        translate_ANF: bool = True,
    ):
        """Insert an equation into the symbolic bag.

        :param equation: A BooleanFunction (from any generator) or
            ndarray coefficient vector (static mode only).
        :type equation: BooleanFunction | np.ndarray
        :param identifier: Optional metadata (clock time, cube, etc.).
        :type identifier: Any, optional
        :param translate_ANF: Whether to convert BooleanFunction to ANF
            before storing. Defaults to True.
        :type translate_ANF: bool, optional
        :return: True (always succeeds — bags don't reject equations).
        :rtype: bool
        """
        if type(equation) == np.ndarray:
            if self.dynamic:
                raise ValueError(
                    "Cannot insert an ndarray into a dynamic equation store. "
                    + "An index mapping must be provided in order to enable this feature."
                )
            if equation.shape != (self.num_vars,):
                raise ValueError(
                    f"Equation store has {self.num_vars} variables, "
                    + f"but the provided ndarray has shape {equation.shape}."
                )
            anf = coef_vector_to_anf(equation, self.idx_to_comb)

        elif isinstance(equation, BooleanFunction):
            for comb in extract_monomials(equation, translate_ANF):
                if comb not in self.comb_to_idx:
                    self._update_known_monomials(comb)

            anf = BooleanANF.from_BooleanFunction(
                equation.translate_ANF() if translate_ANF else equation
            )
        else:
            raise TypeError(
                f"Equation must be a BooleanFunction or ndarray, not {type(equation)}."
            )

        self.equation_ids[self.num_eqs] = identifier
        self.equations.append(anf)
        self.num_eqs += 1
        return True
