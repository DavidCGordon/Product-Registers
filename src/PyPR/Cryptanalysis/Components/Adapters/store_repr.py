"""Store-level representation converters.

Extract equations from any equation store as either a coefficient matrix
(ndarray) or a list of BooleanANF, regardless of the store's native
representation.
"""
import numpy as np
from numpy.typing import NDArray

from PyPR.BooleanLogic import BooleanANF
from PyPR.Cryptanalysis.Components.Adapters.equation_repr import coef_vector_to_anf


def to_coef_matrix(store) -> tuple[NDArray[np.uint8], dict, dict]:
    """Extract a coefficient matrix from any equation store.

    Returns ``(matrix, comb_to_idx, idx_to_comb)`` where ``matrix`` has
    shape ``(num_eqs, num_vars)`` with one equation per row.

    Handles each store type by duck-typing:

    - **EqStore** (has ``.equations`` as ndarray + ``.comb_to_idx``):
      direct extraction.
    - **LUEqStore** (has ``.upper_matrix`` + ``.solved_for``):
      extract rows from ``upper_matrix`` where ``solved_for[i] == 1``.
    - **GroebnerEqStore** (has ``.equations`` as list of BooleanANF,
      no ``.comb_to_idx``): enumerate all monomials, build index maps,
      construct matrix.
    - **SymbolicEqStore** (has ``.equations`` as list of BooleanANF
      + ``.comb_to_idx``): use existing index maps, construct matrix.

    :param store: Any equation store.
    :return: ``(matrix, comb_to_idx, idx_to_comb)``
    :rtype: tuple[NDArray[np.uint8], dict, dict]
    """
    # LUEqStore path: upper_matrix rows where solved_for == 1
    if hasattr(store, 'upper_matrix') and hasattr(store, 'solved_for'):
        c2i = dict(store.comb_to_idx)
        i2c = dict(store.idx_to_comb)
        n = store.num_vars

        # For consistent stores the const_idx row is an axiom (CONST=1),
        # not a real equation — drop it and let solvers/guessing handle it.
        const_idx: int | None = c2i.get(tuple())

        pivot_rows = [i for i in range(n) if store.solved_for[i]]
        if store.consistent and const_idx is not None:
            pivot_rows = [i for i in pivot_rows if i != const_idx]

        if pivot_rows:
            matrix = np.stack([store.upper_matrix[i, :n] for i in pivot_rows])
            # _LU_reduction_consistent zeros U[idx, const_idx] for idx > const_idx
            # and stores the accumulated constant in _extra_constants[idx].
            # Restore so extracted equations are complete.
            if store.consistent and const_idx is not None:
                for row_idx, i in enumerate(pivot_rows):
                    if i > const_idx:
                        matrix[row_idx, const_idx] = store._extra_constants[i]
        else:
            matrix = np.zeros((0, n), dtype=np.uint8)
        return matrix, c2i, i2c

    # EqStore path: equations is an ndarray
    if (
        hasattr(store, 'equations') and
        isinstance(store.equations, np.ndarray) and
        hasattr(store, 'comb_to_idx')
    ):
        c2i = dict(store.comb_to_idx)
        i2c = dict(store.idx_to_comb)
        matrix = store.equations[:store.num_eqs, :store.num_vars].copy()
        return matrix, c2i, i2c

    # Symbolic path: equations is a list of BooleanANF
    if hasattr(store, 'equations') and isinstance(store.equations, list):
        # Use existing index maps if available (SymbolicEqStore)
        if hasattr(store, 'comb_to_idx') and store.comb_to_idx:
            c2i = dict(store.comb_to_idx)
            i2c = dict(store.idx_to_comb)
        else:
            # Build index maps from scratch (GroebnerEqStore)
            all_monomials: set[frozenset[int]] = set()
            for eq in store.equations:
                all_monomials |= eq.terms
            sorted_monomials = sorted(
                all_monomials,
                key=lambda m: (len(m), tuple(sorted(m)))
            )
            c2i = {tuple(sorted(m)): i for i, m in enumerate(sorted_monomials)}
            i2c = {i: tuple(sorted(m)) for i, m in enumerate(sorted_monomials)}

        n = len(c2i)
        matrix = np.zeros((len(store.equations), n), dtype=np.uint8)
        for row, eq in enumerate(store.equations):
            for term in eq.terms:
                comb = tuple(sorted(term))
                if comb in c2i:
                    matrix[row, c2i[comb]] = 1
        return matrix, c2i, i2c

    raise TypeError(f"Cannot extract coefficient matrix from {type(store)}")

def to_anf_list(store) -> list[BooleanANF]:
    """Extract equations as a list of BooleanANF from any equation store.

    Handles each store type by duck-typing:

    - **GroebnerEqStore / SymbolicEqStore** (equations is list of
      BooleanANF): direct extraction.
    - **EqStore** (equations is ndarray + has idx_to_comb): convert
      each row to BooleanANF.
    - **LUEqStore** (upper_matrix + solved_for + idx_to_comb): convert
      pivot rows to BooleanANF.

    :param store: Any equation store.
    :return: List of BooleanANF equations.
    :rtype: list[BooleanANF]
    """
    # Symbolic path: already BooleanANF
    if hasattr(store, 'equations') and isinstance(store.equations, list):
        return list(store.equations)

    # LUEqStore path
    if hasattr(store, 'upper_matrix') and hasattr(store, 'solved_for'):
        i2c = store.idx_to_comb
        n = store.num_vars
        const_idx: int | None = store.comb_to_idx.get(tuple())

        result = []
        for i in range(n):
            if not store.solved_for[i]:
                continue
            # Skip the const axiom row for consistent stores
            if store.consistent and const_idx is not None and i == const_idx:
                continue

            row = store.upper_matrix[i, :n].copy()
            # Restore zeroed-out constant for post-const_idx rows
            if store.consistent and const_idx is not None and i > const_idx:
                row[const_idx] = store._extra_constants[i]

            result.append(coef_vector_to_anf(row, i2c))
        return result

    # EqStore path: ndarray rows
    if (
        hasattr(store, 'equations') and
        isinstance(store.equations, np.ndarray) and
        hasattr(store, 'idx_to_comb')
    ):
        i2c = store.idx_to_comb
        n = store.num_vars
        return [
            coef_vector_to_anf(store.equations[row, :n], i2c)
            for row in range(store.num_eqs)
        ]

    raise TypeError(f"Cannot extract BooleanANF list from {type(store)}")
