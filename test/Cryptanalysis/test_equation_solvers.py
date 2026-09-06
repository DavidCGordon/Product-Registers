"""Cross-representation solver tests.

Verifies that each solver accepts stores of different types via
automatic conversion, producing consistent results across
representations.
"""
import numpy as np
import pytest

from PyPR.BooleanLogic import VAR, XOR, AND, CONST

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import SymbolicEqStore
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore

import PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver as LU_Solver
import PyPR.Cryptanalysis.Components.EquationSolving.GaussElim as GaussElim
import PyPR.Cryptanalysis.Components.EquationSolving.Grob_Solver as Grob_Solver


_COMB_TO_IDX = {
    (): 0,
    (0,): 1, (1,): 2, (2,): 3,
    (0, 1): 4, (0, 2): 5, (1, 2): 6,
}


# ── LU_Solver: native (consistent) ──────────────────────────────────────────

def test_lu_solver_native_consistent():
    """LU_Solver on its native LUEqStore with consistency mode."""
    store = LUEqStore(consistent=True)
    store.insert_equation(XOR(VAR(0), CONST(1)))
    store.insert_equation(VAR(1))
    store.insert_equation(XOR(VAR(2), CONST(1)))
    solution = LU_Solver.reduce(store)
    assert solution[store.comb_to_idx[(0,)]] == 1
    assert solution[store.comb_to_idx[(1,)]] == 0
    assert solution[store.comb_to_idx[(2,)]] == 1


# ── LU_Solver: cross-repr preserves conversion ──────────────────────────────

def test_lu_solver_from_eq_store_accepts():
    """LU_Solver accepts EqStore input without raising."""
    store = EqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    solution = LU_Solver.reduce(store)
    assert solution.shape == (len(_COMB_TO_IDX),)


def test_lu_solver_from_symbolic_store_accepts():
    """LU_Solver accepts SymbolicEqStore input without raising."""
    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    solution = LU_Solver.reduce(store)
    assert solution.shape == (len(_COMB_TO_IDX),)


def test_lu_cross_repr_matches_native():
    """LU_Solver via EqStore conversion gives the same result as native LUEqStore.

    Uses a homogeneous system (no constants) so consistency mode is irrelevant.
    """
    eqs = [VAR(0), VAR(1), XOR(VAR(0), VAR(2))]

    lu_store = LUEqStore(comb_to_idx=_COMB_TO_IDX)
    eq_store = EqStore(comb_to_idx=_COMB_TO_IDX)
    for eq in eqs:
        lu_store.insert_equation(eq)
        eq_store.insert_equation(eq)

    native_sol = LU_Solver.reduce(lu_store)
    cross_sol = LU_Solver.reduce(eq_store)
    np.testing.assert_array_equal(native_sol, cross_sol)


# ── GaussElim: reduce_matrix ────────────────────────────────────────────────

def test_reduce_matrix_identity():
    """reduce_matrix on an already-reduced identity matrix."""
    matrix = np.eye(3, dtype=np.uint8)
    rref, free = GaussElim.reduce_matrix(matrix)
    assert rref.shape[0] == 3
    assert np.array_equal(free, np.array([0, 0, 0], dtype=np.uint8))


def test_reduce_matrix_upper_triangular():
    """reduce_matrix on an upper triangular matrix with rank 3."""
    matrix = np.array([
        [1, 1, 0],
        [0, 1, 1],
        [0, 0, 1],
    ], dtype=np.uint8)
    rref, free = GaussElim.reduce_matrix(matrix)
    assert rref.shape[0] == 3
    np.testing.assert_array_equal(free, np.array([0, 0, 0], dtype=np.uint8))
    np.testing.assert_array_equal(rref, np.eye(3, dtype=np.uint8))


def test_reduce_matrix_dependent():
    """reduce_matrix on a rank-2 system (row 3 = row 1 XOR row 2)."""
    matrix = np.array([
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, 1],
    ], dtype=np.uint8)
    rref, free = GaussElim.reduce_matrix(matrix)
    assert rref.shape[0] == 2
    assert sum(free) == 1


# ── GaussElim: cross-repr acceptance ─────────────────────────────────────────

def test_gauss_accepts_eq_store():
    """GaussElim.reduce accepts EqStore natively."""
    store = EqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    solution = GaussElim.reduce(store)
    assert len(solution) == len(_COMB_TO_IDX)


def test_gauss_accepts_lu_store():
    """GaussElim.reduce accepts LUEqStore via conversion."""
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    solution = GaussElim.reduce(store)
    assert len(solution) == store.num_vars


def test_gauss_accepts_symbolic_store():
    """GaussElim.reduce accepts SymbolicEqStore via conversion."""
    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    solution = GaussElim.reduce(store)
    assert len(solution) == len(_COMB_TO_IDX)


# ── Grob_Solver: native + cross-representation ──────────────────────────────

def test_grob_from_eq_store():
    """Grob_Solver converts EqStore equations to BooleanANF."""
    store = EqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(XOR(VAR(0), CONST(1)))
    store.insert_equation(VAR(1))
    store.insert_equation(XOR(VAR(2), CONST(1)))
    gb = Grob_Solver.reduce(store)
    assert gb.solved_vars.get(0) == 1
    assert gb.solved_vars.get(1) == 0
    assert gb.solved_vars.get(2) == 1


def test_grob_from_symbolic_store():
    """Grob_Solver converts SymbolicEqStore equations to BooleanANF."""
    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(XOR(VAR(0), CONST(1)))
    store.insert_equation(VAR(1))
    store.insert_equation(XOR(VAR(2), CONST(1)))
    gb = Grob_Solver.reduce(store)
    assert gb.solved_vars.get(0) == 1
    assert gb.solved_vars.get(1) == 0
    assert gb.solved_vars.get(2) == 1


def test_grob_native_passthrough():
    """Grob_Solver returns the GroebnerEqStore as-is when given one."""
    store = GroebnerEqStore(simplify_mode=None)
    store.insert_equation(XOR(VAR(0), CONST(1)))
    store.insert_equation(VAR(1))
    result = Grob_Solver.reduce(store)
    assert result is store


def test_grob_from_lu_store():
    """Grob_Solver converts LUEqStore via to_anf_list."""
    store = LUEqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(XOR(VAR(0), CONST(1)))
    store.insert_equation(VAR(1))
    store.insert_equation(XOR(VAR(2), CONST(1)))
    gb = Grob_Solver.reduce(store)
    assert gb.solved_vars.get(0) == 1
    assert gb.solved_vars.get(1) == 0
    assert gb.solved_vars.get(2) == 1
