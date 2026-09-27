"""Cross-representation solver tests, plus the GF(2) elimination underneath them.

Verifies that each solver accepts stores of different types via
automatic conversion, producing consistent results across
representations.

The final section tests `GaussElim.reduce_matrix` directly against the
definition of a reduced row echelon form: the row space and rank are
determined by the input, so they cannot depend on elimination order or pivot
choice.  Each property is recomputed with an independent rank routine.
"""
import random

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
from PyPR.Cryptanalysis.Components.EquationSolving.GaussElim import reduce_matrix


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


# ── reduce_matrix as a GF(2) reduced row echelon form ─────────────────────────

@pytest.mark.parametrize("trial", range(10))
def test_reduce_matrix_returns_one_row_per_pivot_and_preserves_the_row_space(trial):
    """The RREF has exactly `rank` rows and spans the same space as the input.

    Row reduction is a change of basis for the row space, so stacking the input
    and its RREF must not raise the rank.  Combined with the row count, this
    says the returned matrix is a basis of the original row space -- the
    property every solver downstream depends on.
    """
    rng = random.Random(700 + trial)

    def gf2_rank(rows):
        rows, pivot_row = [int(r) for r in rows], 0
        for bit in range(max((r.bit_length() for r in rows), default=0)):
            pivot = next((i for i in range(pivot_row, len(rows)) if (rows[i] >> bit) & 1), None)
            if pivot is None:
                continue
            rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
            for i in range(len(rows)):
                if i != pivot_row and (rows[i] >> bit) & 1:
                    rows[i] ^= rows[pivot_row]
            pivot_row += 1
        return pivot_row

    height, width = rng.randint(2, 7), rng.randint(2, 7)
    matrix = np.array(
        [[rng.randint(0, 1) for _ in range(width)] for _ in range(height)], dtype=np.uint8
    )
    as_masks = [sum(int(matrix[i][j]) << j for j in range(width)) for i in range(height)]
    rank = gf2_rank(as_masks)

    reduced, _ = reduce_matrix(matrix.copy())
    reduced_masks = [
        sum(int(reduced[i][j]) << j for j in range(width)) for i in range(reduced.shape[0])
    ]

    assert reduced.shape[0] == rank, (
        f"trial {trial}: RREF has {reduced.shape[0]} rows for a rank-{rank} matrix"
    )
    assert gf2_rank(as_masks + reduced_masks) == rank, (
        f"trial {trial}: row reduction changed the row space"
    )


@pytest.mark.parametrize("trial", range(10))
def test_reduce_matrix_output_is_in_reduced_row_echelon_form(trial):
    """Pivots advance strictly left to right, and each pivot column is a unit vector.

    These two conditions are what makes the form *reduced* rather than merely
    triangular, and they are what lets a caller read a solution straight off
    the rows.  A repeated or out-of-order pivot column would silently produce
    wrong back-substitutions.
    """
    rng = random.Random(800 + trial)
    height, width = rng.randint(2, 7), rng.randint(2, 7)
    matrix = np.array(
        [[rng.randint(0, 1) for _ in range(width)] for _ in range(height)], dtype=np.uint8
    )

    reduced, _ = reduce_matrix(matrix.copy())
    pivots = [int(np.argmax(reduced[i])) for i in range(reduced.shape[0])]

    assert pivots == sorted(pivots) and len(set(pivots)) == len(pivots), (
        f"trial {trial}: pivot columns {pivots} are not strictly increasing"
    )
    for pivot in pivots:
        assert int(reduced[:, pivot].sum()) == 1, (
            f"trial {trial}: pivot column {pivot} is not a unit column"
        )


@pytest.mark.parametrize("trial", range(10))
def test_free_variables_are_exactly_the_non_pivot_columns(trial):
    """`free_vars` marks a column iff no row pivots there.

    The free/pivot split is the dimension count of the solution space: with
    `rank` pivot columns out of `width`, the system has `width - rank` free
    variables and 2^(width-rank) solutions.  Miscounting here would make a
    solver report the wrong number of candidate keys.
    """
    rng = random.Random(900 + trial)
    height, width = rng.randint(2, 7), rng.randint(2, 7)
    matrix = np.array(
        [[rng.randint(0, 1) for _ in range(width)] for _ in range(height)], dtype=np.uint8
    )

    reduced, free_vars = reduce_matrix(matrix.copy())
    pivots = {int(np.argmax(reduced[i])) for i in range(reduced.shape[0])}

    assert len(free_vars) == width, "free_vars should carry one flag per column"
    assert [j for j in range(width) if free_vars[j]] == [
        j for j in range(width) if j not in pivots
    ], f"trial {trial}: free_vars does not complement the pivot columns"
    assert int(free_vars.sum()) == width - reduced.shape[0], (
        f"trial {trial}: free variable count should be width - rank"
    )
