"""Tests for representation converters between coef vectors, BooleanANF, and BooleanFunction.

The converter module bridges the three equation representations used
by the cryptanalysis pipeline:
  - Coefficient vectors (ndarray uint8)
  - BooleanANF (frozenset of frozensets)
  - BooleanFunction (DAG of gates)

Round-trip correctness is the core invariant: converting in one direction
and back should yield an equivalent equation.
"""
import numpy as np
import pytest

from PyPR.BooleanLogic import BooleanFunction, BooleanANF, XOR, AND, VAR, CONST
from PyPR.Cryptanalysis.Components.Adapters import (
    extract_monomials,
    boolean_function_to_coef_vector,
    coef_vector_to_anf,
    coef_vector_to_boolean_function,
    to_coef_matrix,
    to_anf_list,
)


# ── Fixtures ─────────────────────────────────────────────────────────────────

# Index map: () → 0, (0,) → 1, (1,) → 2, (2,) → 3, (0,1) → 4, (0,2) → 5, (1,2) → 6
_COMB_TO_IDX = {
    (): 0,
    (0,): 1, (1,): 2, (2,): 3,
    (0, 1): 4, (0, 2): 5, (1, 2): 6,
}
_IDX_TO_COMB = {v: k for k, v in _COMB_TO_IDX.items()}
_NUM_VARS = len(_COMB_TO_IDX)


# ── extract_monomials ────────────────────────────────────────────────────────

def test_extract_monomials_linear():
    """Extract monomials from a linear function XOR(VAR(0), VAR(1))."""
    f = XOR(VAR(0), VAR(1))
    monos = extract_monomials(f)
    assert set(monos) == {(0,), (1,)}


def test_extract_monomials_with_constant():
    """Constant 1 in the function appears as the empty tuple."""
    f = XOR(VAR(0), CONST(1))
    monos = extract_monomials(f)
    assert set(monos) == {(0,), ()}


def test_extract_monomials_degree_2():
    """AND(VAR(0), VAR(1)) XOR VAR(2) has monomials (0,1) and (2,)."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    monos = extract_monomials(f)
    assert set(monos) == {(0, 1), (2,)}


def test_extract_monomials_const_zero_dropped():
    """CONST(0) terms are silently dropped."""
    f = XOR(VAR(0), CONST(0))
    monos = extract_monomials(f)
    assert set(monos) == {(0,)}


# ── boolean_function_to_coef_vector ──────────────────────────────────────────

def test_bf_to_coef_linear():
    """x0 XOR x1 → [0, 1, 1, 0, 0, 0, 0]."""
    f = XOR(VAR(0), VAR(1))
    vec = boolean_function_to_coef_vector(f, _COMB_TO_IDX, _NUM_VARS)
    expected = np.array([0, 1, 1, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_array_equal(vec, expected)


def test_bf_to_coef_with_constant():
    """x0 XOR 1 → [1, 1, 0, 0, 0, 0, 0]."""
    f = XOR(VAR(0), CONST(1))
    vec = boolean_function_to_coef_vector(f, _COMB_TO_IDX, _NUM_VARS)
    expected = np.array([1, 1, 0, 0, 0, 0, 0], dtype=np.uint8)
    np.testing.assert_array_equal(vec, expected)


def test_bf_to_coef_degree_2():
    """x0*x1 XOR x2 → [0, 0, 0, 1, 1, 0, 0]."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    vec = boolean_function_to_coef_vector(f, _COMB_TO_IDX, _NUM_VARS)
    expected = np.array([0, 0, 0, 1, 1, 0, 0], dtype=np.uint8)
    np.testing.assert_array_equal(vec, expected)


# ── coef_vector_to_anf ──────────────────────────────────────────────────────

def test_coef_to_anf_linear():
    """[0, 1, 1, 0, 0, 0, 0] → x0 XOR x1 as BooleanANF."""
    vec = np.array([0, 1, 1, 0, 0, 0, 0], dtype=np.uint8)
    anf = coef_vector_to_anf(vec, _IDX_TO_COMB)
    assert anf.terms == frozenset({frozenset({0}), frozenset({1})})


def test_coef_to_anf_with_constant():
    """[1, 1, 0, 0, 0, 0, 0] → 1 XOR x0 as BooleanANF."""
    vec = np.array([1, 1, 0, 0, 0, 0, 0], dtype=np.uint8)
    anf = coef_vector_to_anf(vec, _IDX_TO_COMB)
    # frozenset() represents the constant 1 term
    assert anf.terms == frozenset({frozenset(), frozenset({0})})


def test_coef_to_anf_zero():
    """All-zero vector → empty BooleanANF (the zero polynomial)."""
    vec = np.zeros(_NUM_VARS, dtype=np.uint8)
    anf = coef_vector_to_anf(vec, _IDX_TO_COMB)
    assert anf.terms == frozenset()


def test_coef_to_anf_degree_2():
    """[0, 0, 0, 1, 1, 0, 0] → x0*x1 XOR x2."""
    vec = np.array([0, 0, 0, 1, 1, 0, 0], dtype=np.uint8)
    anf = coef_vector_to_anf(vec, _IDX_TO_COMB)
    assert anf.terms == frozenset({frozenset({2}), frozenset({0, 1})})


# ── coef_vector_to_boolean_function ──────────────────────────────────────────

def test_coef_to_bf_evaluates_correctly():
    """Round-trip: coef_vector → BooleanFunction evaluates the same
    as directly constructing the function.

    x0*x1 XOR x2 XOR 1 should evaluate identically for all 8 inputs.
    """
    vec = np.array([1, 0, 0, 1, 1, 0, 0], dtype=np.uint8)
    bf = coef_vector_to_boolean_function(vec, _IDX_TO_COMB)
    direct = XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1))
    for x0 in range(2):
        for x1 in range(2):
            for x2 in range(2):
                inp = [x0, x1, x2]
                assert bf.eval(inp) == direct.eval(inp), (
                    f"Mismatch at input {inp}"
                )


# ── Round-trip: BooleanFunction → coef_vector → BooleanFunction ─────────────

def test_round_trip_bf_coef_bf():
    """Converting BooleanFunction → coef_vector → BooleanFunction
    yields a function with identical truth table.
    """
    original = XOR(AND(VAR(0), VAR(2)), VAR(1), CONST(1))
    vec = boolean_function_to_coef_vector(original, _COMB_TO_IDX, _NUM_VARS)
    reconstructed = coef_vector_to_boolean_function(vec, _IDX_TO_COMB)

    for x0 in range(2):
        for x1 in range(2):
            for x2 in range(2):
                inp = [x0, x1, x2]
                assert original.eval(inp) == reconstructed.eval(inp), (
                    f"Round-trip mismatch at input {inp}"
                )


def test_round_trip_coef_anf_coef():
    """Converting coef_vector → BooleanANF → coef_vector is lossless."""
    original_vec = np.array([1, 0, 1, 0, 1, 0, 1], dtype=np.uint8)
    anf = coef_vector_to_anf(original_vec, _IDX_TO_COMB)

    reconstructed = np.zeros(_NUM_VARS, dtype=np.uint8)
    for term in anf.terms:
        comb = tuple(sorted(term))
        reconstructed[_COMB_TO_IDX[comb]] = 1

    np.testing.assert_array_equal(reconstructed, original_vec)


# ── to_coef_matrix / to_anf_list with real stores ────────────────────────────

def test_to_coef_matrix_from_eq_store():
    """to_coef_matrix extracts the correct matrix from an EqStore."""
    from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore

    store = EqStore(comb_to_idx=_COMB_TO_IDX)
    f1 = XOR(VAR(0), VAR(1))
    f2 = XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1))
    store.insert_equation(f1)
    store.insert_equation(f2)

    matrix, c2i, i2c = to_coef_matrix(store)
    assert matrix.shape == (2, _NUM_VARS)
    assert c2i == _COMB_TO_IDX

    # f1 = x0 + x1 → row should have 1s at indices 1, 2
    assert matrix[0, c2i[(0,)]] == 1
    assert matrix[0, c2i[(1,)]] == 1
    assert matrix[0, c2i[()]  ] == 0

    # f2 = x0*x1 + x2 + 1 → row should have 1s at indices 0, 3, 4
    assert matrix[1, c2i[()]    ] == 1
    assert matrix[1, c2i[(2,)]  ] == 1
    assert matrix[1, c2i[(0, 1)]] == 1


def test_to_coef_matrix_from_lu_store():
    """to_coef_matrix extracts pivot rows from an LUEqStore."""
    from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

    store = LUEqStore(comb_to_idx=_COMB_TO_IDX)
    f1 = XOR(VAR(0), VAR(1))
    f2 = XOR(VAR(1), VAR(2))
    store.insert_equation(f1)
    store.insert_equation(f2)

    matrix, c2i, _ = to_coef_matrix(store)
    # LUEqStore rejects linearly dependent equations; we inserted 2
    # independent ones so expect 2 rows (though the upper_matrix
    # stores the reduced form, not the original equations)
    assert matrix.shape[0] == 2
    assert matrix.shape[1] == _NUM_VARS


def test_to_anf_list_from_eq_store():
    """to_anf_list converts EqStore matrix rows to BooleanANF."""
    from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore

    store = EqStore(comb_to_idx=_COMB_TO_IDX)
    f = XOR(VAR(0), VAR(1))
    store.insert_equation(f)

    anf_list = to_anf_list(store)
    assert len(anf_list) == 1
    assert anf_list[0].terms == frozenset({frozenset({0}), frozenset({1})})


def test_to_anf_list_from_grobner_store():
    """to_anf_list extracts equations directly from GroebnerEqStore."""
    from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import (
        GroebnerEqStore,
    )

    store = GroebnerEqStore(simplify_mode=None)
    f = XOR(VAR(0), VAR(1))
    store.insert_equation(f)

    anf_list = to_anf_list(store)
    assert len(anf_list) >= 1
    # The Gröbner store may have reduced the equation, but the
    # variable indices should be preserved
    all_vars = set()
    for anf in anf_list:
        for term in anf.terms:
            all_vars |= term
    assert all_vars <= {0, 1}


def test_to_coef_matrix_from_grobner_store():
    """to_coef_matrix builds index maps and matrix from a GroebnerEqStore."""
    from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import (
        GroebnerEqStore,
    )

    store = GroebnerEqStore(simplify_mode=None)
    f1 = XOR(VAR(0), VAR(1))
    f2 = XOR(VAR(1), VAR(2))
    store.insert_equation(f1)
    store.insert_equation(f2)

    matrix, c2i, i2c = to_coef_matrix(store)
    # Should have built index maps covering the monomials present
    assert matrix.shape[0] >= 1
    assert len(c2i) == matrix.shape[1]
    # Each row should be binary
    assert set(np.unique(matrix)) <= {0, 1}


# ── Cross-representation consistency ────────────────────────────────────────

def test_to_coef_matrix_from_symbolic_store():
    """to_coef_matrix uses existing index maps from a SymbolicEqStore."""
    from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import (
        SymbolicEqStore,
    )

    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    f1 = XOR(VAR(0), VAR(1))
    f2 = XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1))
    store.insert_equation(f1)
    store.insert_equation(f2)

    matrix, c2i, i2c = to_coef_matrix(store)
    assert matrix.shape == (2, _NUM_VARS)
    assert c2i == _COMB_TO_IDX

    assert matrix[0, c2i[(0,)]] == 1
    assert matrix[0, c2i[(1,)]] == 1
    assert matrix[0, c2i[()]] == 0

    assert matrix[1, c2i[()]] == 1
    assert matrix[1, c2i[(2,)]] == 1
    assert matrix[1, c2i[(0, 1)]] == 1


def test_to_coef_matrix_from_symbolic_store_dynamic():
    """to_coef_matrix works on dynamic SymbolicEqStore (index maps built on insertion)."""
    from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import (
        SymbolicEqStore,
    )

    store = SymbolicEqStore()
    store.insert_equation(XOR(VAR(0), VAR(1)))
    store.insert_equation(XOR(VAR(1), VAR(2)))

    matrix, c2i, i2c = to_coef_matrix(store)
    assert matrix.shape[0] == 2
    assert len(c2i) == matrix.shape[1]
    assert set(np.unique(matrix)) <= {0, 1}


def test_to_anf_list_from_symbolic_store():
    """to_anf_list extracts equations directly from SymbolicEqStore."""
    from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import (
        SymbolicEqStore,
    )

    store = SymbolicEqStore()
    f = XOR(VAR(0), VAR(1))
    store.insert_equation(f)

    anf_list = to_anf_list(store)
    assert len(anf_list) == 1
    assert anf_list[0].terms == frozenset({frozenset({0}), frozenset({1})})


def test_symbolic_store_round_trip_via_matrix():
    """SymbolicEqStore → to_coef_matrix → reconstruct BooleanANF matches original."""
    from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import (
        SymbolicEqStore,
    )

    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    f1 = XOR(AND(VAR(0), VAR(1)), VAR(2))
    f2 = XOR(VAR(0), CONST(1))
    store.insert_equation(f1)
    store.insert_equation(f2)

    matrix, c2i, i2c = to_coef_matrix(store)
    anf_list = to_anf_list(store)

    reconstructed = np.zeros_like(matrix)
    for row, anf in enumerate(anf_list):
        for term in anf.terms:
            comb = tuple(sorted(term))
            reconstructed[row, c2i[comb]] = 1
    np.testing.assert_array_equal(reconstructed, matrix)


def test_eq_store_round_trip_via_anf():
    """EqStore → to_anf_list → reconstruct coef vectors matches original."""
    from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore

    store = EqStore(comb_to_idx=_COMB_TO_IDX)
    f1 = XOR(AND(VAR(0), VAR(1)), VAR(2))
    f2 = XOR(VAR(0), CONST(1))
    store.insert_equation(f1)
    store.insert_equation(f2)

    original_matrix, c2i, i2c = to_coef_matrix(store)
    anf_list = to_anf_list(store)

    # Reconstruct matrix from ANF list
    reconstructed = np.zeros_like(original_matrix)
    for row, anf in enumerate(anf_list):
        for term in anf.terms:
            comb = tuple(sorted(term))
            reconstructed[row, c2i[comb]] = 1

    np.testing.assert_array_equal(reconstructed, original_matrix)
