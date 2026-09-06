"""Sanity tests for equation stores.

Tests use small systems with 2-4 variables to keep runtime low while covering
the core contracts:
- LUEqStore: independent equations increase rank, dependent ones don't
- SymbolicEqStore: accepts all inputs, accumulates without reduction
- EqStore: bag of coefficient vectors
- Cross-store linking propagates monomials
"""
import numpy as np
from PyPR.BooleanLogic import BooleanANF, VAR, XOR, AND, CONST
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import SymbolicEqStore


# ── Empty store ───────────────────────────────────────────────────────────────

def test_empty_store_rank_zero():
    assert LUEqStore().rank == 0


# ── Single insertion ──────────────────────────────────────────────────────────

def test_insert_single_equation_increases_rank():
    store = LUEqStore()
    store.insert_equation(VAR(0))
    assert store.rank == 1

def test_insert_returns_true_for_independent():
    store = LUEqStore()
    result = store.insert_equation(VAR(0))
    assert result is True


# ── Independence detection ────────────────────────────────────────────────────

def test_duplicate_equation_does_not_increase_rank():
    """Inserting the same equation twice should keep rank at 1."""
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(0))
    assert store.rank == 1

def test_duplicate_returns_false():
    store = LUEqStore()
    store.insert_equation(VAR(0))
    result = store.insert_equation(VAR(0))
    assert result is False

def test_xor_combination_is_dependent():
    """x0 XOR x1 is linearly dependent once x0 and x1 are in the store."""
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    result = store.insert_equation(XOR(VAR(0), VAR(1)))
    assert result is False
    assert store.rank == 2

def test_xor_combination_with_constant_is_dependent():
    """x0 XOR x1 XOR 1 is dependent once x0, x1, and CONST(1) are in the store."""
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    store.insert_equation(CONST(1))
    result = store.insert_equation(XOR(VAR(0), VAR(1), CONST(1)))
    assert result is False
    assert store.rank == 3


# ── Multiple independent equations ────────────────────────────────────────────

def test_three_independent_vars_rank_three():
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    store.insert_equation(VAR(2))
    assert store.rank == 3

def test_rank_increments_monotonically():
    """Rank should be non-decreasing as equations are inserted."""
    store = LUEqStore()
    prev_rank = 0
    for i in range(4):
        store.insert_equation(VAR(i))
        assert store.rank >= prev_rank
        prev_rank = store.rank
    assert store.rank == 4

def test_nonlinear_equation_is_independent():
    """AND(x0, x1) is linearly independent from x0 and x1 individually."""
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    result = store.insert_equation(AND(VAR(0), VAR(1)))
    assert result is True
    assert store.rank == 3

def test_and_then_xor_independent():
    """x0*x1 XOR x0 is independent of x0 and x1 alone."""
    store = LUEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    result = store.insert_equation(XOR(AND(VAR(0), VAR(1)), VAR(0)))
    assert result is True
    assert store.rank == 3


# ── SymbolicEqStore ──────────────────────────────────────────────────────────

_COMB_TO_IDX = {
    (): 0,
    (0,): 1, (1,): 2, (2,): 3,
    (0, 1): 4, (0, 2): 5, (1, 2): 6,
}

def test_symbolic_empty():
    store = SymbolicEqStore()
    assert store.num_eqs == 0
    assert store.equations == []


def test_symbolic_insert_boolean_function_dynamic():
    """Dynamic mode: discover monomials on insertion."""
    store = SymbolicEqStore()
    store.insert_equation(XOR(VAR(0), VAR(1)))
    assert store.num_eqs == 1
    assert len(store.equations) == 1
    assert isinstance(store.equations[0], BooleanANF)
    assert store.equations[0].terms == frozenset({frozenset({0}), frozenset({1})})


def test_symbolic_insert_multiple():
    store = SymbolicEqStore()
    store.insert_equation(VAR(0))
    store.insert_equation(VAR(1))
    store.insert_equation(XOR(AND(VAR(0), VAR(1)), VAR(2)))
    assert store.num_eqs == 3
    assert len(store.equations) == 3


def test_symbolic_accepts_duplicates():
    """Bags don't reject — duplicate equations are kept."""
    store = SymbolicEqStore()
    store.insert_equation(VAR(0))
    r = store.insert_equation(VAR(0))
    assert r is True
    assert store.num_eqs == 2


def test_symbolic_insert_returns_true():
    store = SymbolicEqStore()
    assert store.insert_equation(VAR(0)) is True


def test_symbolic_static_ndarray():
    """Static mode: accept ndarray input."""
    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    vec = np.array([1, 1, 0, 0, 0, 0, 0], dtype=np.uint8)
    store.insert_equation(vec)
    assert store.num_eqs == 1
    assert store.equations[0].terms == frozenset({frozenset(), frozenset({0})})


def test_symbolic_static_boolean_function():
    """Static mode: accept BooleanFunction input."""
    store = SymbolicEqStore(comb_to_idx=_COMB_TO_IDX)
    store.insert_equation(XOR(VAR(0), VAR(1)))
    assert store.num_eqs == 1
    assert store.equations[0].terms == frozenset({frozenset({0}), frozenset({1})})


def test_symbolic_dynamic_rejects_ndarray():
    """Dynamic mode cannot accept ndarray (no index map)."""
    store = SymbolicEqStore()
    import pytest
    with pytest.raises(ValueError, match="Cannot insert an ndarray"):
        store.insert_equation(np.array([1, 0, 1], dtype=np.uint8))


def test_symbolic_discovers_monomials_dynamic():
    """Dynamic mode should build index maps from inserted equations."""
    store = SymbolicEqStore()
    store.insert_equation(XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1)))
    assert (0, 1) in store.comb_to_idx
    assert (2,) in store.comb_to_idx
    assert () in store.comb_to_idx


def test_symbolic_equation_ids():
    """Identifiers are stored per equation."""
    store = SymbolicEqStore()
    store.insert_equation(VAR(0), identifier="t=0")
    store.insert_equation(VAR(1), identifier="t=1")
    assert store.equation_ids[0] == "t=0"
    assert store.equation_ids[1] == "t=1"


def test_symbolic_nonlinear():
    """Nonlinear terms are preserved as BooleanANF."""
    store = SymbolicEqStore()
    store.insert_equation(AND(VAR(0), VAR(1)))
    assert store.equations[0].terms == frozenset({frozenset({0, 1})})


# ── Cross-store linking ──────────────────────────────────────────────────────

def test_link_symbolic_to_eq_store():
    """Monomials discovered in SymbolicEqStore propagate to a linked EqStore."""
    sym = SymbolicEqStore()
    eq = EqStore()
    sym.link(eq)

    sym.insert_equation(XOR(VAR(0), VAR(1)))
    assert (0,) in eq.comb_to_idx
    assert (1,) in eq.comb_to_idx


def test_link_eq_store_to_symbolic():
    """Monomials discovered in EqStore propagate to a linked SymbolicEqStore."""
    eq = EqStore()
    sym = SymbolicEqStore()
    eq.link(sym)

    eq.insert_equation(XOR(VAR(0), VAR(1)))
    assert (0,) in sym.comb_to_idx
    assert (1,) in sym.comb_to_idx


def test_link_symbolic_to_lu_store():
    """Monomials discovered in SymbolicEqStore propagate to a linked LUEqStore."""
    sym = SymbolicEqStore()
    lu = LUEqStore()
    sym.link(lu)

    sym.insert_equation(XOR(AND(VAR(0), VAR(1)), VAR(2)))
    assert (0, 1) in lu.comb_to_idx
    assert (2,) in lu.comb_to_idx
