"""Sanity tests for equation stores.

Tests use small systems with 2-4 variables to keep runtime low while covering
the core contracts:
- LUEqStore: independent equations increase rank, dependent ones don't
- SymbolicEqStore: accepts all inputs, accumulates without reduction
- EqStore: bag of coefficient vectors
- Cross-store linking propagates monomials
- LUEqStore.rank matches an independently computed GF(2) rank
"""
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, CONST, VAR, XOR, BooleanANF

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.FilteringEqStore import (
    FilteringEqStore,
)
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore2 import (
    GroebnerEqStore2,
)
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
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


# ── Store rank against an independent GF(2) rank ──────────────────────────────

@pytest.mark.parametrize("trial", range(10))
def test_store_rank_equals_the_rank_of_the_system_it_holds(trial):
    """`LUEqStore.rank` is the GF(2) rank of every equation inserted so far.

    The store maintains an LU-style factorization incrementally, accepting an
    equation only when it is independent of the ones already held.  The running
    rank must therefore match the rank of the whole batch computed offline --
    if it drifts, an attack will believe it has enough equations to solve when
    it does not.
    """
    rng = random.Random(600 + trial)

    def gf2_rank(rows):
        """Rank over GF(2) of rows given as integer bitmasks."""
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

    variables = rng.randint(2, 6)
    store, masks = LUEqStore(), []
    for _ in range(rng.randint(1, 10)):
        support = [i for i in range(variables) if rng.random() < 0.5]
        store.insert_equation(XOR(*[VAR(i) for i in support]) if support else CONST(0))
        masks.append(sum(1 << i for i in support))

    assert store.rank == gf2_rank(masks), (
        f"trial {trial}: store reports rank {store.rank} for a system of rank "
        f"{gf2_rank(masks)}"
    )


# ── FilteringEqStore contract ────────────────────────────────────────────────
# A filtering store converges toward a solution; a non-filtering one only
# accumulates. The type is the authority on which a store is, and every
# filtering store reports progress as num_determined, whatever it reduces with.

def test_filtering_flag_agrees_with_the_type():
    """`store.filtering` and isinstance(store, FilteringEqStore) are one question."""
    stores = [EqStore(), SymbolicEqStore(), LUEqStore(),
              GroebnerEqStore(None), GroebnerEqStore2()]
    for store in stores:
        assert store.filtering == isinstance(store, FilteringEqStore), (
            f"{type(store).__name__} reports filtering={store.filtering} but "
            f"isinstance says {isinstance(store, FilteringEqStore)}"
        )

def test_only_filtering_stores_report_progress():
    assert not hasattr(EqStore(), "num_determined")
    assert not hasattr(SymbolicEqStore(), "num_determined")
    assert LUEqStore().num_determined == 0
    assert GroebnerEqStore(None).num_determined == 0

def test_lu_progress_is_the_rank():
    """LU counts pivot columns, so progress moves only on an independent row."""
    store = LUEqStore()
    store.insert_equation(XOR(VAR(0), VAR(1)))
    store.insert_equation(XOR(VAR(1), VAR(2)))
    assert store.num_determined == store.rank == 2
    store.insert_equation(XOR(VAR(0), VAR(2)))  # the XOR of the first two
    assert store.num_determined == 2, "a dependent row must not register as progress"

@pytest.mark.parametrize("make_store", [lambda: GroebnerEqStore(None), GroebnerEqStore2],
                         ids=["GroebnerEqStore", "GroebnerEqStore2"])
def test_grobner_num_vars_counts_variables_ever_seen(make_store):
    """num_vars must not fall when reduction moves a variable to solved_vars.

    unknown_vars shrinks as variables are pinned, so num_vars is maintained at
    the one place the set grows. Solved variables are composed out before
    idxs_used() runs, so a variable reappearing in a later equation is not
    counted twice -- which the x0 in the fourth and sixth equations exercises.
    """
    store = make_store()
    equations = [
        XOR(VAR(0), VAR(1)),                        # new: 0, 1
        XOR(VAR(1), VAR(2), CONST(1)),              # new: 2
        XOR(VAR(0), CONST(1)),                      # pins 0, cascades to 1 and 2
        XOR(VAR(0), VAR(3)),                        # 0 already solved; 3 is new
        XOR(VAR(1), VAR(2), CONST(1)),              # nothing new
        XOR(AND(VAR(4), VAR(5)), VAR(0), CONST(1)), # new: 4, 5
    ]
    for equation in equations:
        store.enqueue_equation(equation)
        store.process_pending()
        ever_seen = store.unknown_vars | set(store.solved_vars)
        assert store.num_vars == len(ever_seen), (
            f"num_vars={store.num_vars} but {len(ever_seen)} variables have been seen"
        )

@pytest.mark.parametrize(
    "make_store",
    [LUEqStore, lambda: GroebnerEqStore(None), GroebnerEqStore2],
    ids=["LUEqStore", "GroebnerEqStore", "GroebnerEqStore2"],
)
def test_progress_reaches_num_vars_exactly_when_determined(make_store):
    """The completion invariant every consumer of num_determined relies on."""
    store = make_store()
    for equation in [XOR(VAR(0), VAR(1)), XOR(VAR(1), VAR(2), CONST(1)),
                     XOR(VAR(0), CONST(1)), XOR(VAR(0), VAR(3))]:
        store.queue_equation(equation)
        store.process_pending()
        assert (store.num_determined == store.num_vars) == store.is_determined, (
            f"{type(store).__name__}: num_determined={store.num_determined}, "
            f"num_vars={store.num_vars}, is_determined={store.is_determined}"
        )


# ── process_pending reports how much it consumed ─────────────────────────────
# BaseEqStore.process_pending returns the number of pending equations consumed,
# zero for a store with nothing deferred, so that a caller may accumulate it
# across stores without checking which kind it holds -- SplitGrob_Solver and the
# online inserter both do. It used to be declared `-> None` with a `pass` body
# while the Groebner stores returned an int, so `count += process_pending()`
# worked or raised TypeError depending on the store.

@pytest.mark.parametrize("make_store", [
    EqStore, SymbolicEqStore, LUEqStore,
    lambda: GroebnerEqStore(None), GroebnerEqStore2,
], ids=["EqStore", "SymbolicEqStore", "LUEqStore", "GroebnerEqStore", "GroebnerEqStore2"])
def test_process_pending_returns_an_int_for_every_store(make_store):
    store = make_store()
    store.queue_equation(XOR(VAR(0), VAR(1)))
    store.queue_equation(XOR(VAR(1), CONST(1)))
    consumed = store.process_pending()
    assert isinstance(consumed, int)
    total = 0
    total += store.process_pending()          # the accumulation callers rely on
    assert total >= 0

def test_deferred_stores_count_what_they_consume():
    """The Groebner stores defer reduction, so the count is the queue they drained."""
    for store in (GroebnerEqStore(None), GroebnerEqStore2()):
        store.queue_equation(XOR(VAR(0), VAR(1)))
        store.queue_equation(XOR(VAR(1), CONST(1)))
        assert store.process_pending() == 2
        assert store.process_pending() == 0, "a drained queue has nothing left to consume"
