"""Tests for the Adapters layer: equation representation and online insertion.

`equation_repr` converts a single equation between the three forms the rest of
the pipeline uses -- a BooleanFunction DAG, a coefficient vector over a monomial
index map, and a BooleanANF. The property that matters is that the conversions
agree: whichever way an equation is carried, it must denote the same function.

`online_insertion` is a factory rather than a converter. It inspects a store
once and returns closures specialised for that store's shape, so the hot loop
carries no type checks. Four variants come out of two independent axes -- is the
store indexed, and can it stop early -- and each has to insert correctly.
"""
import itertools

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, CONST, VAR, XOR

from PyPR.Cryptanalysis.Components.Adapters.equation_repr import (
    boolean_function_to_coef_vector,
    coef_vector_to_anf,
    coef_vector_to_boolean_function,
    extract_monomials,
)
from PyPR.Cryptanalysis.Components.Adapters.online_insertion import make_online_inserter
from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore
from PyPR.Cryptanalysis.Components.EquationStores.IndexedEqStore import IndexedEqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import SymbolicEqStore

# monomials over three variables: the constant, each single, and one product
_COMB_TO_IDX = {(): 0, (0,): 1, (1,): 2, (2,): 3, (0, 1): 4}
_IDX_TO_COMB = {v: k for k, v in _COMB_TO_IDX.items()}


# ── equation_repr: the three representations must agree ──────────────────────

def test_extract_monomials_reads_an_anf_term_by_term():
    """An ANF's terms are products of variables; a constant 1 is the empty tuple."""
    equation = XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1))
    monomials = extract_monomials(equation)

    assert sorted(monomials) == [(), (0, 1), (2,)]

def test_extract_monomials_drops_a_zero_constant():
    """CONST(0) contributes nothing to an ANF and must not become a term."""
    assert extract_monomials(XOR(VAR(0), CONST(0))) == [(0,)]

@pytest.mark.parametrize("equation", [
    VAR(0),
    XOR(VAR(0), VAR(1)),
    XOR(AND(VAR(0), VAR(1)), VAR(2)),
    XOR(AND(VAR(0), VAR(1)), CONST(1)),
    CONST(1),
], ids=["single", "xor", "with-product", "with-constant", "constant"])
def test_coefficient_vector_round_trips_through_every_representation(equation):
    """DAG -> coefficient vector -> ANF and -> DAG must all denote one function."""
    vector = boolean_function_to_coef_vector(equation, _COMB_TO_IDX, len(_COMB_TO_IDX))
    as_fn = coef_vector_to_boolean_function(vector, _IDX_TO_COMB)
    # BooleanANF is a term set rather than something evaluable, so it is checked
    # through the DAG it converts back to
    from_anf = coef_vector_to_anf(vector, _IDX_TO_COMB).to_BooleanFunction()

    for bits in itertools.product([0, 1], repeat=3):
        state = list(bits)
        expected = equation.eval(state)
        assert as_fn.eval(state) == expected, f"DAG round trip differs at {state}"
        assert from_anf.eval(state) == expected, f"ANF round trip differs at {state}"

def test_the_coefficient_vector_marks_exactly_the_terms_present():
    equation = XOR(AND(VAR(0), VAR(1)), VAR(2))
    vector = boolean_function_to_coef_vector(equation, _COMB_TO_IDX, len(_COMB_TO_IDX))

    assert list(vector) == [0, 0, 0, 1, 1], "expected only the (2,) and (0,1) columns set"


# ── online_insertion: four variants, two axes ────────────────────────────────
# `is_indexed` decides whether coefficient vectors go in directly or are
# converted to BooleanFunctions first; `eager and filtering` decides whether the
# loop may stop early once the store is determined. Every store below picks a
# different corner, and all of them must accept the equations they are handed.

@pytest.mark.parametrize(("make_store", "indexed", "can_stop"), [
    (lambda: LUEqStore(_COMB_TO_IDX, consistent=True), True, True),
    (lambda: EqStore(_COMB_TO_IDX), True, False),
    (lambda: SymbolicEqStore(_COMB_TO_IDX), True, False),
    (lambda: GroebnerEqStore(simplify_mode=None), False, False),
], ids=["LUEqStore", "EqStore", "SymbolicEqStore", "GroebnerEqStore"])
def test_online_inserter_accepts_equations_for_each_store_shape(make_store, indexed, can_stop):
    store = make_store()
    # the two axes the factory dispatches on, pinned so a store that changes
    # shape moves the case rather than silently retargeting it
    assert isinstance(store, IndexedEqStore) == indexed
    assert (store.eager and store.filtering) == can_stop, (
        "this store no longer sits in the corner this case was written for"
    )

    insert_eq, finalize = make_online_inserter(
        store, _IDX_TO_COMB, total_eqs=3, num_vars=len(_COMB_TO_IDX),
        verbose=False,
    )
    # x0 = 1, x1 = 0, x2 = 1, as coefficient vectors over _COMB_TO_IDX
    vectors = [
        np.array([1, 1, 0, 0, 0], dtype=np.uint8),
        np.array([0, 0, 1, 0, 0], dtype=np.uint8),
        np.array([1, 0, 0, 1, 0], dtype=np.uint8),
    ]
    for eq_idx, vector in enumerate(vectors):
        stop = insert_eq(vector, eq_idx)
        assert isinstance(stop, bool), "insert_fn must report whether to stop early"
        if stop:
            break
    finalize()

    if not store.eager:
        store.process_pending()
    assert store.num_eqs > 0, "no equation reached the store"

def test_online_inserter_never_stops_early_for_an_accumulator():
    """A pure accumulator has no notion of being determined, so it never stops."""
    store = EqStore(_COMB_TO_IDX)
    insert_eq, finalize = make_online_inserter(
        store, _IDX_TO_COMB, total_eqs=2, num_vars=len(_COMB_TO_IDX), verbose=False,
    )
    stops = [insert_eq(np.array([1, 1, 0, 0, 0], dtype=np.uint8), 0),
             insert_eq(np.array([0, 0, 1, 0, 0], dtype=np.uint8), 1)]
    finalize()

    assert stops == [False, False]

def test_online_inserter_stops_exactly_when_a_filtering_store_is_determined():
    """An eager filtering store reports stop on the insertion that determines it.

    Over the monomials 1, x0, x1, x2 a consistent LU store starts with the
    constant row already pinned, so three independent equations bring its rank
    to the full four columns -- and not one insertion earlier.
    """
    comb_to_idx = {(): 0, (0,): 1, (1,): 2, (2,): 3}
    idx_to_comb = {v: k for k, v in comb_to_idx.items()}
    store = LUEqStore(comb_to_idx, consistent=True)
    insert_eq, finalize = make_online_inserter(
        store, idx_to_comb, total_eqs=3, num_vars=len(comb_to_idx), verbose=False,
    )
    # x0 = 1, x1 = 0, x2 = 1
    vectors = [[1, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 1]]
    stops = [insert_eq(np.array(v, dtype=np.uint8), k) for k, v in enumerate(vectors)]
    finalize()

    assert stops == [False, False, True]
    assert store.is_determined
