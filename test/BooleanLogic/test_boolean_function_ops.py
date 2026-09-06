"""Tests for BooleanFunction DAG operations: index queries, manipulation, and composition.

These tests cover the structural and semantic operations on BooleanFunction DAGs
that are not exercised by the gate-evaluation or ANF tests: variable index
tracking, deep copy semantics, subfunction detection, argument mutation, index
remapping, composition, and the agreement between eval() and eval_ANF().
"""
from itertools import product as iter_product

import pytest

from PyPR.BooleanLogic import AND, OR, XOR, NOT, VAR, CONST
from PyPR.BooleanLogic.BooleanANF import BooleanANF


# ── idxs_used ────────────────────────────────────────────────────────────────

def test_idxs_used_const():
    """CONST nodes reference no variables."""
    assert CONST(0).idxs_used() == set()

def test_idxs_used_var():
    """VAR(i) references only index i."""
    assert VAR(5).idxs_used() == {5}

def test_idxs_used_composite():
    """A composite function collects indices from all reachable VAR leaves."""
    fn = XOR(VAR(0), AND(VAR(2), VAR(4)))
    assert fn.idxs_used() == {0, 2, 4}

def test_idxs_used_shared_node():
    """Shared DAG nodes are not double-counted."""
    shared = VAR(3)
    fn = XOR(shared, AND(shared, VAR(1)))
    assert fn.idxs_used() == {1, 3}


# ── max_idx ───────────────────────────────────────────────────────────────────

def test_max_idx_const_is_negative_one():
    """CONST has max_idx() == -1 (no variable indices)."""
    assert CONST(1).max_idx() == -1

def test_max_idx_var():
    """VAR(i).max_idx() == i."""
    assert VAR(7).max_idx() == 7

def test_max_idx_composite():
    """Composite function returns the highest index among all VAR nodes."""
    fn = AND(VAR(1), VAR(5), VAR(3))
    assert fn.max_idx() == 5


# ── is_leaf ───────────────────────────────────────────────────────────────────

def test_const_is_leaf():
    assert CONST(0).is_leaf() is True

def test_var_is_leaf():
    assert VAR(0).is_leaf() is True

def test_gate_is_not_leaf():
    assert XOR(VAR(0), VAR(1)).is_leaf() is False


# ── condense_idxs ─────────────────────────────────────────────────────────────

def test_condense_idxs_produces_consecutive():
    """condense_idxs() remaps sparse indices to {0, 1, ..., k-1}."""
    fn = XOR(VAR(0), VAR(5), VAR(10))
    condensed = fn.condense_idxs()
    assert condensed.idxs_used() == {0, 1, 2}

def test_condense_idxs_preserves_evaluation():
    """condense_idxs() produces a function that evaluates identically on mapped inputs."""
    fn = XOR(VAR(0), VAR(3))
    condensed = fn.condense_idxs()
    # Original: x0 XOR x3.  After condensing: x0 XOR x1.
    for a, b in iter_product([0, 1], repeat=2):
        orig_val = fn.eval([a, 0, 0, b])
        cond_val = condensed.eval([a, b])
        assert orig_val == cond_val, f"Mismatch at ({a},{b})"

def test_condense_idxs_already_consecutive_is_noop():
    """condense_idxs() on {0,1,2} leaves indices unchanged."""
    fn = XOR(VAR(0), VAR(1), VAR(2))
    condensed = fn.condense_idxs()
    assert condensed.idxs_used() == {0, 1, 2}


# ── compose ───────────────────────────────────────────────────────────────────

def test_compose_replaces_variable():
    """compose({i: g}) substitutes g for VAR(i)."""
    fn = XOR(VAR(0), VAR(1))
    result = fn.compose({0: CONST(1)})
    # (1 XOR x1)
    assert result.eval([0, 0]) == 1
    assert result.eval([0, 1]) == 0

def test_compose_leaves_unmapped_variables():
    """compose() leaves VAR nodes not in the map unchanged."""
    fn = XOR(VAR(0), VAR(1))
    result = fn.compose({0: CONST(0)})
    # VAR(1) should still be live
    assert 1 in result.idxs_used()

def test_compose_with_another_function():
    """Composing with a non-trivial function correctly chains evaluation.

    outer = AND(VAR(0), VAR(1)); VAR(0) is replaced by XOR(VAR(2), VAR(3)).
    Result: AND(XOR(VAR(2), VAR(3)), VAR(1)) — eval uses bits[1], bits[2], bits[3].
    """
    outer = AND(VAR(0), VAR(1))
    inner = XOR(VAR(2), VAR(3))
    composed = outer.compose({0: inner})
    for bits in iter_product([0, 1], repeat=4):
        # composed = (bits[2] XOR bits[3]) AND bits[1]
        expected = (bits[2] ^ bits[3]) & bits[1]
        assert composed.eval(list(bits)) == expected, f"Mismatch at {bits}"


# ── remap_indices ──────────────────────────────────────────────────────────────

def test_remap_indices_renumbers_correctly():
    """remap_indices({i: j}) shifts variable i to variable j."""
    fn = VAR(0)
    remapped = fn.remap_indices({0: 4})
    assert 4 in remapped.idxs_used()
    assert 0 not in remapped.idxs_used()

def test_remap_indices_composite():
    """remap_indices remaps all matching VAR leaves in a composite function."""
    fn = XOR(VAR(0), VAR(1))
    remapped = fn.remap_indices({0: 2, 1: 3})
    assert remapped.idxs_used() == {2, 3}


# ── eval vs eval_ANF agreement ────────────────────────────────────────────────

def _check_eval_agreement(fn, n_vars):
    """Assert eval and eval_ANF agree on all 2^n_vars inputs."""
    for bits in iter_product([0, 1], repeat=n_vars):
        inp = list(bits)
        assert fn.eval(inp) == fn.eval_ANF(inp), (
            f"eval/eval_ANF mismatch at {inp}"
        )

def test_eval_eval_anf_agree_xor():
    """eval and eval_ANF agree for a 3-input XOR."""
    _check_eval_agreement(XOR(VAR(0), VAR(1), VAR(2)), 3)

def test_eval_eval_anf_agree_and():
    """eval and eval_ANF agree for a 2-input AND."""
    _check_eval_agreement(AND(VAR(0), VAR(1)), 2)

def test_eval_eval_anf_agree_or():
    """eval and eval_ANF agree for a 2-input OR."""
    _check_eval_agreement(OR(VAR(0), VAR(1)), 2)

def test_eval_eval_anf_agree_not():
    """eval and eval_ANF agree for NOT."""
    _check_eval_agreement(NOT(VAR(0)), 1)

def test_eval_eval_anf_agree_composite():
    """eval and eval_ANF agree for a nontrivial mixed-gate function."""
    fn = XOR(AND(VAR(0), VAR(1)), OR(VAR(1), VAR(2)))
    _check_eval_agreement(fn, 3)


# ── copy ──────────────────────────────────────────────────────────────────────

def test_copy_is_functionally_equivalent():
    """A copied function evaluates identically to the original."""
    fn = XOR(AND(VAR(0), VAR(1)), VAR(2))
    fn_copy = fn.__copy__()
    for bits in iter_product([0, 1], repeat=3):
        assert fn.eval(list(bits)) == fn_copy.eval(list(bits))

def test_copy_is_structurally_separate():
    """Modifying the copy's args does not affect the original."""
    fn = XOR(VAR(0), VAR(1))
    fn_copy = fn.__copy__()
    fn_copy.add_arguments(VAR(2))
    assert len(fn.args) == 2, "Original should still have 2 args"
    assert len(fn_copy.args) == 3


# ── subfunctions ──────────────────────────────────────────────────────────────

def test_subfunctions_empty_for_tree():
    """A pure tree (no shared nodes) has no subfunctions."""
    fn = XOR(AND(VAR(0), VAR(1)), AND(VAR(2), VAR(3)))
    assert fn.subfunctions() == []

def test_subfunctions_detects_shared_node():
    """A node referenced by two parents appears in subfunctions()."""
    shared = AND(VAR(0), VAR(1))
    fn = XOR(shared, shared)
    subs = fn.subfunctions()
    assert shared in subs

def test_subfunctions_topologically_sorted():
    """subfunctions() are in topological (bottom-up) order.

    Build a DAG where 'a' is shared by both 'b' and the top-level function,
    and 'b' is itself shared by the top-level function twice.  Since 'b'
    depends on 'a', 'a' must appear first in the subfunction list.
    """
    a = AND(VAR(0), VAR(1))         # shared: appears in b and top-level fn
    b = XOR(a, VAR(2))              # shared: appears twice in top-level fn
    fn = OR(b, AND(a, b))           # both a and b referenced more than once
    subs = fn.subfunctions()
    assert a in subs, "shared node 'a' should be a subfunction"
    assert b in subs, "shared node 'b' should be a subfunction"
    assert subs.index(a) < subs.index(b), (
        "'a' must appear before 'b' (a is a dependency of b)"
    )


# ── inputs ────────────────────────────────────────────────────────────────────

def test_inputs_returns_all_leaves():
    """inputs() returns all leaf nodes (VAR and CONST) in DFS order."""
    v0, v1 = VAR(0), VAR(1)
    fn = AND(v0, v1)
    leaves = fn.inputs()
    assert v0 in leaves
    assert v1 in leaves

def test_inputs_for_const():
    """CONST node returns itself as the only input."""
    c = CONST(1)
    fn = XOR(VAR(0), c)
    leaves = fn.inputs()
    assert c in leaves


# ── add_arguments / remove_arguments ─────────────────────────────────────────

def test_add_arguments_increases_arity():
    """add_arguments appends new children."""
    fn = XOR(VAR(0))
    fn.add_arguments(VAR(1), VAR(2))
    assert len(fn.args) == 3

def test_add_arguments_respects_arg_limit():
    """add_arguments raises ValueError when arg_limit would be exceeded."""
    fn = NOT(VAR(0))  # NOT has arg_limit=1
    with pytest.raises(ValueError):
        fn.add_arguments(VAR(1))

def test_remove_arguments_specific():
    """remove_arguments with arguments removes only those children."""
    v0, v1, v2 = VAR(0), VAR(1), VAR(2)
    fn = XOR(v0, v1, v2)
    fn.remove_arguments(v1)
    assert v1 not in fn.args
    assert v0 in fn.args and v2 in fn.args

def test_remove_arguments_none_clears_all():
    """remove_arguments() with no args clears the entire argument list."""
    fn = XOR(VAR(0), VAR(1), VAR(2))
    fn.remove_arguments()
    assert fn.args == tuple()

def test_remove_arguments_missing_is_silent():
    """Removing an argument not present does not raise an error."""
    v0, v1 = VAR(0), VAR(1)
    fn = XOR(v0)
    fn.remove_arguments(v1)  # v1 not in fn.args — should not raise
    assert len(fn.args) == 1


# ── binarize ──────────────────────────────────────────────────────────────────

def test_binarize_preserves_xor_semantics():
    """binarize() produces a function semantically identical to the original XOR."""
    fn = XOR(VAR(0), VAR(1), VAR(2), VAR(3))
    binary = fn.binarize()
    for bits in iter_product([0, 1], repeat=4):
        assert fn.eval(list(bits)) == binary.eval(list(bits)), (
            f"binarize mismatch at {bits}"
        )

def test_binarize_produces_binary_tree():
    """After binarize(), every gate node has at most 2 arguments."""
    fn = XOR(VAR(0), VAR(1), VAR(2), VAR(3))
    binary = fn.binarize()

    def max_arity(node, visited=None):
        if visited is None:
            visited = set()
        if id(node) in visited or node.is_leaf():
            return 0
        visited.add(id(node))
        return max(len(node.args), *(max_arity(c, visited) for c in node.args))

    assert max_arity(binary) <= 2
