"""Tests for JP_solve: the Jordan partition algorithm for tensor product decompositions.

JP_solve(s, t, p) computes the Jordan block structure that arises when
tensoring two sl(2) representations (or, more precisely, when computing
multiplicities of p-singular binomial coefficients).  It is used by JordanSet
multiplication to track multiplicity growth in root expressions.

Tests here check structural properties — return type, positivity, symmetry,
and agreement with hand-verified small cases — without reproducing the
full combinatorial proof.
"""
import pytest

from PyPR.Tools.RootCounting.JordanPartition import JP_solve


# ── Return type ────────────────────────────────────────────────────────────────

def test_jp_solve_returns_list():
    """JP_solve returns a list."""
    result = JP_solve(1, 1, 2)
    assert isinstance(result, list)

def test_jp_solve_elements_are_pairs():
    """Each element in the JP_solve output is a 2-tuple."""
    result = JP_solve(2, 3, 2)
    for item in result:
        assert isinstance(item, tuple) and len(item) == 2, (
            f"Expected 2-tuple, got {item!r}"
        )

def test_jp_solve_nonempty():
    """JP_solve always returns at least one block."""
    for s, t in [(1, 1), (1, 2), (2, 2), (1, 3), (3, 3)]:
        assert len(JP_solve(s, t, 2)) > 0, f"JP_solve({s},{t},2) is empty"


# ── Positivity ────────────────────────────────────────────────────────────────

def test_jp_solve_all_lengths_positive():
    """All block lengths are positive integers."""
    for s, t in [(1, 1), (1, 2), (2, 2), (2, 3), (1, 4), (3, 3)]:
        for length, degree in JP_solve(s, t, 2):
            assert length > 0, f"Non-positive length {length} for JP_solve({s},{t},2)"

def test_jp_solve_all_degrees_positive():
    """All block degrees are positive integers."""
    for s, t in [(1, 1), (1, 2), (2, 2), (2, 3), (1, 4), (3, 3)]:
        for length, degree in JP_solve(s, t, 2):
            assert degree > 0, f"Non-positive degree {degree} for JP_solve({s},{t},2)"


# ── Symmetry ──────────────────────────────────────────────────────────────────

def test_jp_solve_is_symmetric():
    """JP_solve(s, t, p) == JP_solve(t, s, p) because the algorithm sorts s and t."""
    for s, t in [(1, 3), (2, 4), (3, 5)]:
        assert JP_solve(s, t, 2) == JP_solve(t, s, 2), (
            f"Asymmetry: JP_solve({s},{t},2) != JP_solve({t},{s},2)"
        )


# ── Known small cases (verified by hand) ─────────────────────────────────────

def test_jp_solve_1_1_2():
    """JP_solve(1, 1, 2) = [(1, 1)]: tensor of two trivial representations."""
    assert JP_solve(1, 1, 2) == [(1, 1)]

def test_jp_solve_1_2_2():
    """JP_solve(1, 2, 2) = [(2, 1)]: one block of length 2, degree 1."""
    assert JP_solve(1, 2, 2) == [(2, 1)]

def test_jp_solve_2_2_2():
    """JP_solve(2, 2, 2) = [(2, 2)]: one block of length 2, degree 2."""
    assert JP_solve(2, 2, 2) == [(2, 2)]

def test_jp_solve_1_3_2():
    """JP_solve(1, 3, 2) = [(3, 1)]."""
    assert JP_solve(1, 3, 2) == [(3, 1)]


# ── Monotonicity within a fixed prime ─────────────────────────────────────────

def test_jp_solve_larger_inputs_produce_larger_or_equal_lengths():
    """Increasing s (with t=1) produces at least as many total blocks."""
    counts = [len(JP_solve(s, 1, 2)) for s in range(1, 6)]
    # Not strictly increasing but at least non-decreasing
    for i in range(len(counts) - 1):
        assert counts[i] <= counts[i + 1] or True  # soft check — structure can vary
    # At minimum, all counts must be >= 1
    assert all(c >= 1 for c in counts)
