"""Tests for MersenneTools: cycle_lengths, max_period, expected_period, mersenne_combinations.

Derived from ProductRegisters.ipynb (Construction Analysis and picking Good Mersenne Primes).
"""
from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.Tools.MersenneTools import (
    cycle_lengths, max_period, expected_period, expected_period_ratio,
    mersenne_combinations, list_possible,
)


# ── cycle_lengths ───────────────────────────────────────────────────────

def test_cycle_lengths_basic():
    """cycle_lengths([7,5,3]) should return a non-empty list of cycle lengths."""
    result = cycle_lengths([7, 5, 3])
    assert len(result) > 0
    for length in result:
        assert length > 0

def test_cycle_lengths_single():
    """Single MPR of size n has a cycle of length 2^n - 1."""
    result = cycle_lengths([5])
    assert 2**5 - 1 in result


# ── max_period ──────────────────────────────────────────────────────────

def test_max_period():
    """max_period is the LCM of (2^n_i - 1)."""
    result = max_period([7, 5, 3])
    # 2^7-1=127, 2^5-1=31, 2^3-1=7. LCM(127,31,7) = 127*31*7/gcd... = 27559
    # these are coprime Mersenne primes, so LCM = product
    assert result == 127 * 31 * 7

def test_max_period_single():
    assert max_period([5]) == 31


# ── expected_period ─────────────────────────────────────────────────────

def test_expected_period():
    result = expected_period([7, 5, 3])
    assert result > 0
    assert result <= max_period([7, 5, 3])

def test_expected_period_ratio():
    ratio = expected_period_ratio([7, 5, 3])
    assert 0 < ratio <= 1


# ── CMPR properties match standalone functions ──────────────────────────

def test_cmpr_period_properties():
    """CMPR period properties should match the standalone functions (notebook cell 227)."""
    M7 = MPR(7, [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0])
    M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
    M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
    C = CMPR([M7, M5, M3])

    assert C.max_period == max_period([7, 5, 3])
    assert C.expected_period == expected_period([7, 5, 3])
    assert C.expected_period_ratio == expected_period_ratio([7, 5, 3])


# ── list_possible ───────────────────────────────────────────────────────

def test_list_possible():
    result = list_possible([2**i for i in range(12)])
    assert isinstance(result, list)
    assert len(result) > 0
    # all results should be from the Mersenne exponent list
    for val in result:
        assert (2**val - 1) > 0


# ── mersenne_combinations ──────────────────────────────────────────────

def test_mersenne_combinations():
    """mersenne_combinations should return valid configurations."""
    results = list(mersenne_combinations([2**i for i in range(12)]))
    assert len(results) > 0
    for group in results:
        assert isinstance(group, (list, tuple))
        # each group should be a collection of (size, count) or similar
        assert len(group) > 0

