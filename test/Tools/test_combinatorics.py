"""Tests for combinatorics utilities: choose, binsum, powerset.

These are building blocks for root counting and monomial profile upper-bound
computations.  Each function has a simple closed-form definition that can be
checked against the standard library or by hand.
"""
import math

from PyPR.Tools.RootCounting.Combinatorics import choose, binsum, powerset


# ── choose ────────────────────────────────────────────────────────────────────

def test_choose_matches_math_comb_small():
    """choose(n, k) equals math.comb(n, k) for small values."""
    for n in range(10):
        for k in range(n + 1):
            assert choose(n, k) == math.comb(n, k), (
                f"choose({n},{k}) = {choose(n,k)}, expected {math.comb(n,k)}"
            )

def test_choose_zero_k_is_one():
    """C(n, 0) == 1 for all n."""
    for n in range(15):
        assert choose(n, 0) == 1

def test_choose_k_equals_n_is_one():
    """C(n, n) == 1 for all n."""
    for n in range(15):
        assert choose(n, n) == 1

def test_choose_symmetry():
    """C(n, k) == C(n, n-k) (Pascal symmetry)."""
    for n in range(12):
        for k in range(n + 1):
            assert choose(n, k) == choose(n, n - k), (
                f"Symmetry violated at choose({n},{k})"
            )

def test_choose_pascals_rule():
    """C(n, k) == C(n-1, k-1) + C(n-1, k) (Pascal's rule)."""
    for n in range(2, 12):
        for k in range(1, n):
            assert choose(n, k) == choose(n - 1, k - 1) + choose(n - 1, k), (
                f"Pascal's rule violated at choose({n},{k})"
            )

def test_choose_large():
    """choose handles moderate values correctly."""
    assert choose(20, 10) == math.comb(20, 10)
    assert choose(30, 5) == math.comb(30, 5)


# ── binsum ────────────────────────────────────────────────────────────────────

def test_binsum_is_sum_of_binomials():
    """binsum(n, d) == sum_{k=1}^{d} C(n, k)."""
    for n in range(10):
        for d in range(1, n + 1):
            expected = sum(math.comb(n, k) for k in range(1, d + 1))
            assert binsum(n, d) == expected, (
                f"binsum({n},{d}) = {binsum(n,d)}, expected {expected}"
            )

def test_binsum_zero_degree_is_zero():
    """binsum(n, 0) == 0 (empty sum)."""
    for n in range(10):
        assert binsum(n, 0) == 0

def test_binsum_full_degree_is_two_to_n_minus_one():
    """binsum(n, n) == 2^n - 1 (sum of all proper subsets)."""
    for n in range(1, 10):
        assert binsum(n, n) == 2**n - 1, (
            f"binsum({n},{n}) = {binsum(n,n)}, expected {2**n - 1}"
        )


# ── powerset ──────────────────────────────────────────────────────────────────

def test_powerset_size_is_two_to_n():
    """powerset of an n-element list has 2^n elements."""
    for n in range(6):
        ls = list(range(n))
        count = sum(1 for _ in powerset(ls))
        assert count == 2**n, f"powerset size for n={n}: {count} != {2**n}"

def test_powerset_contains_empty_set():
    """The empty tuple () is always in powerset."""
    ps = list(powerset([1, 2, 3]))
    assert () in ps

def test_powerset_contains_full_set():
    """The full set is always in powerset."""
    ls = [1, 2, 3]
    ps = list(powerset(ls))
    assert tuple(ls) in ps

def test_powerset_all_subsets_are_subsets():
    """Every element of powerset is a subset of the original list."""
    ls = [10, 20, 30]
    ls_set = set(ls)
    for subset in powerset(ls):
        assert set(subset).issubset(ls_set), f"{subset} is not a subset of {ls}"

def test_powerset_empty_list():
    """powerset of an empty list yields only the empty tuple."""
    result = list(powerset([]))
    assert result == [()]
