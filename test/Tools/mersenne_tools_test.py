"""Tests for MersenneTools: cycle_lengths, max_period, expected_period, mersenne_combinations.

Derived from ProductRegisters.ipynb (Construction Analysis and picking Good Mersenne Primes).

`expected_period` and `expected_period_ratio` are closed forms "proved in the
paper"; `expected_period_brute_force` and `epr_brute_force` are the
definitions they replace, evaluated by enumerating all 2^k cycle lengths.
Testing the fast form against its own definition is the point -- the closed
form is an optimization, and an algebra slip in it would not show up anywhere
else.
"""
import math
from functools import reduce

import pytest

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.Tools.MersenneTools import (
    cycle_lengths, max_period, expected_period, expected_period_ratio,
    mersenne_combinations, list_possible,
    epr_brute_force, expected_period_brute_force,
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


# ── Closed forms against their brute-force definitions ────────────────────────

BLOCK_SIZES = [[3, 4], [3, 5], [4, 5], [3, 4, 5], [2, 3, 5], [3, 5, 7], [2, 3, 5, 7]]


PRIMITIVE = {
    2: [1, 1, 1],
    3: [1, 1, 0, 1],
    4: [1, 1, 0, 0, 1],
    5: [1, 0, 1, 0, 0, 1],
    7: [1, 1, 0, 0, 0, 0, 0, 1],
}


@pytest.mark.parametrize("sizes", BLOCK_SIZES, ids=str)
def test_closed_form_expected_period_matches_its_brute_force_definition(sizes):
    """The product formula equals sum(c^2)/sum(c) over all cycle lengths c.

    The definition weights each cycle by the probability of landing on it
    (proportional to its length), giving E[period] = sum(c^2)/sum(c).  The
    closed form prod((2^s - 1)^2 + 1) / prod(2^s) is the factorization of that
    sum across independent blocks.
    """
    assert math.isclose(
        expected_period(sizes), expected_period_brute_force(sizes), rel_tol=1e-12
    ), f"{sizes}: {expected_period(sizes)} != {expected_period_brute_force(sizes)}"


@pytest.mark.parametrize("sizes", BLOCK_SIZES, ids=str)
def test_closed_form_expected_period_ratio_matches_its_brute_force_definition(sizes):
    """The ratio formula equals sum((c/total)^2) over all cycle lengths c."""
    assert math.isclose(
        expected_period_ratio(sizes), epr_brute_force(sizes), rel_tol=1e-12
    ), f"{sizes}: {expected_period_ratio(sizes)} != {epr_brute_force(sizes)}"


@pytest.mark.parametrize("sizes", BLOCK_SIZES, ids=str)
def test_expected_period_ratio_is_the_expected_period_over_the_state_count(sizes):
    """ratio = expected_period / prod(2^s), the full state-space size.

    Both brute-force forms divide by `total = sum(cycle_lengths)`, and that sum
    telescopes to prod(1 + (2^s - 1)) = prod(2^s) -- the number of states of the
    whole register.  So the ratio is normalized against the state count, which
    is strictly larger than `max_period = lcm(2^s - 1)`.
    """
    state_count = reduce(lambda a, b: a * b, [2 ** s for s in sizes])
    assert math.isclose(
        expected_period_ratio(sizes) * state_count, expected_period(sizes), rel_tol=1e-9
    ), f"{sizes}: ratio is not expected_period / {state_count}"


@pytest.mark.parametrize("sizes", BLOCK_SIZES, ids=str)
def test_cycle_lengths_enumerate_every_subset_of_blocks(sizes):
    """There are 2^k cycle lengths, summing to prod(2^s), with max = max_period.

    A cycle length is the product of (2^s - 1) over some subset of blocks --
    the subset whose registers are in a nonzero state.  Hence 2^k of them; the
    sum telescopes to prod(2^s); and the full subset gives the longest cycle,
    which is what `max_period` names.
    """
    lengths = cycle_lengths(sizes)
    assert len(lengths) == 2 ** len(sizes), (
        f"{sizes}: expected one cycle length per subset of blocks"
    )
    assert sum(lengths) == reduce(lambda a, b: a * b, [2 ** s for s in sizes]), (
        f"{sizes}: subset products should telescope to the state count"
    )
    assert max(lengths) == max_period(sizes), (
        f"{sizes}: the longest cycle should be max_period"
    )


@pytest.mark.parametrize("sizes", BLOCK_SIZES, ids=str)
def test_max_period_is_the_lcm_of_the_component_periods(sizes):
    """For pairwise-coprime sizes, max_period is the lcm of the block periods.

    Each block of size s cycles with period 2^s - 1, so a state with every block
    nonzero repeats only when all of them do at once.  Coprime periods combine
    multiplicatively, and since gcd(2^a - 1, 2^b - 1) = 2^gcd(a,b) - 1 that is
    exactly the pairwise-coprime case -- which is the library's operating domain,
    since MPR block sizes are Mersenne exponents and therefore prime.

    The doubling branch, which applies when a block cannot contribute a fresh
    factor, is covered separately below.
    """
    for i, a in enumerate(sizes):
        for b in sizes[i + 1:]:
            assert math.gcd(a, b) == 1, "this test only covers pairwise-coprime sizes"

    expected = reduce(lambda a, b: a * b // math.gcd(a, b), [2 ** s - 1 for s in sizes])
    assert max_period(sizes) == expected, f"{sizes}: {max_period(sizes)} != lcm {expected}"


@pytest.mark.parametrize("sizes,expected", [
    ([3, 3], 7 * 2),              # repeat: the second copy adds no new factor
    ([5, 5, 5], 31 * 4),          # two repeats, two doublings
    ([3, 5, 5], 7 * 31 * 2),      # one distinct pair plus one repeat
    ([1], 2),                     # size-1 block: 2^1 - 1 = 1 contributes nothing
    ([1, 1], 4),
    ([3, 1], 7 * 2),
], ids=str)
def test_blocks_that_add_no_new_factor_contribute_a_doubling(sizes, expected):
    """A repeated size, or a size-1 block, multiplies the period by 2 instead.

    Such a block cannot extend the cycle by a fresh coprime factor -- a repeat's
    period is already counted, and a size-1 block's own period is 2^1 - 1 = 1.
    Each contributes a factor of 2, which for the size-1 case reflects that it
    behaves as a T-function-style bit driven by its chaining input rather than
    as an MPR cycling on its own.

    Note this is strictly greater than the lcm of the component periods, so
    these are the cases where max_period and that lcm part ways.
    """
    assert max_period(sizes) == expected, (
        f"{sizes}: got {max_period(sizes)}, expected {expected}"
    )


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [4, 5], [3, 4, 5], [2, 3, 5]], ids=str)
def test_observed_cmpr_period_divides_max_period(sizes):
    """Every state's actual period divides lcm(2^s - 1).

    With no chaining the blocks are independent, so a state's period is the lcm
    of its blocks' individual periods -- each of which divides that block's
    2^s - 1.  A period that failed to divide max_period would mean the
    simulation had left the product-of-cycles structure entirely.
    """
    fn = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])
    fn.compile()
    size = len(fn)
    bound = max_period(sizes)

    for seed in [1, (1 << size) - 1, 12345 % (1 << size)]:
        result = FeedbackRegister(seed, fn).period(limit=2 ** 20)
        assert result is not None, f"{sizes} seed={seed}: no period found"
        period, _ = result
        assert bound % period == 0, (
            f"{sizes} seed={seed}: period {period} does not divide max_period {bound}"
        )
