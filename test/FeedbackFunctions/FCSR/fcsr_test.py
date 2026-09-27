"""Tests for FCSR: construction, fromSeq round-trip, state_from_frac, 2-adic structure.

Derived from ProductRegisters.ipynb (Library Basics > FCSRs, and
Applications > Berlekamp-Massey Variants for Register Synthesis).

An FCSR is the 2-adic analogue of an LFSR: it generates the binary expansion
of a rational p/q in the 2-adic integers, with q the connection integer.  The
final section tests the consequences -- clocking multiplies the expansion by
2, so cycle lengths are governed by the multiplicative order of 2 modulo q,
and the carry bits make the update non-bijective so most states sit on a
transient.
"""
import random

import pytest

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import FCSR
from PyPR.Tools.RegisterSynthesis.fcsrSynthesis import BM_FCSR, FCSR_size, fcsr_eval


# ── Construction ────────────────────────────────────────────────────────

def test_fcsr_construction():
    """FCSR(5, 27) should create a register with size = 2*5 - 1 = 9."""
    F = FCSR(5, 27)
    assert F.size == 9, f"Expected size 9, got {F.size}"
    assert F.connection_int == 27

def test_fcsr_produces_binary_output():
    F = FCSR(5, 27)
    reg = FeedbackRegister(2**F.size - 1, F)
    output = [state[0] for state in reg.run(compiled=False, limit=50)]
    assert all(b in (0, 1) for b in output)


# ── state_from_frac ────────────────────────────────────────────────────

def test_state_from_frac_edge_cases():
    """Edge cases: 0/1 and 1/1.

    0/1 is the one fraction needing no register at all -- a single cell holding
    0, emitting zeros forever.

    1/1 is not symmetric with it. A positive numerator has a terminating
    expansion (1,0,0,...), which the register can only produce with its feedback
    left open, i.e. den < 2**size - 1. At size 1 that fails with equality, and
    the register instead emits the expansion of -1/1, all ones. So 1/1 needs two
    value cells where 0/1 needs one; the expansion itself is checked in
    test_state_from_frac_emits_the_expansion_of_that_fraction below.
    """
    size, state = FCSR.state_from_frac(0, 1)
    assert size == 1
    assert state == [0]

    size, state = FCSR.state_from_frac(1, 1)
    assert size == 2, "1/1 needs the feedback left open, which size 1 does not give"
    assert len(state) == 2 * size - 1

def test_state_from_frac_positive_num():
    size, state = FCSR.state_from_frac(3, 7)
    assert size >= 1
    assert len(state) == 2 * size - 1

def test_state_from_frac_negative_num():
    size, state = FCSR.state_from_frac(-5, 7)
    assert size >= 1
    assert len(state) == 2 * size - 1


# ── fromSeq round-trip ──────────────────────────────────────────────────

def test_fcsr_fromSeq_round_trip():
    """FCSR.fromSeq should recover a register that reproduces the input sequence."""
    F = FCSR(3, 5)
    reg = FeedbackRegister(2**F.size - 1, F)
    seq = [state[0] for state in reg.run(compiled=False, limit=80)]

    init_state, recovered = FCSR.fromSeq(seq)
    rec_reg = FeedbackRegister(init_state, recovered)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=80)]
    assert seq == seq2, "FCSR fromSeq round-trip failed"

def test_fcsr_fromSeq_random():
    """fromSeq on a random sequence should reproduce the sequence prefix."""
    random.seed(42)
    seq = [random.randint(0, 1) for _ in range(60)]

    init_state, recovered = FCSR.fromSeq(seq)
    rec_reg = FeedbackRegister(init_state, recovered)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=60)]
    assert seq == seq2, "FCSR fromSeq on random sequence failed"


# ── BM_FCSR ─────────────────────────────────────────────────────────────

def test_bm_fcsr_returns_valid():
    """BM_FCSR should return (size, numerator, denominator)."""
    random.seed(7)
    seq = [random.randint(0, 1) for _ in range(60)]
    size, num, den = BM_FCSR(seq)
    assert isinstance(size, int) and size >= 1
    assert isinstance(num, int)
    assert isinstance(den, int) and den >= 1

def test_bm_fcsr_all_zeros():
    """All-zero sequence should return minimal complexity."""
    seq = [0] * 30
    size, num, den = BM_FCSR(seq)
    assert size == 1
    assert num == 0
    assert den == 1


# ── carries / values properties ─────────────────────────────────────────

def test_carries_values():
    F = FCSR(4, 11)
    assert len(F.carries) == F.size // 2
    assert len(F.values) == F.size // 2


# ── 2-adic period structure ───────────────────────────────────────────────────

CONNECTIONS = [(4, 13), (5, 19), (5, 29), (5, 37), (6, 43)]


@pytest.mark.parametrize("complexity,q", CONNECTIONS, ids=lambda v: str(v))
def test_longest_cycle_has_length_equal_to_the_order_of_two_mod_q(complexity, q):
    """The maximum period over all states is exactly ord_q(2), and it is attained.

    This is the property that makes q the FCSR's analogue of a primitive
    polynomial: it fixes the length of the main cycle.  Asserting the maximum
    rather than every state's period is the precise claim -- degenerate states
    exist (see the spectrum test) and sit on shorter cycles.
    """
    order, value = 1, 2 % q
    while value != 1:
        value = (value * 2) % q
        order += 1

    fn = FCSR(complexity, q)
    fn.compile()
    width = 2 * complexity - 1

    periods = set()
    for seed in range(2 ** width):
        result = FeedbackRegister(seed, fn).period(limit=2 ** 16)
        assert result is not None, f"q={q} seed={seed}: Brent's algorithm found no cycle"
        periods.add(result[0])

    assert max(periods) == order, (
        f"q={q}: longest cycle is {max(periods)}, expected ord_{q}(2) = {order}"
    )


@pytest.mark.parametrize("complexity,q", CONNECTIONS, ids=lambda v: str(v))
def test_period_spectrum_is_exactly_one_and_the_order_of_two(complexity, q):
    """For prime q the only cycle lengths are 1 and ord_q(2).

    A reachable 2-adic value p/q either has gcd(p, q) = 1, giving period
    ord_q(2), or reduces to an integer, giving a constant expansion of period
    1.  With q prime there is nothing in between.  An intermediate period would
    mean the register had reached a value whose denominator is a proper
    divisor of q -- impossible unless the connection integer is not the
    denominator the implementation actually realizes.
    """
    order, value = 1, 2 % q
    while value != 1:
        value = (value * 2) % q
        order += 1

    fn = FCSR(complexity, q)
    fn.compile()
    width = 2 * complexity - 1

    periods = set()
    for seed in range(2 ** width):
        result = FeedbackRegister(seed, fn).period(limit=2 ** 16)
        assert result is not None, f"q={q} seed={seed}: Brent's algorithm found no cycle"
        periods.add(result[0])
    assert periods == {1, order}, (
        f"q={q}: observed cycle lengths {sorted(periods)}, expected {{1, {order}}}"
    )


@pytest.mark.parametrize("complexity,q", CONNECTIONS, ids=lambda v: str(v))
def test_every_period_divides_the_order_of_two_mod_q(complexity, q):
    """No state's period exceeds or fails to divide ord_q(2).

    Shifting is multiplication by 2 modulo the realized denominator, which
    divides q, so every cycle length divides ord_q(2).  This is the weaker
    statement that survives for composite q as well, where the spectrum test
    above would not.
    """
    order, value = 1, 2 % q
    while value != 1:
        value = (value * 2) % q
        order += 1

    fn = FCSR(complexity, q)
    fn.compile()
    width = 2 * complexity - 1

    for seed in range(2 ** width):
        result = FeedbackRegister(seed, fn).period(limit=2 ** 16)
        assert result is not None, f"q={q} seed={seed}: Brent's algorithm found no cycle"
        period, _ = result
        assert order % period == 0, (
            f"q={q} seed={seed}: period {period} does not divide ord_{q}(2) = {order}"
        )


@pytest.mark.parametrize("complexity,q", CONNECTIONS, ids=lambda v: str(v))
def test_register_width_is_twice_the_complexity_minus_one(complexity, q):
    """An FCSR of 2-adic complexity d occupies 2d - 1 bits: d values, d - 1 carries.

    The carry cells are what distinguish an FCSR from an LFSR of the same
    length, and they are why the state count exceeds the number of distinct
    output phases -- which is what makes the update map non-injective.
    """
    assert len(FCSR(complexity, q)) == 2 * complexity - 1


@pytest.mark.parametrize("complexity,q", CONNECTIONS, ids=lambda v: str(v))
def test_some_states_lie_on_a_transient_into_the_cycle(complexity, q):
    """At least one seed has a nonzero preperiod.

    Carry propagation discards information, so the update map is not a
    bijection and the state graph is "rho"-shaped rather than a disjoint union
    of cycles.  A register reporting preperiod 0 for every state would be
    behaving like an LFSR, meaning the carries were doing nothing.
    """
    fn = FCSR(complexity, q)
    fn.compile()
    width = 2 * complexity - 1

    preperiods = []
    for seed in range(2 ** width):
        result = FeedbackRegister(seed, fn).period(limit=2 ** 16)
        assert result is not None, f"q={q} seed={seed}: Brent's algorithm found no cycle"
        preperiods.append(result[1])
    assert any(p > 0 for p in preperiods), (
        f"q={q}: no state had a transient, so the update map looks bijective"
    )


# -- state_from_frac agrees with the 2-adic expansion -------------------------
#
# state_from_frac takes the true fraction: p may be either sign, q is the
# positive connection integer. The register built from it must emit the 2-adic
# expansion of p/q exactly. Sign selects the regime -- a non-positive numerator
# runs with the feedback closed, a positive one has a terminating expansion and
# needs it left open (q < 2**d - 1), which is what the sizing arranges.
#
# The sweep matters here: p = q = 1 is the single input where the log-based size
# is tight against that bound, and two fixed round-trip cases will not find it.

@pytest.mark.parametrize("num,den", [
    (1, 1),     # the tight case: needs size 2, not 1
    (-1, 1),    # all ones -- distinct from 1/1, and the pair is the regression
    (0, 1),
    (2, 1), (3, 1), (5, 1), (9, 1),
    (-2, 1), (-3, 1),
    (-1, 7), (-3, 7), (-7, 7),
    (-8, 5), (5, 3), (-13, 11),
])
def test_state_from_frac_emits_the_expansion_of_that_fraction(num, den):
    """The register built from p/q emits the 2-adic expansion of p/q."""
    # expansion of p/q, low bit first: s = p mod 2, then p <- (p - s*q)/2
    expansion, p = [], num
    for _ in range(14):
        bit = p % 2
        expansion.append(bit)
        p = (p - bit * den) // 2

    size, state = FCSR.state_from_frac(num, den)
    register = FCSR(size, den)
    emitted = [int(s[0]) for s
               in FeedbackRegister(state, register).run(compiled=False, limit=14)]

    assert emitted == expansion, (
        f"{num}/{den}: register emitted {emitted}, expansion is {expansion}"
    )


def test_positive_numerators_leave_the_feedback_open():
    """A positive numerator needs q < 2**d - 1, or its expansion cannot terminate.

    With the feedback closed every reachable numerator is non-positive, so the
    emitted sequence is never eventually zero. The sizing in state_from_frac is
    what keeps the positive case on the open side of that bound.
    """
    for num in range(1, 30):
        size, _ = FCSR.state_from_frac(num, 1)
        assert 1 < 2 ** size - 1, (
            f"num={num}: size {size} closes the feedback, so +{num} cannot be emitted"
        )


# -- FCSR_size agrees with state_from_frac -----------------------------------
#
# Two functions compute the register size for a fraction: FCSR_size (which
# BM_FCSR returns) and state_from_frac (which fromSeq actually uses). They must
# agree, or BM_FCSR hands a size to callers that does not fit the state
# state_from_frac would build. p = q = 1 is the case where the logarithmic
# sizing is tight and both special-case it -- a positive numerator needs the
# feedback left open, and size 1 with q = 1 closes the loop.

def test_fcsr_size_agrees_with_state_from_frac():
    """The two sizing routines return the same d for every fraction."""
    disagreements = []
    for num in range(-40, 41):
        for den in range(1, 22, 2):
            if FCSR_size(num, den) != FCSR.state_from_frac(num, den)[0]:
                disagreements.append((num, den))
    assert not disagreements, (
        f"FCSR_size and state_from_frac disagree on {disagreements[:5]}"
    )


def test_fcsr_size_leaves_the_feedback_open_for_positive_numerators():
    """A positive numerator needs q < 2**d - 1, which is what forces d = 2 at 1/1.

    A terminating expansion is only reachable with the feedback open. At
    p = q = 1 the logarithmic sizing gives d = 1, where q = 2**d - 1 closes the
    loop and the register emits -1/1 instead.
    """
    for num in range(1, 30):
        assert 1 < 2 ** FCSR_size(num, 1) - 1, (
            f"num={num}: FCSR_size closes the feedback, so +{num} cannot be emitted"
        )

