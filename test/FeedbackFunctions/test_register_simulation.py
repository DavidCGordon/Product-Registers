"""Tests for FeedbackRegister simulation: clock, run, reset, seed, period.

Derived from ProductRegisters.ipynb (Library Basics > Simulating a Function in a Register).

The final section checks the redundant execution paths against each other.
`FeedbackRegister` carries four period algorithms (compiled x safe) and two
clocking paths (numba-compiled vs. interpreted ANF evaluation).  These are
performance variants of one mathematical object, so any disagreement between
them is a bug in one of them -- there is no design freedom there.
"""
import numpy as np
import pytest

from PyPR.BooleanLogic import AND, VAR
from PyPR.BooleanLogic.ChainingGeneration.Templates import fast_template

from PyPR.FeedbackFunctions import (
    CMPR,
    FCSR,
    MPR,
    CrossJoin,
    Fibonacci,
    Galois,
    TFunction,
)
from PyPR.FeedbackRegister import FeedbackRegister

# ── Fixtures ────────────────────────────────────────────────────────────

def make_small_cmpr():
    M7 = MPR(7, [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0])
    M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
    M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
    M2 = MPR(2, [1, 1, 1], [1, 1])
    C = CMPR([M7, M5, M3, M2])
    C.generateChaining(template=fast_template())
    return C


# ── Basic simulation ───────────────────────────────────────────────────

def test_clock_changes_state():
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    state_before = F[:].copy()
    F.clock(compiled=False)
    state_after = F[:].copy()
    assert not (state_before == state_after).all(), "State should change after clock"

def test_run_produces_correct_count():
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    states = [s[:].copy() for s in F.run(compiled=False, limit=10)]
    assert len(states) == 10, f"Expected 10 states, got {len(states)}"

def test_run_states_are_distinct():
    """Consecutive states should be distinct (non-degenerate seed)."""
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    states = [tuple(s[:].copy()) for s in F.run(compiled=False, limit=20)]
    assert len(set(states)) == 20, "Expected 20 distinct states"

def test_run_returns_reference():
    """run() yields a reference — if you don't copy, all entries are the same."""
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    state_list = [state for state in F.run(compiled=False, limit=5)]
    # all entries are the same object (run yields self)
    for s in state_list:
        assert s is state_list[-1], "All refs should point to same register object"


# ── Reset and seed ─────────────────────────────────────────────────────

def test_reset_returns_to_seed():
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    initial = F[:].copy()
    for _ in F.run(compiled=False, limit=20):
        pass
    F.reset()
    assert (F[:] == initial).all(), "reset() should restore the seed state"

def test_seed_changes_initial_state():
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    initial_1 = F[:].copy()
    F.seed(1)
    F.reset()
    initial_2 = F[:].copy()
    assert not (initial_1 == initial_2).all(), "Different seeds should give different states"

@pytest.mark.parametrize("bits", [[1], [1, 0, 1, 0]])
def test_seed_and_state_reject_wrong_length(bits):
    """A bit list must cover the register exactly -- a short one used to be
    accepted and left a state array shorter than the register."""
    F = Fibonacci(3, [1, 1, 0, 1])
    reg = FeedbackRegister(1, F)
    with pytest.raises(ValueError, match="exactly 3 bits"):
        reg.seed(bits)
    with pytest.raises(ValueError, match="exactly 3 bits"):
        reg.set_state(bits)

def test_seed_array_is_copied():
    """The register owns its seed: mutating the caller's array afterwards
    must not change what reset() restores."""
    F = Fibonacci(3, [1, 1, 0, 1])
    bits = np.array([1, 0, 0], dtype=np.uint8)
    reg = FeedbackRegister(bits, F)
    bits[0] = 0
    reg.reset()
    assert reg[:].tolist() == [1, 0, 0]

def test_negative_int_seed_rejected():
    F = Fibonacci(3, [1, 1, 0, 1])
    with pytest.raises(ValueError, match="outside register capacity"):
        FeedbackRegister(-1, F)

def test_compiled_clock_on_uncompiled_fn_raises_value_error():
    """FeedbackFunction.__init__ sets _compiled = None, so the guard must test
    for None rather than for the attribute's presence."""
    F = TFunction(4)
    reg = FeedbackRegister(1, F)
    with pytest.raises(ValueError, match="not compiled"):
        reg.clock()
    with pytest.raises(ValueError, match="not compiled"):
        next(reg.run(1))
    with pytest.raises(ValueError, match="not compiled"):
        reg.period()


# ── Output stream and filtering ────────────────────────────────────────

def test_output_stream_is_binary():
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    output = [state[0] for state in F.run(compiled=False, limit=50)]
    assert all(b in (0, 1) for b in output)

def test_filtered_output():
    C = make_small_cmpr()
    F = FeedbackRegister(2**C.size - 1, C)
    filt = AND(VAR(0), VAR(1))
    filtered = [filt.eval(state) for state in F.run(compiled=False, limit=50)]
    assert all(b in (0, 1) for b in filtered)


# ── Period ──────────────────────────────────────────────────────────────

def test_period_fibonacci():
    """Fibonacci LFSR with primitive poly of degree 3 has period 2^3 - 1 = 7."""
    F = Fibonacci(3, [1, 1, 0, 1])
    reg = FeedbackRegister(1, F)
    # period() reports None when it hits its limit without closing a cycle;
    # these registers always close, so pin that before unpacking
    result = reg.period(compiled=False)
    assert result is not None
    period, _ = result
    assert period == 7, f"Expected period 7, got {period}"

def test_period_galois():
    """Galois LFSR with same polynomial should also have period 7."""
    G = Galois(3, [1, 1, 0, 1])
    reg = FeedbackRegister(1, G)
    # period() reports None when it hits its limit without closing a cycle;
    # these registers always close, so pin that before unpacking
    result = reg.period(compiled=False)
    assert result is not None
    period, _ = result
    assert period == 7, f"Expected period 7, got {period}"

def test_period_mpr():
    """MPR with primitive poly of degree 3 has period 2^3 - 1 = 7."""
    M = MPR(3, [1, 0, 1, 1])
    reg = FeedbackRegister(1, M)
    # period() reports None when it hits its limit without closing a cycle;
    # these registers always close, so pin that before unpacking
    result = reg.period(compiled=False)
    assert result is not None
    period, _ = result
    assert period == 7, f"Expected period 7, got {period}"


# ── Agreement between the redundant execution paths ───────────────────────────

PRIMITIVE = {
    3: [1, 1, 0, 1],
    4: [1, 1, 0, 0, 1],
    5: [1, 0, 1, 0, 0, 1],
    6: [1, 1, 0, 0, 0, 0, 1],
}


def every_feedback_family():
    """One instance of each concrete FeedbackFunction, paired with a usable seed.

    A function rather than a module constant because CrossJoin has to be mutated
    after construction, and a function rather than two inline lists because both
    tests below parametrize over the same set -- duplicating a seven-entry
    construction would be the less readable option.
    """
    crossjoin = CrossJoin(6, PRIMITIVE[6])
    crossjoin.generateNonlinearity(2)
    return [
        ("MPR", MPR(5, PRIMITIVE[5]), 31),
        ("Fibonacci", Fibonacci(5, "12"), 31),
        ("Galois", Galois(5, "12"), 31),
        ("FCSR", FCSR(5, 37), 100),
        ("TFunction", TFunction(6), 21),
        ("CMPR", CMPR([MPR(3, PRIMITIVE[3]), MPR(4, PRIMITIVE[4])]), 100),
        ("CrossJoin", crossjoin, 21),
    ]


@pytest.mark.parametrize(("name", "fn", "seed"), every_feedback_family(), ids=lambda v: v if isinstance(v, str) else "")
def test_compiled_and_uncompiled_runs_agree(name, fn, seed):
    """The numba path and the ANF-evaluation path produce identical state sequences.

    `_clock_compiled` swaps state pointers and writes through a numba kernel,
    while `_clock_uncompiled` evaluates each bit's BooleanFunction DAG.  A
    divergence means the compiled kernel and the DAG disagree about the
    feedback function -- e.g. a stale compile, or an aliasing error in the
    pointer swap.
    """
    fn.compile()
    interpreted = [int(s) for s in FeedbackRegister(seed, fn).run(compiled=False, limit=60)]
    compiled = [int(s) for s in FeedbackRegister(seed, fn).run(compiled=True, limit=60)]
    assert interpreted == compiled, f"{name}: compiled and uncompiled runs diverged"


@pytest.mark.parametrize(("name", "fn", "_seed"), every_feedback_family(), ids=lambda v: v if isinstance(v, str) else "")
@pytest.mark.parametrize("seed_kind", ["one", "all_ones", "arbitrary"])
def test_period_variants_agree_where_the_unsafe_search_is_valid(name, fn, _seed, seed_kind):
    """The four period algorithms obey the contract their docstring states.

    `safe=True` runs Brent's cycle-finding, which handles a state that sits on
    a transient leading into the cycle.  `safe=False` is a naive search for a
    repeat of the *initial* state, so it cannot terminate when that state is
    not itself on the cycle -- exactly the non-bijective case the docstring
    warns about.  The contract is therefore three claims, not one equality:

      1. compiling changes speed, never the answer, on both algorithms;
      2. the unsafe search returns None precisely when the preperiod is
         nonzero (the limit below is orders of magnitude above every period
         here, which is what lets a None be read as "the start state is not on
         the cycle" rather than "ran out of iterations");
      3. whenever it does return, it agrees with Brent's answer.

    FCSR is the register that exercises the None branch: its carry bits make
    the update non-bijective, so most seeds land on a transient.
    """
    fn.compile()
    size = len(fn)
    seed = {"one": 1, "all_ones": (1 << size) - 1, "arbitrary": 100 % (1 << size)}[seed_kind]

    reg = FeedbackRegister(seed, fn)
    results = {}
    for compiled in (True, False):
        for safe in (True, False):
            reg.reset()
            results[(compiled, safe)] = reg.period(compiled=compiled, safe=safe, limit=2 ** 12)

    safe_compiled, safe_plain = results[(True, True)], results[(False, True)]
    unsafe_compiled, unsafe_plain = results[(True, False)], results[(False, False)]

    assert safe_compiled == safe_plain, (
        f"{name}/{seed_kind}: compiling changed the safe period: {results}"
    )
    assert unsafe_compiled == unsafe_plain, (
        f"{name}/{seed_kind}: compiling changed the unsafe period: {results}"
    )

    assert safe_compiled is not None, f"{name}/{seed_kind}: Brent's algorithm found no cycle"
    period, preperiod = safe_compiled

    assert (unsafe_compiled is None) == (preperiod != 0), (
        f"{name}/{seed_kind}: preperiod is {preperiod} but the unsafe search "
        f"returned {unsafe_compiled}"
    )
    if unsafe_compiled is not None:
        assert unsafe_compiled == (period, preperiod), (
            f"{name}/{seed_kind}: unsafe {unsafe_compiled} != safe {safe_compiled}"
        )
