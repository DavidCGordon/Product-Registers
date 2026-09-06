"""Tests for FeedbackRegister simulation: clock, run, reset, seed, period.

Derived from ProductRegisters.ipynb (Library Basics > Simulating a Function in a Register).
"""
from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, CMPR, Fibonacci, Galois
from PyPR.BooleanLogic import AND, VAR
from PyPR.BooleanLogic.ChainingGeneration.Templates import fast_template


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
    period, _ = reg.period(compiled=False)
    assert period == 7, f"Expected period 7, got {period}"

def test_period_galois():
    """Galois LFSR with same polynomial should also have period 7."""
    G = Galois(3, [1, 1, 0, 1])
    reg = FeedbackRegister(1, G)
    period, _ = reg.period(compiled=False)
    assert period == 7, f"Expected period 7, got {period}"

def test_period_mpr():
    """MPR with primitive poly of degree 3 has period 2^3 - 1 = 7."""
    M = MPR(3, [1, 0, 1, 1])
    reg = FeedbackRegister(1, M)
    period, _ = reg.period(compiled=False)
    assert period == 7, f"Expected period 7, got {period}"
