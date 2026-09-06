"""Tests for TFunction: construction, binary counter behavior, CMPR interface, period.

Derived from ProductRegisters.ipynb (Library Basics > Other Function Families,
and ANF Iteration > An interesting observation with a TFunction).
"""
from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import TFunction, CMPR, MPR
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template


# ── Construction ────────────────────────────────────────────────────────

def test_tfunction_construction():
    T = TFunction(5)
    assert T.size == 5

def test_tfunction_is_cmpr():
    """TFunction subclasses CMPR."""
    T = TFunction(5)
    assert isinstance(T, CMPR)

def test_induction_order():
    T = TFunction(5)
    assert T.induction_order == [4, 3, 2, 1, 0]


# ── Binary counter behavior ────────────────────────────────────────────

def test_binary_counter_period():
    """Default TFunction(n) is a binary counter with period 2^n."""
    T = TFunction(4)
    reg = FeedbackRegister(0, T)
    period, _ = reg.period(compiled=False)
    assert period == 2**4, f"Expected period 16, got {period}"

def test_binary_counter_full_orbit():
    """The counter visits all 2^n states."""
    T = TFunction(4)
    reg = FeedbackRegister(0, T)
    states = set()
    for state in reg.run(compiled=False, limit=16):
        states.add(tuple(state[:].copy()))
    assert len(states) == 16, f"Expected 16 distinct states, got {len(states)}"

def test_binary_counter_increments():
    """Verify the counter increments by 1 each clock (with bit n-1 as LSB)."""
    T = TFunction(4)
    reg = FeedbackRegister(0, T)
    for expected in range(16):
        state = reg[:].copy()
        # bit n-1 is LSB, bit 0 is MSB
        value = sum(state[T.size - 1 - i] * (2**i) for i in range(T.size))
        assert value == expected, f"At step {expected}, got value {value}"
        reg.clock(compiled=False)


# ── Chaining on TFunction ──────────────────────────────────────────────

def test_tfunction_with_chaining():
    """TFunction can accept chaining like any CMPR (notebook cell 111)."""
    M7 = MPR(7, [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0])
    M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
    M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
    T = TFunction(5)
    C = CMPR([T, M7, M5, M3])
    C.generateChaining(template=arman_template())
    reg = FeedbackRegister(2**C.size - 1, C)
    output = [state[0] for state in reg.run(compiled=False, limit=50)]
    assert all(b in (0, 1) for b in output)
    assert len(output) == 50


# ── Root expressions / monomial profiles ────────────────────────────────

def test_tfunction_monomial_profiles():
    """TFunction should be able to compute monomial profiles as a CMPR."""
    T = TFunction(5)
    M7 = MPR(7, [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0])
    M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
    M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
    C = CMPR([T, M7, M5, M3])
    C.generateChaining(template=arman_template())
    MPs = C.monomial_profiles()
    assert len(MPs) == C.size

