"""Tests for FCSR: construction, fromSeq round-trip, state_from_frac.

Derived from ProductRegisters.ipynb (Library Basics > FCSRs, and
Applications > Berlekamp-Massey Variants for Register Synthesis).
"""
import random

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import FCSR
from PyPR.Tools.RegisterSynthesis.fcsrSynthesis import BM_FCSR


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
    """Edge cases: 0/1 and 1/1."""
    size, state = FCSR.state_from_frac(0, 1)
    assert size == 1
    assert state == [0]

    size, state = FCSR.state_from_frac(1, 1)
    assert size == 1
    assert state == [1]

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

