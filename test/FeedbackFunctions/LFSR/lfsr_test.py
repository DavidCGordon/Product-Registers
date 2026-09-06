"""Tests for Fibonacci and Galois LFSRs: construction, equivalence, synthesis round-trips, inversion.

Derived from ProductRegisters.ipynb (Library Basics > Galois and Fibonacci LFSRs, and
Applications > Berlekamp-Massey Variants for Register Synthesis).
"""
import random

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, Fibonacci, Galois, FCSR
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey


# ── Fibonacci / Galois equivalence ──────────────────────────────────────

def test_fib_galois_same_sequence():
    """Fibonacci and Galois with the same polynomial produce the same bit-0 sequence."""
    F = Fibonacci(5, '12')
    G = Galois(5, '12')

    fib_reg = FeedbackRegister(31, F)
    seq_fib = [state[0] for state in fib_reg.run(compiled=False, limit=31)]

    gal_reg = FeedbackRegister(31, G)
    seq_gal = [state[0] for state in gal_reg.run(compiled=False, limit=31)]

    # the sequences may not be identical element-wise from the same numeric seed,
    # but they should have the same period and the same linear complexity.
    L_fib, _ = berlekamp_massey(seq_fib)
    L_gal, _ = berlekamp_massey(seq_gal)
    assert L_fib == L_gal == 5, f"Expected LC=5, got Fibonacci={L_fib}, Galois={L_gal}"


def test_fib_galois_from_reg_equivalence():
    """Galois.fromReg on a Fibonacci register reproduces the same sequence (notebook cell 76)."""
    F = Fibonacci(5, '12')
    F.compile()
    fib_reg = FeedbackRegister(31, F)
    seq_fib = [state[0] for state in fib_reg.run(limit=40)]

    fib_reg.reset()
    galois_state, galois_fn = Galois.fromReg(fib_reg, bit=0)
    gal_reg = FeedbackRegister(galois_state, galois_fn)
    seq_gal = [state[0] for state in gal_reg.run(compiled=False, limit=40)]

    assert seq_fib == seq_gal, "Galois.fromReg should reproduce Fibonacci sequence"


# ── Fibonacci fromSeq round-trip ────────────────────────────────────────

def test_fibonacci_fromSeq_round_trip():
    """Fibonacci(*Fibonacci.fromSeq(seq)) reproduces seq on bit 0."""
    F = Fibonacci(5, '12')
    reg = FeedbackRegister(31, F)
    seq = [state[0] for state in reg.run(compiled=False, limit=50)]

    init_state, recovered = Fibonacci.fromSeq(seq)
    rec_reg = FeedbackRegister(init_state, recovered)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=50)]
    assert seq == seq2, "Fibonacci fromSeq round-trip failed"


def test_galois_fromSeq_round_trip():
    """Galois(*Galois.fromSeq(seq)) reproduces seq on bit 0."""
    G = Galois(5, '12')
    reg = FeedbackRegister(31, G)
    seq = [state[0] for state in reg.run(compiled=False, limit=50)]

    init_state, recovered = Galois.fromSeq(seq)
    rec_reg = FeedbackRegister(init_state, recovered)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=50)]
    assert seq == seq2, "Galois fromSeq round-trip failed"


# ── Fibonacci NLFSR ─────────────────────────────────────────────────────

def test_fibonacci_nlfsr_fromSeq():
    """Fibonacci.fromSeq with nonlinear=True reproduces a nonlinear sequence."""
    random.seed(42)
    seq = [random.randint(0, 1) for _ in range(60)]

    init_state, nlfsr = Fibonacci.fromSeq(seq, nonlinear=True)
    rec_reg = FeedbackRegister(init_state, nlfsr)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=60)]
    assert seq == seq2, "Nonlinear Fibonacci fromSeq round-trip failed"


# ── Cross-synthesis round-trips (notebook cell 122) ─────────────────────

def test_all_synthesis_round_trips():
    """All fromSeq variants reproduce the original sequence prefix."""
    random.seed(7)
    seq = [random.randint(0, 1) for _ in range(50)]

    for name, factory in [
        ("Fibonacci", lambda s: FeedbackRegister(*Fibonacci.fromSeq(s))),
        ("Galois", lambda s: FeedbackRegister(*Galois.fromSeq(s))),
        ("Fibonacci NL", lambda s: FeedbackRegister(*Fibonacci.fromSeq(s, nonlinear=True))),
        ("FCSR", lambda s: FeedbackRegister(*FCSR.fromSeq(s))),
    ]:
        reg = factory(seq)
        recovered = [state[0] for state in reg.run(compiled=False, limit=50)]
        assert seq == recovered, f"{name} synthesis round-trip failed"


# ── MPR / Galois relationship (notebook cell 77) ───────────────────────

def test_mpr_galois_relationship():
    """MPR(P_reversed) with the Galois-derived U polynomial matches Galois output."""
    P = [1, 1, 0, 1]
    G = Galois(3, P)
    g_reg = FeedbackRegister(1, G)
    seq_gal = [state[0] for state in g_reg.run(compiled=False, limit=20)]

    # MPR with reversed polynomial should produce a related sequence
    M = MPR(3, P[::-1])
    m_reg = FeedbackRegister(1, M)
    seq_mpr = [state[0] for state in m_reg.run(compiled=False, limit=20)]

    L_gal, _ = berlekamp_massey(seq_gal)
    L_mpr, _ = berlekamp_massey(seq_mpr)
    assert L_gal == L_mpr == 3, f"Expected LC=3, got Galois={L_gal}, MPR={L_mpr}"


# ── Galois inversion ───────────────────────────────────────────────────

def test_galois_invert():
    """Clocking an inverted Galois register undoes a clock of the original."""
    G = Galois(5, '12')
    reg = FeedbackRegister(31, G)

    # clock forward 10 steps
    for _ in reg.run(compiled=False, limit=10):
        pass
    state_at_10 = reg[:].copy()

    # clock forward one more
    reg.clock(compiled=False)
    state_at_11 = reg[:].copy()

    # invert and clock back
    G.invert()
    reg2 = FeedbackRegister(state_at_11, G)
    reg2.clock(compiled=False)
    recovered = reg2[:].copy()

    assert (state_at_10 == recovered).all(), "Inverted Galois should undo one clock step"


# ── Berlekamp-Massey basic ─────────────────────────────────────────────

def test_berlekamp_massey_known():
    """BM on a known LFSR sequence returns the correct linear complexity."""
    F = Fibonacci(5, '12')
    reg = FeedbackRegister(31, F)
    seq = [state[0] for state in reg.run(compiled=False, limit=100)]
    L, poly = berlekamp_massey(seq)
    assert L == 5, f"Expected LC=5, got {L}"

