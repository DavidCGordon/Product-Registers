"""Tests for JSON serialization round-trips: FeedbackFunction and FeedbackRegister.

Covers to_JSON/from_JSON and to_file/from_file for MPR, CMPR, Fibonacci, Galois, FCSR,
CrossJoin, and TFunction.  Each round-trip is verified by comparing the output bit
sequences of the original and reconstructed registers.
"""
import os
import random
import tempfile

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

# ── Helpers ──────────────────────────────────────────────────────────────────

def _seq(reg, limit=20):
    return [int(state[0]) for state in reg.run(compiled=False, limit=limit)]

def _ff_json_round_trip(fn, seed, limit=20):
    """Serialize a FeedbackFunction to JSON and back; compare output sequences."""
    reg1 = FeedbackRegister(seed, fn)
    seq1 = _seq(reg1, limit)

    json_data = fn.to_JSON()
    fn2 = type(fn).from_JSON(json_data)

    reg2 = FeedbackRegister(seed, fn2)
    seq2 = _seq(reg2, limit)
    assert seq1 == seq2, f"{type(fn).__name__} JSON round-trip changed the output sequence"


# ── FeedbackFunction JSON round-trips ─────────────────────────────────────────

def test_mpr_json_round_trip():
    _ff_json_round_trip(MPR(5, "12"), seed=1)

def test_fibonacci_json_round_trip():
    _ff_json_round_trip(Fibonacci(5, "12"), seed=31)

def test_galois_json_round_trip():
    _ff_json_round_trip(Galois(5, "12"), seed=31)

def test_fcsr_json_round_trip():
    F = FCSR(3, 5)
    _ff_json_round_trip(F, seed=2**F.size - 1)

def test_cmpr_json_round_trip():
    random.seed(0)
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    C.generateChaining(template=fast_template())
    _ff_json_round_trip(C, seed=2**C.size - 1)

def test_crossjoin_json_round_trip():
    random.seed(0)
    CJ = CrossJoin(7, "65")
    CJ.generateNonlinearity(maxAnds=3)
    _ff_json_round_trip(CJ, seed=2**7 - 1)

def test_tfunction_json_round_trip():
    T = TFunction(4)
    _ff_json_round_trip(T, seed=0, limit=16)


# ── FeedbackRegister JSON round-trip ──────────────────────────────────────────

def test_register_json_round_trip():
    """FeedbackRegister.to_JSON/from_JSON preserves the output sequence."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    seq1 = _seq(reg)

    reg.reset()
    json_data = reg.to_JSON()
    reg2 = FeedbackRegister.from_JSON(json_data)
    seq2 = _seq(reg2)

    assert seq1 == seq2

def test_register_json_round_trip_cmpr():
    """FeedbackRegister wrapping a CMPR serializes and restores correctly."""
    random.seed(0)
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    C = CMPR([M7, M5])
    C.generateChaining(template=fast_template())
    seed = 2**C.size - 1

    reg = FeedbackRegister(seed, C)
    seq1 = _seq(reg, limit=30)

    reg.reset()
    json_data = reg.to_JSON()
    reg2 = FeedbackRegister.from_JSON(json_data)
    seq2 = _seq(reg2, limit=30)

    assert seq1 == seq2


# ── FeedbackRegister file round-trip ──────────────────────────────────────────

def test_register_file_round_trip():
    """to_file/from_file preserves the output sequence."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    seq1 = _seq(reg)
    reg.reset()

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = f.name
    try:
        reg.to_file(path)
        reg2 = FeedbackRegister.from_file(path)
        seq2 = _seq(reg2)
        assert seq1 == seq2
    finally:
        os.unlink(path)

def test_feedbackfunction_file_round_trip():
    """FeedbackFunction.to_file/from_file preserves the output sequence."""
    M = MPR(7, "65")
    reg1 = FeedbackRegister(1, M)
    seq1 = _seq(reg1, limit=30)

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = f.name
    try:
        M.to_file(path)
        M2 = MPR.from_file(path)
        reg2 = FeedbackRegister(1, M2)
        seq2 = _seq(reg2, limit=30)
        assert seq1 == seq2
    finally:
        os.unlink(path)


# ── Round-trip equality on the whole state space ─────────────────────────────

def test_json_round_trip_agrees_on_every_state():
    """A deserialized feedback function agrees with the original on all 2^n states.

    The sequence comparisons above follow one orbit from one seed, so they only
    exercise the bits that orbit happens to visit.  Comparing every bit's update
    function at every state is the full contract: it catches a bit whose DAG was
    dropped or rewired during serialization even when the bit-0 sequence is
    unaffected, which a single-orbit check cannot distinguish from success.

    Sizes are kept small because the comparison is exhaustive in the state space.
    """
    cases = [
        MPR(5, "12"),
        Fibonacci(5, "12"),
        Galois(5, "12"),
        FCSR(4, 13),
        TFunction(5),
        CMPR([MPR(3, [1, 1, 0, 1]), MPR(4, [1, 1, 0, 0, 1])]),
    ]

    for fn in cases:
        restored = type(fn).from_JSON(fn.to_JSON())
        size = len(fn)
        assert len(restored) == size, f"{type(fn).__name__}: round-trip changed the size"

        for value in range(2 ** size):
            state = [(value >> k) & 1 for k in range(size)]
            original_bits = [bit_fn.eval(state) for bit_fn in fn.fn_list]
            restored_bits = [bit_fn.eval(state) for bit_fn in restored.fn_list]
            assert original_bits == restored_bits, (
                f"{type(fn).__name__}: round-trip differs at state {state} "
                f"({original_bits} vs {restored_bits})"
            )
