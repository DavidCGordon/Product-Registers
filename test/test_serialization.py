"""Tests for JSON serialization round-trips: FeedbackFunction and FeedbackRegister.

Covers to_JSON/from_JSON and to_file/from_file for MPR, CMPR, Fibonacci, Galois, FCSR,
CrossJoin, and TFunction.  Each round-trip is verified by comparing the output bit
sequences of the original and reconstructed registers.
"""
import os
import random
import tempfile

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, CMPR, Fibonacci, Galois, FCSR, CrossJoin, TFunction
from PyPR.BooleanLogic.ChainingGeneration.Templates import fast_template, arman_template


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
