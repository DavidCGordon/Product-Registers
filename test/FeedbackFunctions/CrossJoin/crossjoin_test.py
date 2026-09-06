"""Tests for CrossJoin: construction, nonlinearity generation, compensation, filter decomposition.

Derived from ProductRegisters.ipynb (Library Basics > Other Function Families,
and Cryptanalysis > Espresso Cryptanalysis).
"""
import random

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import Fibonacci, CrossJoin
from PyPR.BooleanLogic import VAR
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey


# ── Construction ────────────────────────────────────────────────────────

def test_crossjoin_construction():
    C = CrossJoin(7, '69')
    assert C.size == 7
    assert C.tau == 6  # before nonlinearity, tau = size - 1

def test_crossjoin_koopman_vs_list():
    """Koopman string and explicit list should produce equivalent registers."""
    C1 = CrossJoin(5, '12')
    C2 = CrossJoin(5, [1, 0, 0, 1, 0, 1])
    # compare fn_list outputs on a test state
    state = [1, 0, 1, 1, 0]
    for i in range(5):
        assert C1[i].eval(state) == C2[i].eval(state), f"Mismatch at bit {i}"


# ── Nonlinearity generation ────────────────────────────────────────────

def test_generate_nonlinearity():
    random.seed(42)
    C = CrossJoin(7, '69')
    C.generateNonlinearity(maxAnds=4, tapDensity=0.75)
    # tau should have been set
    assert C.tau <= int(0.75 * 7)
    # the register should still produce binary output
    reg = FeedbackRegister(2**7 - 1, C)
    output = [state[0] for state in reg.run(compiled=False, limit=50)]
    assert all(b in (0, 1) for b in output)


# ── Compensation list ───────────────────────────────────────────────────

def test_compensation_list_length():
    random.seed(42)
    C = CrossJoin(7, '69')
    C.generateNonlinearity()
    comp = C.compensation_list()
    assert len(comp) == C.size

def test_compensation_preserves_sequence():
    """Applying compensation to the base LFSR should match the crossjoin output."""
    random.seed(42)
    C = CrossJoin(7, '69')
    C.generateNonlinearity(maxAnds=3, tapDensity=0.75)

    seed = 2**7 - 1
    cj_reg = FeedbackRegister(seed, C)
    seq_cj = [state[0] for state in cj_reg.run(compiled=False, limit=31)]

    # decompose into filter generator
    base_lfsr, filters = C.filter_generator()
    base_reg = FeedbackRegister(C.convert_state(cj_reg._seed), base_lfsr)
    seq_filtered = [filters[0].eval(state) for state in base_reg.run(compiled=False, limit=31)]

    assert seq_cj == seq_filtered, "Compensation-based filter should match crossjoin output"


# ── Filter generator ───────────────────────────────────────────────────

def test_filter_generator_returns_pair():
    C = CrossJoin(5, '12')
    base, filters = C.filter_generator()
    assert isinstance(base, Fibonacci)
    assert len(filters) == 5


# ── convert_state ───────────────────────────────────────────────────────

def test_convert_state():
    random.seed(42)
    C = CrossJoin(5, '12')
    C.generateNonlinearity()
    state = [1, 0, 1, 1, 0]
    converted = C.convert_state(state)
    assert len(converted) == 5
    assert all(b in (0, 1) for b in converted)


# ── blocks property ─────────────────────────────────────────────────────

def test_blocks():
    C = CrossJoin(7, '69')
    assert C.blocks == [list(range(7))]


# ── Root expressions / monomial profiles ────────────────────────────────

def test_root_expressions():
    random.seed(42)
    C = CrossJoin(7, '69')
    C.generateNonlinearity()
    REs = C.root_expressions()
    assert len(REs) == 7

def test_monomial_profiles():
    random.seed(42)
    C = CrossJoin(7, '69')
    C.generateNonlinearity()
    MPs = C.monomial_profiles()
    assert len(MPs) == 7


# ── Linear complexity ──────────────────────────────────────────────────

def test_crossjoin_higher_lc_than_lfsr():
    """CrossJoin should have higher linear complexity than its base LFSR."""
    random.seed(42)
    C = CrossJoin(7, '69')
    C.generateNonlinearity(maxAnds=3, tapDensity=0.75)

    reg = FeedbackRegister(2**7 - 1, C)
    seq = [state[0] for state in reg.run(compiled=False, limit=500)]
    lc_cj, _ = berlekamp_massey(seq)

    assert lc_cj >= 7, f"CrossJoin LC ({lc_cj}) should be at least base LFSR size (7)"

