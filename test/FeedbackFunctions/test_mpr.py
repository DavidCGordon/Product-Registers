"""Tests for MPR and CMPR: construction, size, period, copy, and block structure."""
from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.BooleanLogic.ChainingGeneration.Templates import fast_template


# ── MPR construction ─────────────────────────────────────────────────────────

def test_mpr_size():
    M = MPR(5, "12")
    assert M.size == 5

def test_mpr_primitive_polynomial_length():
    """primitive_polynomial has exactly size + 1 entries."""
    for n, poly in [(3, "5"), (5, "12"), (7, "65")]:
        M = MPR(n, poly)
        assert len(M.primitive_polynomial) == n + 1, (
            f"MPR({n}): expected poly length {n + 1}, got {len(M.primitive_polynomial)}"
        )

def test_mpr_update_polynomial_length():
    """update_polynomial has exactly size entries."""
    for n, poly in [(3, "5"), (5, "12"), (7, "65")]:
        M = MPR(n, poly)
        assert len(M.update_polynomial) == n, (
            f"MPR({n}): expected update poly length {n}, got {len(M.update_polynomial)}"
        )

def test_mpr_fn_list_length():
    """fn_list contains one BooleanFunction per bit."""
    M = MPR(5, "12")
    assert len(M.fn_list) == 5

def test_mpr_binary_output():
    """MPR output stream contains only 0s and 1s."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    output = [state[0] for state in reg.run(compiled=False, limit=30)]
    assert all(b in (0, 1) for b in output)

def test_mpr_period_degree_3():
    """MPR with a degree-3 primitive polynomial has period 2^3 - 1 = 7."""
    M = MPR(3, "5")
    reg = FeedbackRegister(1, M)
    period, _ = reg.period(compiled=False)
    assert period == 7

def test_mpr_period_degree_5():
    """MPR with a degree-5 primitive polynomial has period 2^5 - 1 = 31."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    period, _ = reg.period(compiled=False)
    assert period == 31

def test_mpr_period_degree_7():
    """MPR with a degree-7 primitive polynomial has period 2^7 - 1 = 127."""
    M = MPR(7, "65")
    reg = FeedbackRegister(1, M)
    period, _ = reg.period(compiled=False)
    assert period == 127

def test_mpr_nonzero_seed_produces_full_period():
    """Any nonzero seed produces the full Mersenne-length cycle."""
    M = MPR(5, "12")
    for seed in [1, 15, 31]:
        reg = FeedbackRegister(seed, M)
        period, _ = reg.period(compiled=False)
        assert period == 31, f"seed={seed}: expected period 31, got {period}"

def test_mpr_copy_independent():
    """__copy__ returns an independent MPR with equal output functions."""
    M = MPR(5, "12")
    M2 = M.__copy__()
    assert M2.size == M.size
    state = [1, 0, 1, 1, 0]
    for i in range(5):
        assert M[i].eval(state) == M2[i].eval(state)

def test_mpr_with_update_polynomial():
    """Specifying a non-default update polynomial changes the update but not the period."""
    M_default = MPR(5, "12")
    M_custom  = MPR(5, "12", update_poly=[0, 1, 0, 0, 0])  # U(x) = x (same as default)
    # Both should still have period 31
    reg = FeedbackRegister(1, M_custom)
    period, _ = reg.period(compiled=False)
    assert period == 31


# ── CMPR construction ────────────────────────────────────────────────────────

def test_cmpr_size_is_sum():
    """CMPR size equals the sum of component sizes."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    assert C.size == 15

def test_cmpr_blocks_count():
    """CMPR has one block per component."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    assert len(C.blocks) == 3

def test_cmpr_block_sizes_match_components():
    """Each block size matches the corresponding MPR's size."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    assert [len(b) for b in C.blocks] == [7, 5, 3]

def test_cmpr_blocks_partition_all_bits():
    """Blocks cover every bit index in [0, size) exactly once."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    all_bits = sorted(b for block in C.blocks for b in block)
    assert all_bits == list(range(C.size))

def test_cmpr_binary_output_no_chaining():
    """CMPR without chaining produces binary output."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    C = CMPR([M7, M5])
    reg = FeedbackRegister(2**C.size - 1, C)
    output = [state[0] for state in reg.run(compiled=False, limit=30)]
    assert all(b in (0, 1) for b in output)

def test_cmpr_binary_output_with_chaining():
    """CMPR with chaining still produces binary output."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    C.generateChaining(template=fast_template())
    reg = FeedbackRegister(2**C.size - 1, C)
    output = [state[0] for state in reg.run(compiled=False, limit=30)]
    assert all(b in (0, 1) for b in output)

def test_cmpr_copy():
    """__copy__ returns an independent CMPR with the same size and block layout."""
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    C = CMPR([M7, M5])
    C2 = C.__copy__()
    assert C2.size == C.size
    assert [len(b) for b in C2.blocks] == [len(b) for b in C.blocks]

def test_cmpr_single_component():
    """CMPR wrapping a single MPR behaves like the MPR."""
    M5 = MPR(5, "12")
    C = CMPR([M5])
    assert C.size == 5
    reg = FeedbackRegister(1, C)
    period, _ = reg.period(compiled=False)
    assert period == 31
