"""Tests for MPR and CMPR: construction, period, block structure, and field arithmetic.

The second half checks that an MPR's bit-level simulation really is arithmetic
in GF(2^n).  An MPR identifies its n-bit state with an element of
GF(2^n) = GF(2)[x]/P(x) and clocks by multiplying by the update polynomial U
modulo P, but it never performs that multiplication at runtime: `MPR.__init__`
unrolls it once into per-bit boolean functions.  Those tests recompute the
field operation independently, so a failure means the unrolled functions no
longer implement multiplication mod P -- not that the ANF is built differently.

Required reading: docs/conventions/Polynomial Conventions.md (MPR section).
"""
import numpy as np
import pytest

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
    result = reg.period(compiled=False)
    assert result is not None, "period() found no cycle"
    period, _ = result
    assert period == 7

def test_mpr_period_degree_5():
    """MPR with a degree-5 primitive polynomial has period 2^5 - 1 = 31."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    result = reg.period(compiled=False)
    assert result is not None, "period() found no cycle"
    period, _ = result
    assert period == 31

def test_mpr_period_degree_7():
    """MPR with a degree-7 primitive polynomial has period 2^7 - 1 = 127."""
    M = MPR(7, "65")
    reg = FeedbackRegister(1, M)
    result = reg.period(compiled=False)
    assert result is not None, "period() found no cycle"
    period, _ = result
    assert period == 127

def test_mpr_nonzero_seed_produces_full_period():
    """Any nonzero seed produces the full Mersenne-length cycle."""
    M = MPR(5, "12")
    for seed in [1, 15, 31]:
        reg = FeedbackRegister(seed, M)
        result = reg.period(compiled=False)
        assert result is not None, "period() found no cycle"
        period, _ = result
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
    result = reg.period(compiled=False)
    assert result is not None, "period() found no cycle"
    period, _ = result
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
    result = reg.period(compiled=False)
    assert result is not None, "period() found no cycle"
    period, _ = result
    assert period == 31


# ── MPR as arithmetic in GF(2^n) ──────────────────────────────────────────────

PRIMITIVE = {
    3: [1, 1, 0, 1],              # 1 + x + x^3
    4: [1, 1, 0, 0, 1],           # 1 + x + x^4
    5: [1, 0, 1, 0, 0, 1],        # 1 + x^2 + x^5
    6: [1, 1, 0, 0, 0, 0, 1],     # 1 + x + x^6
}


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_clock_is_multiplication_by_update_poly_mod_P(n):
    """One clock cycle sends state s to U*s mod P in GF(2)[x]/P(x).

    This is the defining property of the register.  The comparison is against
    schoolbook GF(2) polynomial arithmetic written out here, so it is
    independent of how MPR chooses to unroll the multiplication into ANF.
    """
    P = PRIMITIVE[n]

    def polymul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    out[i + j] ^= bj
        return out

    def polymod(a, p):
        # Reduce degree-high terms downward using x^deg(p) = p(x) - x^deg(p).
        a, deg = list(a), len(p) - 1
        for i in range(len(a) - 1, deg - 1, -1):
            if a[i]:
                for j in range(deg + 1):
                    a[i - deg + j] ^= p[j]
        return (a[:deg] + [0] * deg)[:deg]

    fn = MPR(n, P)
    U = fn.update_polynomial

    for seed in range(2 ** n):
        reg = FeedbackRegister(seed, fn)
        # state[i] is the coefficient of x^i, matching the seed's bit order
        before = [int(b) for b in reg._state]
        reg.clock(compiled=False)
        after = [int(b) for b in reg._state]
        assert after == polymod(polymul(U, before), P), (
            f"n={n} seed={seed}: clocked to {after}, "
            f"but U*s mod P = {polymod(polymul(U, before), P)}"
        )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_P_tail_represents_x_inverse_mod_P(n):
    """P[1:] is the coefficient list of x^{-1} in GF(2)[x]/P(x).

    Polynomial Conventions.md builds the shift-down (Galois-equivalent) MPR as
    `MPR(n, P, update_poly=P[1:])`, justified by: P is primitive so P(0)=1, and
    writing P(x) = 1 + x*Q(x) gives x*Q(x) = P(x) + 1 == 1 (mod P).  Q is
    exactly P[1:].  This test checks that identity directly.
    """
    P = PRIMITIVE[n]

    def polymul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    out[i + j] ^= bj
        return out

    def polymod(a, p):
        a, deg = list(a), len(p) - 1
        for i in range(len(a) - 1, deg - 1, -1):
            if a[i]:
                for j in range(deg + 1):
                    a[i - deg + j] ^= p[j]
        return (a[:deg] + [0] * deg)[:deg]

    assert polymod(polymul([0, 1], P[1:]), P) == [1] + [0] * (n - 1), (
        f"x * P[1:] should be 1 mod P for n={n}"
    )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_update_poly_P_tail_inverts_the_default_register(n):
    """MPR(n, P, P[1:]) clocks the exact inverse transition of MPR(n, P).

    The default register multiplies by x; the P[1:] variant multiplies by
    x^{-1}.  Composing them must be the identity map on every state, including
    the all-zero state, which is a fixed point of both.
    """
    P = PRIMITIVE[n]
    forward = MPR(n, P)
    backward = MPR(n, P, update_poly=P[1:])

    for seed in range(2 ** n):
        f_reg = FeedbackRegister(seed, forward)
        f_reg.clock(compiled=False)
        b_reg = FeedbackRegister([int(b) for b in f_reg._state], backward)
        b_reg.clock(compiled=False)
        original = [int(x) for x in format(seed, f"0{n}b")[::-1]]
        assert [int(b) for b in b_reg._state] == original, (
            f"n={n}: clocking forward then backward did not return state {seed}"
        )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_orbit_of_nonzero_state_is_every_nonzero_state(n):
    """With a primitive P and U(x)=x, the nonzero states form one 2^n-1 cycle.

    Multiplication by x is multiplication by a generator of the cyclic group
    GF(2^n)*, so iterating it from any nonzero element visits all 2^n-1 of
    them.  Checking set equality (not just the count) also rules out an orbit
    that revisits states while still totalling 2^n-1 yielded values.
    """
    reg = FeedbackRegister(1, MPR(n, PRIMITIVE[n]))
    visited = {int(state) for state in reg.run(compiled=False, limit=2 ** n - 1)}
    assert visited == set(range(1, 2 ** n)), (
        f"n={n}: orbit covered {len(visited)} of {2 ** n - 1} nonzero states"
    )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_zero_state_is_a_fixed_point(n):
    """The zero element of the field is fixed: 0 * U = 0, so the orbit is trivial.

    This is the boundary case that separates "primitive polynomial" from
    "permutation of all 2^n states" -- the register is a bijection, but the
    long cycle covers only the 2^n-1 nonzero states.
    """
    reg = FeedbackRegister(0, MPR(n, PRIMITIVE[n]))
    reg.clock(compiled=False)
    assert int(reg) == 0


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_primitive_polynomials_are_closed_under_reversal(n):
    """Reversing a primitive polynomial yields another primitive polynomial.

    Reversal inverts the roots, and the inverse of a generator of GF(2^n)* is
    again a generator -- so the reversed polynomial is primitive too.
    Primitivity is checked here the direct way: build an MPR on the reversed
    polynomial and confirm its nonzero states form one full 2^n - 1 cycle.
    """
    reversed_poly = PRIMITIVE[n][::-1]
    assert reversed_poly[0] == 1 and reversed_poly[-1] == 1, (
        "a primitive polynomial has nonzero constant and leading terms, so its "
        "reversal is still a degree-n polynomial"
    )

    fn = MPR(n, reversed_poly)
    fn.compile()
    result = FeedbackRegister(1, fn).period(limit=2 ** 16)
    assert result is not None, "period() found no cycle"
    period, preperiod = result
    assert (period, preperiod) == (2 ** n - 1, 0), (
        f"n={n}: reversal of a primitive polynomial gave period {period}"
    )


# ── CMPR block update matrices ────────────────────────────────────────────────

@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [4, 5]])
def test_cmpr_block_update_matrices_reproduce_the_clock(sizes):
    """Each CMPR block's update matrix drives that block's bits, in isolation.

    With no chaining installed, the component feedback is the whole feedback,
    so `update_matrices[b] @ block_state` must equal the block's next state.
    Blocks are indexed high-to-low (see Notation and Terminology), and each
    matrix is offset to its own block, so this also pins the indexing.
    """
    fn = CMPR([MPR(k, PRIMITIVE[k]) for k in sizes])
    matrices, blocks = fn.update_matrices, fn.blocks

    for seed in range(1, 40):
        reg = FeedbackRegister(seed, fn)
        before = [int(b) for b in reg._state]
        reg.clock(compiled=False)
        after = [int(b) for b in reg._state]
        for index, block in enumerate(blocks):
            offset = min(block)
            current = np.array([before[i] for i in sorted(block)], dtype=np.uint8)
            expected = (matrices[index] @ current) % 2
            actual = np.array([after[i] for i in sorted(block)], dtype=np.uint8)
            assert np.array_equal(actual, expected), (
                f"sizes={sizes} seed={seed} block {index} (offset {offset}): "
                f"got {actual}, expected {expected}"
            )
