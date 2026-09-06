"""Tests for BooleanGF: rational polynomials over GF(2).

BooleanGF represents elements of GF(2)(D) — the field of rational functions
in the formal delay operator D over GF(2).  It is used to express generating
functions for binary sequences.  Tests here verify the basic field laws
(identity, absorption, GF(2) additive self-cancellation) and the integration
with Berlekamp-Massey via from_seq.
"""
import numpy as np

from PyPR.BooleanLogic.BooleanGF import BooleanGF
from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR


# ── Helpers ───────────────────────────────────────────────────────────────────

def _eq(a: BooleanGF, b: BooleanGF) -> bool:
    """Equality after reducing both sides to lowest terms."""
    return a.simplify() == b.simplify()


# ── Named element constructors ────────────────────────────────────────────────

def test_one_has_unit_numerator_and_denominator():
    """BooleanGF.one() is the polynomial 1/1."""
    one = BooleanGF.one()
    assert str(one.num) == "1"
    assert str(one.den) == "1"

def test_zero_has_zero_numerator():
    """BooleanGF.zero() has a zero numerator."""
    zero = BooleanGF.zero()
    assert str(zero.num) == "0"

def test_delay_represents_d_operator():
    """BooleanGF.delay() is x/1 (the D operator), distinct from one() and zero()."""
    D = BooleanGF.delay()
    assert not _eq(D, BooleanGF.one())
    assert not _eq(D, BooleanGF.zero())


# ── Additive identity ─────────────────────────────────────────────────────────

def test_add_zero_is_identity():
    """a + zero() simplifies to a for a nontrivial element."""
    D = BooleanGF.delay()
    assert _eq(D + BooleanGF.zero(), D)

def test_add_zero_is_identity_for_one():
    """one() + zero() == one()."""
    assert _eq(BooleanGF.one() + BooleanGF.zero(), BooleanGF.one())


# ── GF(2) self-cancellation ───────────────────────────────────────────────────

def test_add_self_is_zero():
    """In GF(2), a + a = 0 for any element."""
    D = BooleanGF.delay()
    assert _eq(D + D, BooleanGF.zero())

def test_add_one_self_is_zero():
    """one() + one() = zero() in GF(2)."""
    assert _eq(BooleanGF.one() + BooleanGF.one(), BooleanGF.zero())


# ── Multiplicative identity ───────────────────────────────────────────────────

def test_mul_one_is_identity():
    """a * one() == a."""
    D = BooleanGF.delay()
    assert _eq(D * BooleanGF.one(), D)

def test_mul_one_on_arbitrary():
    """Multiplying a rational element by one() leaves it unchanged."""
    a = BooleanGF([1, 1], [1, 0, 1])  # (1+D) / (1+D^2)
    assert _eq(a * BooleanGF.one(), a)


# ── Zero multiplication ───────────────────────────────────────────────────────

def test_mul_zero_is_absorbing():
    """a * zero() == zero() for any element."""
    D = BooleanGF.delay()
    assert _eq(D * BooleanGF.zero(), BooleanGF.zero())


# ── Power ─────────────────────────────────────────────────────────────────────

def test_power_zero_is_one():
    """Any element to the power 0 is one()."""
    D = BooleanGF.delay()
    assert _eq(D ** 0, BooleanGF.one())

def test_power_one_is_self():
    """Any element to the power 1 equals itself."""
    D = BooleanGF.delay()
    assert _eq(D ** 1, D)

def test_power_two_is_square():
    """D^2 * D^2 == D^4 (sequential squaring)."""
    D = BooleanGF.delay()
    assert _eq(D ** 4, (D ** 2) * (D ** 2))


# ── Division round-trip ───────────────────────────────────────────────────────

def test_div_mul_round_trip():
    """(a / b) * b == a after simplification."""
    a = BooleanGF([1, 1], [1])   # 1 + D
    b = BooleanGF([1, 0, 1], [1])  # 1 + D^2
    assert _eq((a / b) * b, a)

def test_one_over_d_times_d_is_one():
    """(1/D) * D == 1 in GF(2)(D)."""
    D = BooleanGF.delay()
    inv_D = BooleanGF.one() / D
    assert _eq(inv_D * D, BooleanGF.one())


# ── Simplify ─────────────────────────────────────────────────────────────────

def test_simplify_produces_coprime_fraction():
    """After simplify(), gcd(num, den) divides 1 (i.e., both are coprime)."""
    import galois as gl
    a = BooleanGF([1, 0, 1], [1])  # 1 + D^2
    b = BooleanGF([1, 1], [1])     # 1 + D
    # (1+D^2) * (1+D) / (1+D) * (1+D) — denominator shares factor (1+D)
    product = a * b
    reduced = product.simplify()
    g = gl.gcd(reduced.num, reduced.den)
    assert str(g) == "1", f"GCD of simplified fraction should be 1, got {g}"


# ── from_seq (Berlekamp-Massey integration) ───────────────────────────────────

def test_from_seq_on_mpr_sequence_has_correct_denominator():
    """from_seq on an MPR bit-0 sequence recovers a degree-n denominator.

    The MPR (shift-up convention) produces a sequence whose connection
    polynomial is the reciprocal of the primitive polynomial P.  from_seq
    uses Berlekamp-Massey, so the denominator degree must equal n and the
    denominator must be the reciprocal of P (i.e. reversed coefficient list).
    """
    M = MPR(3, "5")  # degree 3, primitive poly [1, 0, 1, 1] (low-degree first)
    reg = FeedbackRegister(2**3 - 1, M)
    seq = [state[0] for state in reg.run(20, compiled=False)]

    gf = BooleanGF.from_seq(seq)

    # Denominator coefficients (low-degree first)
    den_coeffs = [int(c) for c in reversed(gf.den.coefficients())]

    # MPR 'shift-up' makes BM recover the reciprocal polynomial
    expected = M.primitive_polynomial[::-1]
    assert len(den_coeffs) == len(expected), (
        f"Denominator length {len(den_coeffs)} != expected {len(expected)}"
    )
    assert den_coeffs == expected, (
        f"Denominator {den_coeffs} != reciprocal of primitive polynomial {expected}"
    )

def test_from_seq_numerator_is_low_degree():
    """The numerator from from_seq has degree strictly less than the denominator.

    Note: galois Poly.degree is a property (int), not a callable.
    """
    M = MPR(5, "12")
    reg = FeedbackRegister(2**5 - 1, M)
    seq = [state[0] for state in reg.run(30, compiled=False)]

    gf = BooleanGF.from_seq(seq)
    # galois Poly.degree is a property
    num_deg = gf.num.degree
    den_deg = gf.den.degree
    assert num_deg < den_deg, (
        f"Numerator degree {num_deg} should be < denominator degree {den_deg}"
    )
