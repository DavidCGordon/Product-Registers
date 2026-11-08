"""Tests for shrinking, self-shrinking, and alternating-step generator constructions.

Verifies:
- SG LC ≤ w * L_B, attained at the tabulated parameters but not in general
- SSG LC ≤ 2^{n-1}, and below it at the tested degrees
- ASG LC ≤ q_A * (L_B + L_C), attained at the tested parameters
- interleaving a sequence with itself gives LC 2L or 2L - 1
- All three work on arbitrary (non-LFSR) sequences

Theory: docs/theory/Interleaving and Clock-Controlled Registers.md
"""
import itertools
import random
from math import gcd

from PyPR.FeedbackFunctions import Fibonacci
from PyPR.FeedbackRegister import FeedbackRegister

from PyPR.Tools.ClockControlled import (
    alternating_step_generator,
    self_shrinking_generator,
    shrinking_generator,
)
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey


def lfsr_sequence(degree, poly_hex, num_bits):
    F = Fibonacci(degree, poly_hex)
    reg = FeedbackRegister(1, F)
    return [state[0] for state in reg.run(compiled=False, limit=num_bits)]


# ── Shrinking Generator ──────────────────────────────────────────────────

def test_sg_basic():
    """SG with trivial control [1,1,1,...] returns the data unchanged."""
    data = [1, 0, 1, 1, 0, 0, 1]
    control = [1] * len(data)
    assert shrinking_generator(control, data) == data


def test_sg_lc_formula_coprime():
    """LC = w * L_B when gcd(q_A, q_B) = 1.

    Control: degree 3, period 7, w = 4 ones.
    Data: degree 4, period 15, L_B = 4.
    gcd(7, 15) = 1, so LC should be exactly 4 * 4 = 16.
    """
    n_a, n_b = 3, 4
    q_a, q_b = 2**n_a - 1, 2**n_b - 1
    assert gcd(q_a, q_b) == 1

    control = lfsr_sequence(n_a, '3', q_a * q_b * 4)
    data = lfsr_sequence(n_b, '3', q_a * q_b * 4)

    w = sum(control[:q_a])
    assert w == 2**(n_a - 1)

    output = shrinking_generator(control, data)
    lc, _ = berlekamp_massey(output)
    assert lc == w * n_b, f"Expected LC={w * n_b}, got {lc}"


def test_sg_lc_formula_multiple_polys():
    """Verify LC = w * L_B across several coprime (n_A, n_B) pairs."""
    cases = [
        (3, '3', 4, '3'),
        (3, '3', 5, '12'),
        (4, '3', 5, '12'),
        (4, '3', 3, '3'),
    ]
    for n_a, poly_a, n_b, poly_b in cases:
        q_a, q_b = 2**n_a - 1, 2**n_b - 1
        assert gcd(q_a, q_b) == 1, f"gcd({q_a},{q_b}) != 1"

        length = q_a * q_b * 4
        control = lfsr_sequence(n_a, poly_a, length)
        data = lfsr_sequence(n_b, poly_b, length)
        w = sum(control[:q_a])

        output = shrinking_generator(control, data)
        lc, _ = berlekamp_massey(output)
        assert lc == w * n_b, (
            f"n_A={n_a}, n_B={n_b}: expected LC={w * n_b}, got {lc}"
        )


def test_sg_arbitrary_sequences():
    """SG works on non-LFSR sequences and produces a valid output."""
    random.seed(42)
    control = [random.randint(0, 1) for _ in range(100)]
    data = [random.randint(0, 1) for _ in range(100)]

    output = shrinking_generator(control, data)
    assert len(output) == sum(control)
    for i, (c, d) in enumerate(zip(control, data)):
        if c == 1:
            assert output.pop(0) == d


# ── Self-Shrinking Generator ─────────────────────────────────────────────

def test_ssg_basic():
    """SSG pairs consecutive bits and keeps second when first is 1."""
    seq = [1, 0, 0, 1, 1, 1, 0, 0]
    # Pairs: (1,0), (0,1), (1,1), (0,0)
    # Keep second when first=1: 0, 1
    assert self_shrinking_generator(seq) == [0, 1]


def test_ssg_lc_bound():
    """SSG LC ≤ 2^{n-1} for m-sequence inputs."""
    for n, poly in [(3, '3'), (4, '3'), (5, '5')]:
        q = 2**n - 1
        seq = lfsr_sequence(n, poly, q * q * 4)
        output = self_shrinking_generator(seq)
        lc, _ = berlekamp_massey(output)
        bound = 2**(n - 1)
        assert lc <= bound, f"n={n}: LC={lc} exceeds bound {bound}"
        assert lc > 0, f"n={n}: LC should be positive"


def test_ssg_arbitrary():
    """SSG works on arbitrary sequences."""
    random.seed(99)
    seq = [random.randint(0, 1) for _ in range(200)]
    output = self_shrinking_generator(seq)
    expected_len = sum(seq[2 * t] for t in range(len(seq) // 2))
    assert len(output) == expected_len


# ── Alternating-Step Generator ────────────────────────────────────────────

def test_asg_basic():
    """ASG with all-ones control only advances B."""
    data_b = [1, 0, 1, 1]
    data_c = [0, 0, 0, 0]
    control = [1, 1, 1, 1]
    # B advances each step, C stays at index 0.
    # y[t] = data_b[t] ^ data_c[0]
    expected = [b ^ data_c[0] for b in data_b]
    assert alternating_step_generator(control, data_b, data_c) == expected


def test_asg_lc_bound():
    """ASG LC ≤ q_A * (L_B + L_C) for coprime LFSR inputs (Günther 1988).

    The sample-and-hold structure means each data register's contribution
    spans all q_A output positions, not just the w (or q_A-w) positions
    where that register advances.
    """
    n_a, n_b, n_c = 3, 4, 5
    q_a = 2**n_a - 1

    length = q_a * (2**n_b - 1) * (2**n_c - 1) * 2
    control = lfsr_sequence(n_a, '3', length)
    data_b = lfsr_sequence(n_b, '3', length)
    data_c = lfsr_sequence(n_c, '12', length)

    bound = q_a * (n_b + n_c)

    output = alternating_step_generator(control, data_b, data_c)
    lc, _ = berlekamp_massey(output)
    assert lc <= bound, f"LC={lc} exceeds bound {bound}"
    assert lc > 0


def test_asg_arbitrary():
    """ASG works on arbitrary sequences."""
    random.seed(77)
    control = [random.randint(0, 1) for _ in range(50)]
    data_b = [random.randint(0, 1) for _ in range(50)]
    data_c = [random.randint(0, 1) for _ in range(50)]

    output = alternating_step_generator(control, data_b, data_c)
    assert len(output) == len(control)


# ── Bounds that are not equalities ────────────────────────────────────────

def test_sg_bound_not_attained_for_general_control():
    """w * L_B is an upper bound: a non-LFSR control with coprime period falls short.

    Control 1110000000 has period 10 and w = 3; the data is an m-sequence of
    period 7, so gcd(10, 7) = 1 and the bound is 3 * 3 = 9. The selected times
    are 10k + j for j = 0, 1, 2, and 10k + j = 3k + j mod 7, which is the
    output index: the output is the data itself, of LC 3.
    """
    data = lfsr_sequence(3, '3', 700)
    control = ([1, 1, 1] + [0] * 7) * 70
    output = shrinking_generator(control, data)
    assert output == data[:len(output)]
    lc, _ = berlekamp_massey(output)
    assert lc == 3


def test_sg_bound_depends_on_polynomial_pair():
    """With m-sequence control and data, 2^{n_A-1} * n_B is not always attained.

    For n_A = 5, n_B = 3 the bound is 48. Whether it is attained depends on
    the pair of primitive polynomials and not on the phases: with coprime
    periods, shifting both registers together shifts the output, and the
    joint shift reaches every pair of phases.
    """
    n = 31 * 7 * 3
    attained = shrinking_generator(lfsr_sequence(5, '12', n), lfsr_sequence(3, '3', n))
    short = shrinking_generator(lfsr_sequence(5, '17', n), lfsr_sequence(3, '3', n))
    assert berlekamp_massey(attained[:140])[0] == 48
    assert berlekamp_massey(short[:140])[0] == 45


def test_ssg_below_bound():
    """The SSG bound 2^{n-1} is not attained at degrees 3, 4, 5."""
    for n, poly, expected in [(3, '3', 2), (4, '3', 5), (5, '5', 10)]:
        q = 2**n - 1
        output = self_shrinking_generator(lfsr_sequence(n, poly, q * q * 4))
        lc, _ = berlekamp_massey(output)
        assert lc == expected
        assert lc < 2**(n - 1)


def test_asg_attains_bound():
    """At (n_A, n_B, n_C) = (3, 4, 5) the ASG meets q_A * (L_B + L_C) = 63."""
    length = 7 * 15 * 31 * 2
    output = alternating_step_generator(
        lfsr_sequence(3, '3', length),
        lfsr_sequence(4, '3', length),
        lfsr_sequence(5, '12', length),
    )
    lc, _ = berlekamp_massey(output)
    assert lc == 63


def test_self_interleaving_lc():
    """Interleaving a with itself gives LC 2L(a) - 1 if 1 is a root of a, else 2L(a).

    Y(D) = (1 + D) A(D^2) = (1 + D) A(D)^2 over F_2. The factor 1 + D sits in
    the numerator, so it can cancel one factor of the squared denominator and
    can never add a root. 1 is a root of a exactly when 1 + D divides its
    minimal polynomial, i.e. when applying 1 + D to a lowers its LC.
    """
    for period in range(1, 8):
        for bits in itertools.product([0, 1], repeat=period):
            if not any(bits):
                continue
            a = list(bits) * (120 // period + 1)
            interleaved = [bit for bit in a for _ in range(2)]
            lc_a, _ = berlekamp_massey(a)
            # (1 + D) a: the termwise sum of a and its shift.
            differenced = [a[t] ^ a[t + 1] for t in range(len(a) - 1)]
            one_is_a_root = berlekamp_massey(differenced)[0] == lc_a - 1
            lc_y, _ = berlekamp_massey(interleaved)
            assert lc_y == (2 * lc_a - 1 if one_is_a_root else 2 * lc_a)
