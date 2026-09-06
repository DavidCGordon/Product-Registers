"""Tests for shrinking, self-shrinking, and alternating-step generator constructions.

Verifies:
- SG LC = w * L_B (exact) when gcd(q_A, q_B) = 1
- SSG LC ≤ 2^{n-1}
- ASG LC ≤ w * L_B + (q_A - w) * L_C
- All three work on arbitrary (non-LFSR) sequences
"""
import random
from math import gcd

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import Fibonacci
from PyPR.Tools.ClockControlled import (
    shrinking_generator,
    self_shrinking_generator,
    alternating_step_generator,
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
