"""Decimation and expansion of periodic binary sequences, measured by Berlekamp-Massey.

Theory: docs/theory/Decimation Expansion and Linear Complexity.md
"""
import itertools
from math import gcd

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey


def test_odd_decimation_coprime_to_period_preserves_lc():
    """Decimation by an odd d coprime to the period leaves the LC unchanged.

    No two roots collide when gcd(d, period) = 1, and for odd d the leading
    coefficient d^r of binom(dt, r) in the basis binom(t, i) is odd, so each
    Jordan block keeps its size. Checked on every sequence of least period up
    to 9, which includes even periods and so Jordan sizes above 1.
    """
    for period in range(1, 10):
        for bits in itertools.product([0, 1], repeat=period):
            least = next(p for p in range(1, period + 1)
                         if period % p == 0 and all(bits[i] == bits[i % p] for i in range(period)))
            if not any(bits) or least != period:
                continue
            s = list(bits) * 120
            lc, _ = berlekamp_massey(s[:120])
            for d in (3, 5, 7):
                if gcd(d, period) == 1:
                    assert berlekamp_massey(s[::d][:120])[0] == lc


def test_decimation_sharing_a_factor_with_the_period_can_preserve_lc():
    """gcd(d, period) > 1 allows the LC to drop but does not force it.

    The m-sequence of period 15 has roots with exponents 1, 2, 4, 8 over 15.
    Decimating by 3 sends them to 3, 6, 12, 9, four distinct roots of order 5:
    the period drops to 5 and the LC stays 4.
    """
    m_sequence = [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1] * 20
    assert berlekamp_massey(m_sequence)[0] == 4
    decimated = m_sequence[::3]
    assert decimated[:5] * 4 == decimated[:20]
    assert berlekamp_massey(decimated)[0] == 4


def test_colliding_roots_add_their_coefficients():
    """When roots collide, the merged root survives iff the coefficients sum to nonzero.

    Both period-3 sequences below have the two roots of order 3, which both
    map to 1 under decimation by 3. One decimates to all-ones (LC 1), the
    other to zero (LC 0): the number of colliding roots does not decide it.
    """
    assert berlekamp_massey(([1, 1, 0] * 40)[::3])[0] == 1
    assert berlekamp_massey(([0, 1, 1] * 40)[::3])[0] == 0


def test_even_decimation_jordan_size_is_a_bound():
    """Decimation by 2 gives Jordan size at most ceil(m / 2), not always equal.

    s[t] = t mod 2 has the single root 1 with Jordan size 2 (LC 2), and
    s[2t] = 0: the size drops to 0, below ceil(2 / 2) = 1.
    """
    s = [0, 1] * 40
    assert berlekamp_massey(s)[0] == 2
    assert berlekamp_massey(s[::2])[0] == 0


def test_expansion_multiplies_lc():
    """Inserting d - 1 zeros between terms of a purely periodic s multiplies LC by d.

    If S(D) = N(D) / M(D) is in lowest terms, so is N(D^d) / M(D^d), and its
    denominator has degree d * deg M. This holds for every d, even or odd.
    """
    for period in range(1, 8):
        for bits in itertools.product([0, 1], repeat=period):
            if not any(bits):
                continue
            s = list(bits) * 40
            lc, _ = berlekamp_massey(s[:80])
            for d in (2, 3, 4, 5):
                expanded = [0] * (len(s) * d)
                expanded[::d] = s
                assert berlekamp_massey(expanded[:2 * d * lc + 8 * d])[0] == d * lc
