"""Tests for AlgClosure.NumberTheory.

Importing this module used to start a 128-bit primitive-polynomial search at
module level and never return, so nothing in it could be tested -- or used --
without first slicing the file apart. That search now runs only as a script.
"""
from math import gcd

import galois
import pytest

import PyPR.Tools.AlgClosure.NumberTheory as number_theory

# x^5 + x^2 + 1, as coefficients [a0, a1, ..., a5] -- the library's convention
_BASE = [1, 0, 1, 0, 0, 1]
_DEGREE = len(_BASE) - 1


def test_importing_the_module_does_not_run_its_exploratory_search():
    """The search bound `base` and `x` at module level; they must not exist."""
    assert not hasattr(number_theory, "base")
    assert not hasattr(number_theory, "x")


@pytest.mark.parametrize(
    "d", [d for d in range(1, 2 ** _DEGREE - 1) if gcd(d, 2 ** _DEGREE - 1) == 1]
)
def test_decimation_gives_the_minimal_polynomial_of_alpha_to_the_d(d):
    """Decimating an m-sequence by d gives an m-sequence of the same degree.

    With alpha a root of the base polynomial, the base sequence is
    s_t = Tr(c * alpha^t), so the decimated one s_{dt} = Tr(c * (alpha^d)^t)
    has alpha^d as a characteristic root, and its connection polynomial is the
    minimal polynomial of alpha^d. Every d coprime to 2^n - 1 keeps alpha^d
    primitive, so that polynomial is primitive of degree n too.

    `galois` lists coefficients highest degree first, so its polynomials are
    reversed into the library's [a0, ..., an] before comparing -- and the
    comparison is against the minimal polynomial itself, not its reciprocal.
    """
    field = galois.GF(2 ** _DEGREE, irreducible_poly=galois.Poly(_BASE[::-1]))
    minimal = (field("x") ** d).minimal_poly()
    expected = [int(c) for c in minimal.coeffs[::-1]]

    decimated = number_theory.fast_decimate(_BASE, d)

    assert decimated == expected
    assert all(type(c) is int for c in decimated), "coefficients must be plain ints"
    assert galois.Poly(decimated[::-1]).is_primitive()


def test_decimation_by_one_is_the_identity():
    assert number_theory.fast_decimate(_BASE, 1) == _BASE


def test_decimation_composes():
    """Decimating by d then by e is decimating by d*e (mod 2^n - 1)."""
    order = 2 ** _DEGREE - 1
    twice = number_theory.fast_decimate(number_theory.fast_decimate(_BASE, 3), 5)
    assert twice == number_theory.fast_decimate(_BASE, (3 * 5) % order)
