"""Sanity tests for annihilator computation: Gaussian and sparse variants.

For each annihilator ann found by either algorithm, AND(ann, f) must be
identically zero on all inputs.  Both algorithms should agree on the
optimal (ann_degree, mult_degree) pair.

Small functions (≤ 4 variables, degree ≤ 2) are used so the tests run quickly.
"""
from itertools import product

from PyPR.BooleanLogic import AND, XOR, OR, VAR, CONST
from PyPR.Cryptanalysis.Components.Annihilators.GaussianAnnihilator import (
    annihilators as gaussian_annihilators,
)
from PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator import (
    annihilators as sparse_annihilators,
)


# ── Gaussian annihilators ────────────────────────────────────────────────────

def test_gaussian_annihilators_found():
    """Gaussian algorithm finds at least one annihilator for a degree-2 function."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    _, anns = gaussian_annihilators(f)
    assert len(anns) > 0

def test_gaussian_annihilators_correct():
    """Every Gaussian annihilator actually annihilates the function.

    By definition, ann is an annihilator of f iff AND(ann, f) ≡ 0 on all inputs
    (i.e. the product is identically the zero function over GF(2)).
    """
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    _, anns = gaussian_annihilators(f)
    for inputs in product([0, 1], repeat=3):
        inp = list(inputs)
        for ann in anns:
            assert AND(ann, f).eval(inp) == 0, (
                f"Annihilator {ann.anf_str()} failed to annihilate {f.anf_str()} "
                f"at input {inp}"
            )

def test_gaussian_degree_pair_is_tuple():
    """The returned degree pair is a length-2 tuple of non-negative integers."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    degrees, _ = gaussian_annihilators(f)
    assert isinstance(degrees, tuple) and len(degrees) == 2
    assert all(isinstance(d, int) and d >= 0 for d in degrees)


# ── Sparse annihilators ──────────────────────────────────────────────────────

def test_sparse_annihilators_found():
    """Sparse algorithm finds at least one annihilator."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    _, anns = sparse_annihilators(f)
    assert len(anns) > 0

def test_sparse_annihilators_correct():
    """Every sparse annihilator actually annihilates the function.

    By definition, ann is an annihilator of f iff AND(ann, f) ≡ 0 on all inputs
    (i.e. the product is identically the zero function over GF(2)).
    """
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    _, anns = sparse_annihilators(f)
    for inputs in product([0, 1], repeat=3):
        inp = list(inputs)
        for ann in anns:
            assert AND(ann, f).eval(inp) == 0, (
                f"Annihilator {ann.anf_str()} failed to annihilate {f.anf_str()} "
                f"at input {inp}"
            )


# ── Algorithm agreement ──────────────────────────────────────────────────────

def test_both_algorithms_same_degree_pair():
    """Gaussian and sparse algorithms report the same optimal degree pair."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))
    g_degrees, _ = gaussian_annihilators(f)
    s_degrees, _ = sparse_annihilators(f)
    assert g_degrees == s_degrees, (
        f"Degree mismatch: Gaussian={g_degrees}, Sparse={s_degrees}"
    )

# ── Theoretical properties ────────────────────────────────────────────────────

def test_annihilator_degree_at_most_function_degree():
    """For a degree-d function, the annihilator degree should be ≤ d."""
    f = XOR(AND(VAR(0), VAR(1)), VAR(2))  # degree 2
    (ann_deg, _), _ = gaussian_annihilators(f)
    assert ann_deg <= f.degree(), (
        f"Annihilator degree {ann_deg} exceeds function degree {f.degree()}"
    )

def test_linear_function_annihilated_by_degree_1():
    """A linear function is annihilated by functions of degree ≤ 1.

    By definition, ann is an annihilator of f iff AND(ann, f) ≡ 0 on all inputs.
    For a degree-1 (linear) function the algebraic immunity is 1, so the optimal
    annihilator has degree ≤ 1.
    """
    f = XOR(VAR(0), VAR(1), VAR(2))  # degree 1 (linear)
    (ann_deg, mult_deg), anns = gaussian_annihilators(f)
    assert len(anns) > 0
    for inputs in product([0, 1], repeat=3):
        inp = list(inputs)
        for ann in anns:
            assert AND(ann, f).eval(inp) == 0, (
                f"Annihilator {ann.anf_str()} failed to annihilate {f.anf_str()} "
                f"at input {inp}"
            )
    # For a linear function, optimal annihilator degree should be ≤ 1
    assert ann_deg <= 1
