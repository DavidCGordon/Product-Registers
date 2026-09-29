"""Tests for annihilator computation: Gaussian and sparse variants.

Both implementations return `(degrees, [g, ...])` where `degrees` is the pair
`(deg g, deg h)` for `h = f*g` -- see Components_Architecture.md section 3.
Only when the second entry is 0 is `g` an annihilator in the strict sense of
`f*g == 0`; otherwise `h` is a nonzero low-degree multiple, which is the
object a fast algebraic attack consumes.  The first sections below use
functions whose pair has a zero multiple degree, so `AND(ann, f)` does vanish
there; the final section tests the general contract, including the case where
it does not.

Both algorithms are exercised on the same inputs because they are independent
implementations of the same specification -- Gaussian goes through
EquationStores and GF(2) RREF, Sparse is a self-contained symbolic search --
so a disagreement localizes a bug to one of them.

Small functions (≤ 4 variables) are used so the tests run quickly.
"""
import contextlib
import io
import random
from itertools import product

import pytest

from PyPR.BooleanLogic import AND, CONST, VAR, XOR, BooleanANF

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
    assert isinstance(degrees, tuple)
    assert len(degrees) == 2
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
    (ann_deg, _mult_deg), anns = gaussian_annihilators(f)
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


# ── The (annihilator, multiple) pair contract ─────────────────────────────────

ALGORITHMS = [("gaussian", gaussian_annihilators), ("sparse", sparse_annihilators)]


def random_xor_of_ands(rng, n):
    """A random XOR of three AND terms over n variables, each of degree 1-3."""
    return XOR(*[
        AND(*[VAR(j) for j in rng.sample(range(n), rng.randint(1, 3))])
        for _ in range(3)
    ])


def search_quietly(algorithm, f):
    """Run an annihilator search with its progress reporting suppressed.

    Both implementations print a live status block to stdout, which pytest
    captures but which makes -s output unreadable; `verbose=False` silences most
    of it and the redirect catches the rest.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        return algorithm(f, verbose=False)


@pytest.mark.parametrize(("name", "algorithm"), ALGORITHMS, ids=[n for n, _ in ALGORITHMS])
@pytest.mark.parametrize("trial", range(6))
def test_returned_degrees_bound_the_pair_they_describe(name, algorithm, trial):
    """`degrees` is (deg g, deg f*g), and both bounds are respected.

    The first entry bounds the annihilator's own degree; the second bounds the
    degree of the product.  If either were reported for the wrong member of the
    pair, an attack sizing its monomial basis from `degrees` would allocate the
    wrong number of unknowns.
    """
    rng = random.Random(100 + trial)
    n = 4
    f = random_xor_of_ands(rng, n)
    (annihilator_degree, multiple_degree), found = search_quietly(algorithm, f)

    assert found, f"{name} trial {trial}: no annihilator returned"
    for g in found:
        assert BooleanANF.from_BooleanFunction(g).degree() <= annihilator_degree, (
            f"{name} trial {trial}: annihilator exceeds its reported degree"
        )
        assert BooleanANF.from_BooleanFunction(AND(g, f)).degree() <= multiple_degree, (
            f"{name} trial {trial}: product f*g exceeds its reported degree"
        )


@pytest.mark.parametrize(("name", "algorithm"), ALGORITHMS, ids=[n for n, _ in ALGORITHMS])
@pytest.mark.parametrize("trial", range(6))
def test_product_vanishes_exactly_when_the_multiple_degree_is_zero(name, algorithm, trial):
    """f*g == 0 on all inputs if and only if `degrees[1] == 0`.

    This is the distinction between a true annihilator and a low-degree pair.
    A zero-degree multiple is the zero function, so the pair degenerates to the
    annihilator case; a positive-degree multiple means `g` does *not* annihilate
    `f`, and any caller treating it as one is computing with a nonzero residue.
    """
    rng = random.Random(200 + trial)
    n = 4
    f = random_xor_of_ands(rng, n)
    (_, multiple_degree), found = search_quietly(algorithm, f)

    vanishes = all(
        AND(g, f).eval([(i >> k) & 1 for k in range(n)]) == 0
        for g in found for i in range(2 ** n)
    )
    assert vanishes == (multiple_degree == 0), (
        f"{name} trial {trial}: multiple degree is {multiple_degree} but "
        f"f*g {'vanishes' if vanishes else 'does not vanish'}"
    )


@pytest.mark.parametrize("trial", range(8))
def test_both_algorithms_find_the_same_optimal_degree_pair(trial):
    """Gaussian and Sparse agree on `degrees` for the same input.

    The degree pair is a property of the function, not of the search strategy:
    it is the minimum over all valid pairs.  Two independent searches that
    disagree mean one of them is not finding the optimum.
    """
    rng = random.Random(300 + trial)
    f = random_xor_of_ands(rng, 4)
    assert search_quietly(gaussian_annihilators, f)[0] == search_quietly(sparse_annihilators, f)[0], (
        f"trial {trial}: the two annihilator algorithms disagree on the degree pair"
    )


@pytest.mark.parametrize(("name", "algorithm"), ALGORITHMS, ids=[n for n, _ in ALGORITHMS])
@pytest.mark.parametrize("trial", range(6))
def test_annihilator_degree_never_exceeds_the_function_degree(name, algorithm, trial):
    """`degrees[0] <= deg(f)`, because 1+f always annihilates f.

    Over GF(2) every function is idempotent (x^2 = x pointwise), so
    f*(1+f) = f + f^2 = 0.  That makes 1+f an annihilator of degree deg(f),
    which caps the optimum the search is allowed to report.  A larger reported
    degree would mean the search missed a candidate it is guaranteed to have.
    """
    rng = random.Random(400 + trial)
    n = 4
    f = random_xor_of_ands(rng, n)
    function_degree = BooleanANF.from_BooleanFunction(f).degree()

    complement = XOR(f, CONST(1))
    for i in range(2 ** n):
        state = [(i >> k) & 1 for k in range(n)]
        assert AND(complement, f).eval(state) == 0, "f*(1+f) must vanish over GF(2)"

    (annihilator_degree, _), _ = search_quietly(algorithm, f)
    assert annihilator_degree <= function_degree, (
        f"{name} trial {trial}: reported degree {annihilator_degree} exceeds "
        f"deg(f) = {function_degree}, but 1+f is always available"
    )


@pytest.mark.parametrize(("name", "algorithm"), ALGORITHMS, ids=[n for n, _ in ALGORITHMS])
@pytest.mark.parametrize("trial", range(6))
def test_returned_annihilators_are_nonzero(name, algorithm, trial):
    """No returned g is the zero function.

    The zero function annihilates everything, so admitting it would make the
    search trivially succeed at degree 0 and report a meaningless pair.  The
    definition requires g != 0.
    """
    rng = random.Random(500 + trial)
    f = random_xor_of_ands(rng, 4)
    _, found = search_quietly(algorithm, f)

    for g in found:
        assert list(BooleanANF.from_BooleanFunction(g)), (
            f"{name} trial {trial}: the zero function was returned as an annihilator"
        )
