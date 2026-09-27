"""Tests for BooleanANF: algebraic structure, identity elements, and round-trip conversion.

BooleanANF represents functions in Algebraic Normal Form as a frozenset of
frozensets (monomials).  Because ANF is a canonical form, equality of
BooleanANF objects implies functional equivalence.  The tests here check
that the algebraic ring laws hold and that conversion to/from BooleanFunction
preserves semantics.

The final section checks the stronger statement behind that canonicity: the
ANF is the image of the truth table under the Mobius (Reed-Muller) transform
on the subset lattice, with the coefficient of monomial S equal to the XOR of
f over every subset of S.  Those tests compute the transform independently
and compare term for term.
"""
import random
from itertools import product as iter_product

import pytest

from PyPR.BooleanLogic import AND, XOR, OR, NOT, VAR, CONST, BooleanFunction
from PyPR.BooleanLogic.BooleanANF import BooleanANF


# ── Fixtures ──────────────────────────────────────────────────────────────────

# x0 XOR x1*x2  (degree 2, 3 variables)
_FN = XOR(VAR(0), AND(VAR(1), VAR(2)))
_A = BooleanANF.from_BooleanFunction(_FN)         # x0 + x1*x2

_B = BooleanANF.from_BooleanFunction(VAR(0))       # x0
_C = BooleanANF.from_BooleanFunction(AND(VAR(1), VAR(2)))  # x1*x2

_ZERO = BooleanANF()                               # additive identity
_ONE  = BooleanANF([True])                         # multiplicative identity


# ── Identity elements ─────────────────────────────────────────────────────────

def test_xor_zero_is_identity():
    """XOR with the empty ANF (zero) leaves a function unchanged."""
    assert (_A ^ _ZERO) == _A

def test_and_one_is_identity():
    """AND with logical-1 ANF leaves a function unchanged."""
    assert (_A & _ONE) == _A

def test_and_zero_is_absorbing():
    """AND with the empty ANF (zero) produces the empty ANF."""
    assert (_A & _ZERO) == _ZERO


# ── Self-cancellation (GF(2) property) ───────────────────────────────────────

def test_xor_with_itself_is_zero():
    """In GF(2), a XOR a = 0."""
    assert (_A ^ _A) == _ZERO

def test_xor_with_itself_is_zero_for_single_var():
    """Self-XOR is zero for a single-variable ANF."""
    assert (_B ^ _B) == _ZERO


# ── Commutativity ─────────────────────────────────────────────────────────────

def test_xor_is_commutative():
    """XOR is commutative: a XOR b == b XOR a."""
    assert (_B ^ _C) == (_C ^ _B)

def test_and_is_commutative():
    """AND is commutative: a AND b == b AND a."""
    assert (_B & _C) == (_C & _B)


# ── Associativity ─────────────────────────────────────────────────────────────

def test_xor_is_associative():
    """XOR is associative: (a XOR b) XOR c == a XOR (b XOR c)."""
    d = BooleanANF.from_BooleanFunction(CONST(1))
    assert ((_A ^ _B) ^ d) == (_A ^ (_B ^ d))

def test_and_is_associative():
    """AND is associative: (a AND b) AND c == a AND (b AND c)."""
    x0 = BooleanANF.from_BooleanFunction(VAR(0))
    x1 = BooleanANF.from_BooleanFunction(VAR(1))
    x2 = BooleanANF.from_BooleanFunction(VAR(2))
    assert ((x0 & x1) & x2) == (x0 & (x1 & x2))


# ── Distributivity ───────────────────────────────────────────────────────────

def test_and_distributes_over_xor():
    """AND distributes over XOR: a AND (b XOR c) == (a AND b) XOR (a AND c)."""
    assert (_B & (_C ^ _ZERO)) == ((_B & _C) ^ (_B & _ZERO))

def test_distributivity_nontrivial():
    """x0 AND (x0 XOR x1*x2) == x0*x0 XOR x0*x1*x2 == x0 XOR x0*x1*x2 in GF(2)."""
    lhs = _B & _A
    rhs = (_B & _B) ^ (_B & _C)
    assert lhs == rhs


# ── Degree ────────────────────────────────────────────────────────────────────

def test_degree_of_constant_is_zero():
    """The constant-1 ANF has degree 0 (empty monomial = degree 0)."""
    assert _ONE.degree() == 0

def test_degree_of_empty_is_zero():
    """The empty (zero) ANF has degree 0 by the max-with-default convention."""
    assert _ZERO.degree() == 0

def test_degree_of_linear_function():
    """A single variable x0 has degree 1."""
    assert _B.degree() == 1

def test_degree_of_product():
    """x1 AND x2 has degree 2."""
    assert _C.degree() == 2

def test_and_degree_additive_for_disjoint_vars():
    """For disjoint variable sets, degree(a AND b) = degree(a) + degree(b)."""
    x0 = BooleanANF.from_BooleanFunction(VAR(0))
    x1 = BooleanANF.from_BooleanFunction(VAR(1))
    assert (x0 & x1).degree() == x0.degree() + x1.degree()


# ── Inversion ─────────────────────────────────────────────────────────────────

def test_inversion_xors_constant():
    """Inverting an ANF XORs in the constant-1 term."""
    assert (~_B) == (_B ^ _ONE)

def test_double_inversion_is_identity():
    """Inverting twice returns the original ANF."""
    assert (~~_A) == _A


# ── Membership ────────────────────────────────────────────────────────────────

def test_constant_term_present_after_inversion():
    """After inverting a non-constant ANF, the constant term (True) is present."""
    assert (True in (~_B))

def test_constant_term_absent_in_pure_variable():
    """A pure variable ANF has no constant term."""
    assert (True not in _B)

def test_variable_term_present():
    """The term {0} (i.e. x0) is present in the x0 ANF."""
    assert ([0] in _B)

def test_variable_term_absent_in_product():
    """x0 alone is not a term in x1*x2."""
    assert ([0] not in _C)

def test_false_term_contains_returns_true_by_convention():
    """False / 0 passed to __contains__ always returns True by API convention.

    _convert_iterable_term maps False/0 to None, and __contains__ returns True
    for None inputs.  This is by design: the zero element is trivially 'present'
    in the sense that the query is vacuous.
    """
    assert (False in _A)
    assert (0 in _A)


# ── Equality and hash ─────────────────────────────────────────────────────────

def test_equal_anfs_have_equal_hash():
    """If a == b, then hash(a) == hash(b) (required by Python contract)."""
    a_copy = BooleanANF(_A.terms, fast_init=True)
    assert _A == a_copy
    assert hash(_A) == hash(a_copy)

def test_unequal_anfs_are_not_equal():
    """Distinct ANFs are not equal."""
    assert _B != _C


# ── Round-trip through BooleanFunction ───────────────────────────────────────

def test_round_trip_is_functionally_equivalent():
    """from_BooleanFunction(anf.to_BooleanFunction()) recovers the same ANF."""
    recovered = BooleanANF.from_BooleanFunction(_A.to_BooleanFunction())
    assert recovered == _A

def test_round_trip_preserves_degree():
    """Degree is preserved after ANF → BooleanFunction → ANF round-trip."""
    fn = _A.to_BooleanFunction()
    recovered = BooleanANF.from_BooleanFunction(fn)
    assert recovered.degree() == _A.degree()

def test_round_trip_agrees_on_all_inputs():
    """The round-trip function evaluates identically to the original on all inputs."""
    fn_orig  = _FN
    fn_round = BooleanANF.from_BooleanFunction(_FN).to_BooleanFunction()
    for bits in iter_product([0, 1], repeat=3):
        inp = list(bits)
        assert fn_orig.eval(inp) == fn_round.eval(inp), (
            f"Mismatch at input {inp}"
        )


# ── ANF as the Mobius transform of the truth table ────────────────────────────

def random_xor_of_ands(rng, n, terms=3, max_degree=3):
    """A random XOR of AND terms over n variables, as an arbitrary test function."""
    return XOR(*[
        AND(*[VAR(j) for j in rng.sample(range(n), rng.randint(1, max_degree))])
        for _ in range(terms)
    ])


@pytest.mark.parametrize("trial", range(8))
def test_anf_terms_are_the_mobius_transform_of_the_truth_table(trial):
    """The ANF's monomials are exactly the subsets with nonzero Mobius coefficient.

    The transform is computed here by the standard in-place butterfly over the
    n bit positions, which is the subset-XOR formula evaluated for every S at
    once.  Comparing the resulting index sets against the library's terms pins
    both the coefficients and the monomial indexing.
    """
    rng = random.Random(1000 + trial)
    n = 4
    f = random_xor_of_ands(rng, n)

    truth_table = [f.eval([(i >> k) & 1 for k in range(n)]) for i in range(2 ** n)]

    coefficients = list(truth_table)
    for k in range(n):
        for i in range(2 ** n):
            if i & (1 << k):
                coefficients[i] ^= coefficients[i ^ (1 << k)]

    expected = {
        frozenset(k for k in range(n) if (i >> k) & 1)
        for i in range(2 ** n) if coefficients[i]
    }
    assert set(BooleanANF.from_BooleanFunction(f)) == expected, (
        f"trial {trial}: ANF terms disagree with the Mobius transform"
    )


@pytest.mark.parametrize("trial", range(8))
def test_anf_reconstructs_the_original_truth_table(trial):
    """Evaluating the ANF reproduces f on every input -- the transform is invertible.

    Mobius inversion over GF(2) is an involution, so rebuilding a function from
    its ANF monomials must return the same truth table.  This is the round trip
    the uniqueness claim rests on, checked exhaustively rather than at sampled
    points.
    """
    rng = random.Random(2000 + trial)
    n = 4
    f = random_xor_of_ands(rng, n)

    terms = [sorted(term) for term in BooleanANF.from_BooleanFunction(f)]
    rebuilt = BooleanFunction.from_ANF(terms) if terms else CONST(0)

    for i in range(2 ** n):
        state = [(i >> k) & 1 for k in range(n)]
        assert rebuilt.eval(state) == f.eval(state), (
            f"trial {trial}: ANF reconstruction differs at input {state}"
        )


@pytest.mark.parametrize("trial", range(8))
def test_eval_and_eval_ANF_agree_on_every_input(trial):
    """`f.eval` and `f.translate_ANF().eval` agree on all 2^n inputs.

    `translate_ANF` rewrites the DAG into a XOR-of-ANDs; it is a structural
    rewrite that must preserve semantics exactly.  Exhaustive comparison is
    cheap at n = 4 and catches sign/ordering errors that sampled inputs miss.
    """
    rng = random.Random(3000 + trial)
    n = 4
    f = XOR(AND(VAR(0), VAR(1)), OR(VAR(2), NOT(VAR(3))), random_xor_of_ands(rng, n, terms=2))
    translated = f.translate_ANF()

    for i in range(2 ** n):
        state = [(i >> k) & 1 for k in range(n)]
        assert translated.eval(state) == f.eval(state), (
            f"trial {trial}: translate_ANF changed the value at {state}"
        )


@pytest.mark.parametrize("trial", range(8))
def test_degree_is_the_largest_monomial_size_and_is_bounded_by_n(trial):
    """`degree()` equals the largest ANF monomial's size, and never exceeds n.

    Over GF(2) every variable satisfies x^2 = x, so no monomial can repeat a
    variable and the degree is capped at the number of variables.  A degree
    above n would mean the ANF had stopped reducing x^2 to x.
    """
    rng = random.Random(4000 + trial)
    n = 5
    f = random_xor_of_ands(rng, n, terms=4, max_degree=4)

    anf = BooleanANF.from_BooleanFunction(f)
    terms = list(anf)
    expected = max((len(term) for term in terms), default=0)

    assert anf.degree() == expected, "ANF degree should be the largest monomial size"
    assert f.degree() == expected, "BooleanFunction.degree should agree with its ANF"
    assert anf.degree() <= n, "degree cannot exceed the number of variables over GF(2)"


@pytest.mark.parametrize("trial", range(8))
def test_xor_of_functions_is_symmetric_difference_of_anf_terms(trial):
    """ANF addition is GF(2) addition: terms cancel in pairs.

    Adding two functions adds their coefficient vectors mod 2, so a monomial
    present in both drops out.  This is what makes the ANF a linear-algebraic
    object and is the step annihilator search relies on.
    """
    rng = random.Random(5000 + trial)
    n = 4
    f, g = random_xor_of_ands(rng, n), random_xor_of_ands(rng, n)

    f_terms = set(BooleanANF.from_BooleanFunction(f))
    g_terms = set(BooleanANF.from_BooleanFunction(g))
    combined = set(BooleanANF.from_BooleanFunction(XOR(f, g)))

    assert combined == f_terms ^ g_terms, (
        f"trial {trial}: XOR of ANFs is not the symmetric difference of their terms"
    )


def test_constant_and_zero_functions_have_the_expected_anf():
    """The constant 1 is the empty monomial; the zero function has no monomials.

    The empty set indexes the degree-0 term (the empty product is 1), so
    CONST(1) must carry exactly `frozenset()` and CONST(0) must carry nothing.
    Getting this backwards would silently shift every function by a constant.
    """
    assert set(BooleanANF.from_BooleanFunction(CONST(1))) == {frozenset()}
    assert set(BooleanANF.from_BooleanFunction(CONST(0))) == set()
    # negation adds the constant term, since NOT(x) = 1 + x over GF(2)
    assert set(BooleanANF.from_BooleanFunction(NOT(VAR(0)))) == {frozenset(), frozenset({0})}
