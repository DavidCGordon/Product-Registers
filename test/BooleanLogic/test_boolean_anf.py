"""Tests for BooleanANF: algebraic structure, identity elements, and round-trip conversion.

BooleanANF represents functions in Algebraic Normal Form as a frozenset of
frozensets (monomials).  Because ANF is a canonical form, equality of
BooleanANF objects implies functional equivalence.  The tests here check
that the algebraic ring laws hold and that conversion to/from BooleanFunction
preserves semantics.
"""
from itertools import product as iter_product

from PyPR.BooleanLogic import AND, XOR, VAR, CONST
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
