"""Tests for BooleanFunction basics: CONST, VAR, gates, ANF construction, and ANF translation.

Derived from the ProductRegisters.ipynb notebook (Library Basics > Boolean Functions).
"""
from itertools import product

from PyPR.BooleanLogic import (
    BooleanFunction, BooleanANF,
    CONST, VAR, AND, OR, XOR, NOT, NAND, NOR, XNOR,
)


# -- Inputs ------------------------------------------------------------------

def test_const():
    state = [1, 0, 1, 1, 1]
    assert CONST(0).eval(state) == 0
    assert CONST(1).eval(state) == 1

def test_var():
    state = [1, 0, 1, 1, 1]
    expected = [1, 0, 1, 1, 1]
    for i in range(5):
        assert VAR(i).eval(state) == expected[i], f"VAR({i}) failed"


# -- Gates -------------------------------------------------------------------

def test_and():
    for a, b in product([0, 1], repeat=2):
        assert AND(CONST(a), CONST(b)).eval([]) == (a & b)

def test_or():
    for a, b in product([0, 1], repeat=2):
        assert OR(CONST(a), CONST(b)).eval([]) == (a | b)

def test_xor():
    for a, b in product([0, 1], repeat=2):
        assert XOR(CONST(a), CONST(b)).eval([]) == (a ^ b)

def test_not():
    assert NOT(CONST(0)).eval([]) == 1
    assert NOT(CONST(1)).eval([]) == 0

def test_nand():
    for a, b in product([0, 1], repeat=2):
        assert NAND(CONST(a), CONST(b)).eval([]) == (1 - (a & b))

def test_nor():
    for a, b in product([0, 1], repeat=2):
        assert NOR(CONST(a), CONST(b)).eval([]) == (1 - (a | b))

def test_xnor():
    for a, b in product([0, 1], repeat=2):
        assert XNOR(CONST(a), CONST(b)).eval([]) == (1 - (a ^ b))

def test_composite_gate():
    """Composite gate from notebook cell 21."""
    fn = XNOR(
        NAND(VAR(0), VAR(1)),
        NOR(VAR(2), VAR(3)),
        AND(VAR(1), VAR(2), CONST(1))
    )
    for inputs in product([0, 1], repeat=4):
        result = fn.eval(list(inputs))
        assert result in (0, 1), f"Unexpected result {result} for inputs {inputs}"


# -- ANF Construction --------------------------------------------------------

def test_anf_from_ANF():
    """BooleanFunction.from_ANF from notebook cell 26."""
    fn = BooleanFunction.from_ANF([[1], [1, 2], [2, 4, 5], [1, 3, 4], True])
    for inputs in product([0, 1], repeat=6):
        result = fn.eval(list(inputs))
        assert result in (0, 1)

def test_anf_translate_zero_function():
    """NOR(AND(1,2,3), CONST(1)) is the zero function (notebook cell 27)."""
    fn = NOR(AND(VAR(1), VAR(2), VAR(3)), CONST(1))
    for inputs in product([0, 1], repeat=4):
        assert fn.eval(list(inputs)) == 0, f"Expected 0 for {inputs}"

    anf = BooleanANF.from_BooleanFunction(fn)
    assert len(anf.terms) == 0, "ANF of zero function should have no terms"

def test_anf_translate_nontrivial():
    """NOR(AND(1,2,3), VAR(4)) has a nontrivial ANF (notebook cell 27)."""
    fn = NOR(AND(VAR(1), VAR(2), VAR(3)), VAR(4))
    anf = BooleanANF.from_BooleanFunction(fn)
    assert len(anf.terms) > 0, "ANF should have terms"

    anf_fn = anf.to_BooleanFunction()
    for inputs in product([0, 1], repeat=5):
        inp = list(inputs)
        assert fn.eval(inp) == anf_fn.eval(inp), f"Mismatch at {inputs}"

def test_anf_round_trip():
    """Converting to ANF and back should preserve the function."""
    fn = XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1))
    anf = BooleanANF.from_BooleanFunction(fn)
    fn2 = anf.to_BooleanFunction()
    for inputs in product([0, 1], repeat=3):
        inp = list(inputs)
        assert fn.eval(inp) == fn2.eval(inp), f"Round-trip mismatch at {inputs}"


# -- Printing ----------------------------------------------------------------

def test_dense_str():
    fn = XOR(VAR(0), VAR(1))
    s = fn.dense_str()
    assert isinstance(s, str) and len(s) > 0

def test_pretty_str():
    fn = AND(VAR(0), OR(VAR(1), VAR(2)))
    s = fn.pretty_str()
    assert isinstance(s, str) and len(s) > 0
