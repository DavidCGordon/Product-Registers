"""Sanity tests for SAT solving: satisfiable, unsatisfiable, solution correctness,
functional equivalence, and model enumeration.

SAT support is patched onto BooleanFunction via PyPR.BooleanLogic.SAT, which is
imported automatically when PyPR.BooleanLogic is loaded.
"""
from PyPR.BooleanLogic import AND, OR, XOR, NOT, VAR, CONST


# ── satisfiable / unsatisfiable ──────────────────────────────────────────────

def test_sat_var_is_satisfiable():
    """VAR(0) is satisfiable — set x0 = True."""
    result = VAR(0).sat()
    assert result is not None
    assert result[0] is True

def test_sat_const_one_is_satisfiable():
    """CONST(1) is always true; sat() should return a (possibly empty) assignment."""
    result = CONST(1).sat()
    assert result is not None

def test_sat_const_zero_is_unsatisfiable():
    """CONST(0) is always false."""
    result = CONST(0).sat()
    assert result is None

def test_sat_contradiction_unsatisfiable():
    """AND(x0, NOT(x0)) is a contradiction."""
    result = AND(VAR(0), NOT(VAR(0))).sat()
    assert result is None

def test_sat_tautology_satisfiable():
    """OR(x0, NOT(x0)) is always true."""
    result = OR(VAR(0), NOT(VAR(0))).sat()
    assert result is not None

def test_sat_conjunction_requires_all_true():
    """AND(x0, x1, x2) is satisfiable only when all three variables are True."""
    fn = AND(VAR(0), VAR(1), VAR(2))
    sol = fn.sat()
    assert sol is not None
    assert sol.get(0) is True
    assert sol.get(1) is True
    assert sol.get(2) is True

def test_sat_solution_satisfies_function():
    """The assignment returned by sat() actually evaluates the function to 1."""
    fn = AND(VAR(0), OR(VAR(1), VAR(2)))
    sol = fn.sat()
    assert sol is not None
    # Build a state vector; unmentioned variables default to 0
    state = [int(sol.get(i, False)) for i in range(3)]
    assert fn.eval(state) == 1

def test_sat_solution_satisfies_xor():
    """XOR(x0, x1) solution correctly evaluates to 1."""
    fn = XOR(VAR(0), VAR(1))
    sol = fn.sat()
    assert sol is not None
    state = [int(sol.get(i, False)) for i in range(2)]
    assert fn.eval(state) == 1


# ── functional equivalence ───────────────────────────────────────────────────

def test_functionally_equivalent_self():
    """Every function is equivalent to itself."""
    fn = XOR(AND(VAR(0), VAR(1)), VAR(2))
    assert fn.functionally_equivalent(fn)

def test_functionally_equivalent_anf_round_trip():
    """A function converted to ANF and back is functionally equivalent to the original."""
    from PyPR.BooleanLogic import BooleanANF
    fn = XOR(AND(VAR(0), VAR(1)), VAR(2), CONST(1))
    fn2 = BooleanANF.from_BooleanFunction(fn).to_BooleanFunction()
    assert fn.functionally_equivalent(fn2)

def test_functionally_not_equivalent():
    """XOR(x0, x1) and AND(x0, x1) are not equivalent."""
    f = XOR(VAR(0), VAR(1))
    g = AND(VAR(0), VAR(1))
    assert not f.functionally_equivalent(g)


# ── model enumeration ────────────────────────────────────────────────────────

def test_enum_models_xor_two_solutions():
    """XOR(x0, x1) has exactly two satisfying assignments: (0,1) and (1,0)."""
    fn = XOR(VAR(0), VAR(1))
    models = list(fn.enum_models())
    assert len(models) == 2
    for sol in models:
        state = [int(sol.get(i, False)) for i in range(2)]
        assert fn.eval(state) == 1

def test_enum_models_and_one_solution():
    """AND(x0, x1) over two bits has exactly one satisfying assignment: (1,1)."""
    fn = AND(VAR(0), VAR(1))
    models = list(fn.enum_models())
    assert len(models) == 1
    assert models[0].get(0) is True
    assert models[0].get(1) is True

def test_enum_models_all_satisfy():
    """Every enumerated model actually satisfies the function."""
    fn = OR(VAR(0), VAR(1), VAR(2))
    for sol in fn.enum_models():
        state = [int(sol.get(i, False)) for i in range(3)]
        assert fn.eval(state) == 1
