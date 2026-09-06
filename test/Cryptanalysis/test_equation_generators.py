"""Tests for all three equation generators: Symbolic, Substitution, and Cube.

All three generators express the register output at time t as a function of the
initial state variables.  They differ in representation:
  - SymbolicEqGenerator  → BooleanFunction (DAG, evaluated via .eval)
  - SubstitutionEqGenerator → BooleanFunction in ANF form (via .translate_ANF)
  - CubeEqGenerator      → numpy coefficient vector indexed by a var_map

For SymbolicEqGenerator, the tests verify:
  - The generator yields exactly limit+1 equations (t = 0 … limit).
  - At t=0 the symbolic expression, when evaluated at the initial state,
    matches the output function evaluated at the same state.
  - For every t, evaluating the symbolic equation at the initial state
    agrees with the actual register output observed at time t.

SubstitutionEqGenerator tests mirror those for SymbolicEqGenerator, and an
additional cross-check verifies both generators agree on all evaluations.

CubeEqGenerator tests verify output shape, binary values, and that evaluating
the returned coefficient vector as an ANF polynomial at the initial state
matches the actual register output.  Note: CubeEqGenerator yields `limit`
equations (range(limit)), not `limit+1` like the other two generators.
"""
import numpy as np

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR
from PyPR.BooleanLogic import VAR, AND, XOR, CONST
from PyPR.Cryptanalysis.Components.EquationGenerators.SymbolicEqGenerator import (
    SymbolicEqGenerator,
)
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import (
    SubstitutionEqGenerator,
)
from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import (
    CubeEqGenerator,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

M3 = MPR(3, "5")   # 3-bit MPR: fast, period 7


# ── SymbolicEqGenerator — equation count ──────────────────────────────────────

def test_generator_yields_limit_plus_one_equations():
    """SymbolicEqGenerator(feedback, output, limit) yields exactly limit+1 equations."""
    output_fn = VAR(0)
    limit = 6
    eqs = list(SymbolicEqGenerator(M3, output_fn, limit))
    assert len(eqs) == limit + 1, (
        f"Expected {limit + 1} equations, got {len(eqs)}"
    )

def test_generator_yields_limit_plus_one_for_const_output():
    """Generator respects limit even for a constant output function."""
    eqs = list(SymbolicEqGenerator(M3, CONST(1), limit=4))
    assert len(eqs) == 5


# ── SymbolicEqGenerator — t=0 correctness ────────────────────────────────────

def test_first_equation_equals_output_at_initial_state():
    """At t=0 the symbolic equation evaluates to output_fn(initial_state).

    With no initialization override, the t=0 symbolic state IS the initial
    state, so the first equation should evaluate identically to output_fn.
    """
    output_fn = VAR(0)
    seed = 0b101  # binary 5
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]

    gen = SymbolicEqGenerator(M3, output_fn, limit=0)
    eq_t0 = next(iter(gen))

    symbolic_val = eq_t0.eval(initial_state)
    direct_val   = output_fn.eval(initial_state)
    assert symbolic_val == direct_val, (
        f"t=0 symbolic equation gave {symbolic_val}, direct eval gave {direct_val}"
    )

def test_first_equation_with_and_output():
    """t=0 equation matches AND(x0, x1) evaluated directly."""
    output_fn = AND(VAR(0), VAR(1))
    seed = 0b110
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]

    gen = SymbolicEqGenerator(M3, output_fn, limit=0)
    eq_t0 = next(iter(gen))

    assert eq_t0.eval(initial_state) == output_fn.eval(initial_state)


# ── SymbolicEqGenerator — multi-step correctness ─────────────────────────────

def test_equations_match_register_output():
    """Each symbolic equation evaluates to the actual register output at that time.

    We run the register from a fixed seed to get the ground-truth output
    sequence, then check that evaluating the t-th symbolic equation at the
    initial state matches that output.
    """
    output_fn = VAR(0)
    seed = 0b101  # initial state = [1, 0, 1] (LSB first)
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]

    limit = 6
    reg = FeedbackRegister(seed, M3.__copy__())
    true_outputs = [output_fn.eval(list(state)) for state in reg.run(limit + 1, compiled=False)]

    gen = SymbolicEqGenerator(M3.__copy__(), output_fn, limit=limit)
    for t, eq in enumerate(gen):
        sym_val = eq.eval(initial_state)
        assert sym_val == true_outputs[t], (
            f"At t={t}: symbolic={sym_val}, actual={true_outputs[t]}"
        )

def test_equations_match_register_output_with_xor_output():
    """Equation correctness holds for a non-trivial XOR output function."""
    output_fn = XOR(VAR(0), VAR(1))
    seed = 0b011
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]

    limit = 5
    reg = FeedbackRegister(seed, M3.__copy__())
    true_outputs = [output_fn.eval(list(state)) for state in reg.run(limit + 1, compiled=False)]

    gen = SymbolicEqGenerator(M3.__copy__(), output_fn, limit=limit)
    for t, eq in enumerate(gen):
        sym_val = eq.eval(initial_state)
        assert sym_val == true_outputs[t], (
            f"At t={t}: symbolic={sym_val}, actual={true_outputs[t]}"
        )


# ── SymbolicEqGenerator — list output functions ───────────────────────────────

def test_list_output_returns_list_of_equations():
    """Passing a list of output functions yields a list of equations per step."""
    output_fns = [VAR(0), VAR(1)]
    gen = SymbolicEqGenerator(M3.__copy__(), output_fns, limit=3)
    for item in gen:
        assert isinstance(item, list), f"Expected list, got {type(item)}"
        assert len(item) == 2, f"Expected 2 equations, got {len(item)}"

def test_list_output_equations_are_correct():
    """Each equation in a list-output generator evaluates correctly."""
    output_fns = [VAR(0), VAR(1)]
    seed = 0b101
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]
    limit = 4

    reg = FeedbackRegister(seed, M3.__copy__())
    true_states = [list(state) for state in reg.run(limit + 1, compiled=False)]

    gen = SymbolicEqGenerator(M3.__copy__(), output_fns, limit=limit)
    for t, eq_list in enumerate(gen):
        for fn_idx, (output_fn, eq) in enumerate(zip(output_fns, eq_list)):
            sym_val  = eq.eval(initial_state)
            true_val = output_fn.eval(true_states[t])
            assert sym_val == true_val, (
                f"fn {fn_idx} at t={t}: symbolic={sym_val}, actual={true_val}"
            )


# ── SubstitutionEqGenerator — equation count ──────────────────────────────────

def test_substitution_yields_limit_plus_one_equations():
    """SubstitutionEqGenerator yields exactly limit+1 equations, same as Symbolic."""
    limit = 6
    eqs = list(SubstitutionEqGenerator(M3.__copy__(), VAR(0), limit))
    assert len(eqs) == limit + 1, (
        f"Expected {limit + 1} equations, got {len(eqs)}"
    )


# ── SubstitutionEqGenerator — t=0 correctness ────────────────────────────────

def test_substitution_first_equation_matches_initial_output():
    """At t=0 the substitution equation evaluates to output_fn(initial_state).

    At t=0 no feedback has been applied, so the substitution generator's equation
    is just the output function composed with identity substitutions — equivalent
    to output_fn itself.
    """
    output_fn = VAR(0)
    seed = 0b101
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]

    eq_t0 = next(iter(SubstitutionEqGenerator(M3.__copy__(), output_fn, limit=0)))
    assert eq_t0.eval(initial_state) == output_fn.eval(initial_state), (
        f"t=0 substitution equation gave {eq_t0.eval(initial_state)}, "
        f"direct eval gave {output_fn.eval(initial_state)}"
    )


# ── SubstitutionEqGenerator — multi-step correctness ─────────────────────────

def test_substitution_equations_match_register_output():
    """Each substitution equation evaluates to the actual register output at that time."""
    output_fn = VAR(0)
    seed = 0b101
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]
    limit = 6

    reg = FeedbackRegister(seed, M3.__copy__())
    true_outputs = [output_fn.eval(list(state)) for state in reg.run(limit + 1, compiled=False)]

    gen = SubstitutionEqGenerator(M3.__copy__(), output_fn, limit)
    for t, eq in enumerate(gen):
        assert eq.eval(initial_state) == true_outputs[t], (
            f"At t={t}: substitution={eq.eval(initial_state)}, actual={true_outputs[t]}"
        )


# ── SubstitutionEqGenerator — agreement with Symbolic ────────────────────────

def test_substitution_agrees_with_symbolic_on_evaluation():
    """SubstitutionEqGenerator and SymbolicEqGenerator evaluate identically.

    Both generators express the output at time t as a function of the initial
    state, using different internal composition directions.  Their evaluations
    at any fixed initial state must agree for all t.
    """
    output_fn = VAR(0)
    seed = 0b110
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]
    limit = 5

    symbolic_eqs = list(SymbolicEqGenerator(M3.__copy__(), output_fn, limit=limit))
    substitution_eqs = list(SubstitutionEqGenerator(M3.__copy__(), output_fn, limit))

    for t, (sym_eq, sub_eq) in enumerate(zip(symbolic_eqs, substitution_eqs)):
        sym_val = sym_eq.eval(initial_state)
        sub_val = sub_eq.eval(initial_state)
        assert sym_val == sub_val, (
            f"At t={t}: symbolic={sym_val} != substitution={sub_val}"
        )


# ── SubstitutionEqGenerator — list output ────────────────────────────────────

def test_substitution_list_output_returns_list_of_equations():
    """Passing a list of output functions yields a list per step."""
    output_fns = [VAR(0), VAR(1)]
    gen = SubstitutionEqGenerator(M3.__copy__(), output_fns, limit=3)
    for item in gen:
        assert isinstance(item, list), f"Expected list, got {type(item)}"
        assert len(item) == 2


# ── CubeEqGenerator — setup ───────────────────────────────────────────────────
#
# CubeEqGenerator takes a var_map: a dict mapping monomial tuples (sorted
# variable index tuples) to contiguous integer indices.  The empty tuple ()
# represents the constant term and MUST be included: the P.I.E. decomposition
# always generates empty sub-cubes (the all-zeros evaluation state).
#
# For a 3-bit linear register, degree-1 monomials are sufficient to represent
# every output exactly — MPR is linear, so all time-evolved outputs are linear
# combinations of the initial bits.
#
# CubeEqGenerator yields `limit` equations (range(limit)), not `limit+1`.

_VAR_MAP_3 = {(): 0, (0,): 1, (1,): 2, (2,): 3}  # constant + degree-1 for n=3


# ── CubeEqGenerator — structural ─────────────────────────────────────────────

def test_cube_yields_limit_equations():
    """CubeEqGenerator yields exactly `limit` equations (not limit+1).

    Unlike Symbolic and Substitution generators (which use range(limit+1)),
    CubeEqGenerator's main loop is range(limit), yielding limit items.
    """
    limit = 5
    eqs = list(CubeEqGenerator(M3.__copy__(), VAR(0), limit, _VAR_MAP_3))
    assert len(eqs) == limit, f"Expected {limit} equations, got {len(eqs)}"

def test_cube_output_is_binary_numpy_array():
    """Each yielded item is a numpy array of dtype uint8 with values in {0, 1}."""
    eqs = list(CubeEqGenerator(M3.__copy__(), VAR(0), 3, _VAR_MAP_3))
    for t, eq in enumerate(eqs):
        assert isinstance(eq, np.ndarray), f"t={t}: expected ndarray, got {type(eq)}"
        assert eq.dtype == np.uint8, f"t={t}: expected uint8, got {eq.dtype}"
        assert eq.shape == (len(_VAR_MAP_3),), (
            f"t={t}: expected shape ({len(_VAR_MAP_3)},), got {eq.shape}"
        )
        assert set(eq.tolist()).issubset({0, 1}), f"t={t}: non-binary values {eq}"


# ── CubeEqGenerator — t=0 coefficient check ──────────────────────────────────

def test_cube_linear_output_coefficients_at_t0():
    """At t=0 the coefficient vector for VAR(0) has exactly x_0 set to 1, rest 0.

    The ANF coefficient of a monomial m is computed by the GF(2) finite
    difference (Möbius inversion): coeff[m] = XOR of output_fn(s) over all
    subsets s of the variables in m.  For output_fn = x_0:
      - constant ():  output([0,0,0]) = 0
      - monomial (0,): output([0,0,0]) XOR output([1,0,0]) = 0 XOR 1 = 1
      - monomial (1,): output([0,0,0]) XOR output([0,1,0]) = 0 XOR 0 = 0
      - monomial (2,): output([0,0,0]) XOR output([0,0,1]) = 0 XOR 0 = 0
    """
    eq_t0 = next(iter(CubeEqGenerator(M3.__copy__(), VAR(0), 1, _VAR_MAP_3)))
    assert int(eq_t0[_VAR_MAP_3[()]]) == 0,   "constant term should be 0"
    assert int(eq_t0[_VAR_MAP_3[(0,)]]) == 1, "x_0 coefficient should be 1"
    assert int(eq_t0[_VAR_MAP_3[(1,)]]) == 0, "x_1 coefficient should be 0"
    assert int(eq_t0[_VAR_MAP_3[(2,)]]) == 0, "x_2 coefficient should be 0"


# ── All generators agree ──────────────────────────────────────────────────────

def test_all_generators_agree_on_evaluation():
    """Symbolic, Substitution, and Cube generators evaluate identically at every t.

    All three generators express the same thing — the output at time t as a
    polynomial in the initial state variables — using different internal
    representations.  Evaluating each at the same initial state must give the
    same bit for all t in [0, limit).

    CubeEqGenerator yields `limit` items (range(limit)); the other two yield
    `limit+1` (range(limit+1)).  We compare over the shared prefix [0, limit).
    For the Cube generator, evaluation uses the ANF formula:
      output(t) = XOR over monomials m of ( c[var_map[m]] * AND(x_i for i in m) )
    MPR is linear, so the degree-0/1 var_map captures all terms exactly.
    """
    output_fn = VAR(0)
    seed = 0b110
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]
    limit = 5

    sym_gen  = SymbolicEqGenerator(M3.__copy__(), output_fn, limit=limit)
    sub_gen  = SubstitutionEqGenerator(M3.__copy__(), output_fn, limit)
    cube_gen = CubeEqGenerator(M3.__copy__(), output_fn, limit, _VAR_MAP_3)

    for t, (sym_eq, sub_eq, coeff_vec) in enumerate(zip(sym_gen, sub_gen, cube_gen)):
        sym_val = sym_eq.eval(initial_state)
        sub_val = sub_eq.eval(initial_state)

        # ANF polynomial evaluation for the cube coefficient vector
        cube_val = 0
        for comb, idx in _VAR_MAP_3.items():
            mon = 1
            for v in comb:
                mon &= initial_state[v]
            cube_val ^= int(coeff_vec[idx]) * mon

        assert sym_val == sub_val == cube_val, (
            f"At t={t}: symbolic={sym_val}, substitution={sub_val}, cube={cube_val}"
        )


# ── CubeEqGenerator — multi-step correctness ─────────────────────────────────

def test_cube_equations_match_register_output():
    """Evaluating the coefficient polynomial at the initial state matches register output.

    The coefficient vector c encodes the ANF of the output at time t:
      output(t) = XOR over monomials m of ( c[var_map[m]] * AND(x_i for i in m) )
    For a linear MPR, all outputs are degree-1 polynomials of the initial bits,
    so the degree-0/1 var_map here is exact (no higher terms are truncated).
    """
    output_fn = VAR(0)
    seed = 0b101
    # initial state as a list, LSB first
    initial_state = [int(b) for b in format(seed, f"0{M3.size}b")[::-1]]
    limit = 5

    reg = FeedbackRegister(seed, M3.__copy__())
    true_outputs = [output_fn.eval(list(state)) for state in reg.run(limit, compiled=False)]

    for t, coeff_vec in enumerate(CubeEqGenerator(M3.__copy__(), output_fn, limit, _VAR_MAP_3)):
        # Evaluate ANF polynomial: XOR over monomials of (coefficient * product of bits)
        poly_val = 0
        for comb, idx in _VAR_MAP_3.items():
            # Monomial value: 1 if all bits in comb are 1, else 0 (empty product = 1)
            mon = 1
            for v in comb:
                mon &= initial_state[v]
            poly_val ^= int(coeff_vec[idx]) * mon

        assert poly_val == true_outputs[t], (
            f"At t={t}: polynomial eval={poly_val}, actual={true_outputs[t]}"
        )
