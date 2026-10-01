"""End-to-end tests for the cube attack on a chained CMPR.

A cube attack treats some state bits as tweakable -- a public IV the attacker
sets freely -- and the rest as a secret key. Summing the output over every
assignment of a set of tweakable bits (a cube) leaves that cube's superpoly, a
polynomial in the key bits, and each keystream position gives one equation. The
offline phase takes cube candidates from the output's monomial profile, lowest
superpoly degree first, and recovers each superpoly exactly up to the degree the
profile allows; the online phase measures the same cube sums on the target and
hands the equations to an equation store and solver, as NAA, RAA and FAA do.

The attack reaches the target only through `access_fns`, which emulates what an
attacker of a real cipher has: `access_fn` runs the target with chosen IV bits,
`test_fn` says whether a candidate state reproduces its keystream, and
`sim_fn` runs the public design from any state. The secret stays inside, so a
test that recovers it cannot have been handed it.

The target has to be one the attack can break. A CMPR whose chaining runs only
from each block to the next -- which is what the shipped templates such as
`arman_template` build -- gives a cube over IV blocks nothing to say about a
separate key block: on the arman CMPRs checked while writing this (block sizes
7/5/3, 7/5/3/2 and 13/7/5/3, four seeds each), no cube candidate lies outside
every block its superpoly involves, at any degree. The degree a superpoly would
need has been traded into higher-degree terms of the next block (the trading is
described in `docs/architecture/Mesh Optimization.md`). So the tests use a
design built to be attacked, fixed rather than sampled: Arman's 17-bit
CMPR(M7, M5, M3, M2) from the cube-attack experiments, whose chaining functions
are written out by hand below and reach across blocks.

Blocks are [10..16] (M7), [5..9] (M5), [2..4] (M3) and [0, 1] (M2), and every
chaining function reads only from blocks before its own -- M7 drives the rest.
"""
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import VAR, BooleanFunction

from PyPR.FeedbackFunctions import CMPR, MPR
from PyPR.FeedbackRegister import FeedbackRegister

from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Cryptanalysis.Attacks.cube_attacks import (
    access_fns,
    cmpr_cube_attack_offline,
    cmpr_cube_summary,
    cube_attack_online,
)
from PyPR.Cryptanalysis.Components.EquationSolving.GaussElim import GaussElimSolver
from PyPR.Cryptanalysis.Components.EquationSolving.Grob_Solver import GrobnerSolver
from PyPR.Cryptanalysis.Components.EquationSolving.SplitGrob_Solver import (
    SplitGrobnerSolver,
)
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore

_INIT_ROUNDS = 100


def _c17():
    """Arman's 17-bit cube-attack target: 3-input XORs and one 4-input AND per
    chaining function, each reading only from earlier blocks."""
    M7 = MPR(7, [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0])
    M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
    M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
    M2 = MPR(2, [1, 1, 1], [1, 1])
    cmpr = CMPR([M7, M5, M3, M2])
    cmpr[0].add_arguments(BooleanFunction.from_ANF([True, [6], [2, 3, 7, 13]]))
    cmpr[1].add_arguments(BooleanFunction.from_ANF([[2], [3], [4, 9, 11, 14]]))
    cmpr[2].add_arguments(BooleanFunction.from_ANF([True, [10], [5, 7, 11, 15]]))
    cmpr[3].add_arguments(BooleanFunction.from_ANF([[5], [7], [8, 9, 14, 15]]))
    cmpr[4].add_arguments(BooleanFunction.from_ANF([[11], [13], [6, 7, 10, 16]]))
    cmpr[5].add_arguments(BooleanFunction.from_ANF([True, [14], [10, 11, 12, 13]]))
    cmpr[7].add_arguments(BooleanFunction.from_ANF([[10], [15], [11, 12, 13, 14]]))
    cmpr[8].add_arguments(BooleanFunction.from_ANF([True, [10], [11, 12, 14, 16]]))
    cmpr[9].add_arguments(BooleanFunction.from_ANF([[11], [12], [13, 14, 15, 16]]))
    return cmpr


# Every store + solver pairing RAA and FAA accept, as factories so each case
# gets fresh instances.
_ONLINE_COMBINATIONS = {
    "defaults-LU": dict,
    "GaussElim": lambda: {"solver": GaussElimSolver()},
    "Groebner-store-and-solver": lambda: {
        "solver": GrobnerSolver(), "online_store": GroebnerEqStore(simplify_mode=None),
    },
    "Groebner-store-split-solver": lambda: {
        "solver": SplitGrobnerSolver(), "online_store": GroebnerEqStore(simplify_mode=None),
    },
}


# ── the candidates the attack is built from ──────────────────────────────────

def test_a_candidate_is_bounded_by_every_term_that_contains_it():
    """The superpoly of a cube collects every output term that contains it.

    For the cube <0:6/7> (six of M7's seven bits), the term <0:7/7, 1:4/5>
    contains it with one more variable in block 0 and four in block 1, so its
    superpoly can reach degree 5 across those blocks. The rule used to compare
    only the blocks the candidate itself lists, missed that term, and called
    the cube linear in block 0 -- the equations it then produced were wrong.
    """
    profile = VAR(0).translate_ANF().remap_constants([
        (0, MonomialProfile.logical_zero()),
        (1, MonomialProfile.logical_one()),
    ]).eval_ANF(_c17().monomial_profiles())
    candidates = {str(candidate): (targets, degree)
                  for candidate, targets, _num_cubes, degree in profile.get_cube_candidates()}

    assert candidates["<0:6/7>"] == ((0, 1, 2), 5)

def test_candidates_with_the_same_sizes_and_counts_are_all_kept():
    """<0:4/7, 1:5/5> and <0:5/7, 1:4/5> have the same multisets of block sizes
    and counts. Candidates were deduplicated by those multisets, so whichever came
    second was dropped; they are different cubes and both must survive."""
    profile = VAR(0).translate_ANF().remap_constants([
        (0, MonomialProfile.logical_zero()),
        (1, MonomialProfile.logical_one()),
    ]).eval_ANF(_c17().monomial_profiles())
    candidates = {str(candidate): degree
                  for candidate, _targets, _num_cubes, degree in profile.get_cube_candidates()}

    assert candidates["<0:4/7, 1:5/5>"] == 1
    assert candidates["<0:5/7, 1:4/5>"] == 2

def test_the_summary_reports_every_candidate(capsys):
    """The summary unpacked a fourth field the candidates no longer carried and
    could not run at all."""
    cmpr_cube_summary(_c17(), VAR(0), list(range(10)), analyze_sources=True)
    report = capsys.readouterr().out

    assert "Summary Finished!" in report
    assert "Superpoly Degree (at most):  1" in report


# ── round trips ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("make_kwargs", list(_ONLINE_COMBINATIONS.values()),
                         ids=list(_ONLINE_COMBINATIONS))
@pytest.mark.parametrize(("tweakable", "expected_guesses"), [
    # The downstream blocks as IV and M7 as key, as in the experiment. Cubes over
    # bits 0-9 give a linear equation for every M7 bit. Along the way the
    # offline phase saturates the target region mid-cube and then skips the
    # candidates left over it -- the branch that used to read `loop_nums`
    # before it was bound.
    (list(range(10)), []),
    # Leaving M2 out of the IV makes bits 0 and 1 secret too. Cubes still solve
    # M7, but nothing yields an equation on M2, so the solver has to recover
    # those two bits by guessing, each guess checked through test_fn.
    (list(range(2, 10)), [0, 1]),
], ids=["cubes-only", "cubes-and-guessing"])
def test_linear_cubes_recover_the_secret_state(tweakable, expected_guesses, make_kwargs):
    # chaining here is fixed, but seeded anyway so nothing depends on order
    random.seed(0)
    np.random.seed(0)
    cmpr = _c17()
    output_fn = VAR(0)
    register = FeedbackRegister(12345, cmpr)
    secret = [int(bit) for bit in register._state]
    access_fn, sim_fn, test_fn = access_fns(register, output_fn, tweakable, init_rounds=_INIT_ROUNDS)
    # the IV is public, and both phases hold it at the target's values
    known_bits = {bit: secret[bit] for bit in tweakable}

    attack_data = cmpr_cube_attack_offline(
        cmpr, output_fn, sim_fn, tweakable, known_bits, max_degree=1, time_limit=60,
    )
    covered = {bit for _cube, _t, monomials, _c in attack_data.equations
               for monomial in monomials for bit in monomial}
    guesses = [bit for bit in range(cmpr.size) if bit not in known_bits and bit not in covered]
    assert guesses == expected_guesses, (
        "the cubes no longer cover the key bits this case was written for"
    )

    recovered = cube_attack_online(cmpr, output_fn, access_fn, test_fn, attack_data, **make_kwargs())

    assert recovered is not None, "the cube attack returned no solution"
    assert [int(bit) for bit in recovered] == secret

@pytest.mark.slow
def test_nonlinear_cubes_recover_the_secret_state():
    """Tweaking M7 instead: every cube over the driving block has a nonlinear
    superpoly, so a linear-only attack finds nothing -- and the old rule, which
    mistook some for linear, produced wrong equations. With the degree bound the
    superpolys are recovered as they are and every solver takes them.
    """
    random.seed(0)
    np.random.seed(0)
    cmpr = _c17()
    output_fn = VAR(0)
    tweakable = list(range(10, 17))
    register = FeedbackRegister(12345, cmpr)
    secret = [int(bit) for bit in register._state]
    access_fn, sim_fn, test_fn = access_fns(register, output_fn, tweakable, init_rounds=_INIT_ROUNDS)
    known_bits = {bit: secret[bit] for bit in tweakable}

    attack_data = cmpr_cube_attack_offline(
        cmpr, output_fn, sim_fn, tweakable, known_bits, time_limit=120,
    )
    degrees = {len(monomial) for _cube, _t, monomials, _c in attack_data.equations
               for monomial in monomials}
    assert max(degrees) > 1, "every equation was linear; this case no longer tests degree"

    for name, make_kwargs in _ONLINE_COMBINATIONS.items():
        recovered = cube_attack_online(cmpr, output_fn, access_fn, test_fn, attack_data, **make_kwargs())
        assert recovered is not None, f"{name} returned no solution"
        assert [int(bit) for bit in recovered] == secret, f"{name} recovered the wrong state"


@pytest.mark.parametrize(("tweakable", "max_degree"), [
    # observed on this design: cubes over M2 alone give nothing
    ([0, 1], None),
    # M7's superpolys are all of degree 3 or more, so a cap of 2 excludes them
    (list(range(10, 17)), 2),
], ids=["M2-only", "M7-capped-at-degree-2"])
def test_an_attack_without_equations_is_rejected(tweakable, max_degree):
    """With no cubes the attack is a brute force with extra steps; it says so."""
    random.seed(0)
    np.random.seed(0)
    cmpr = _c17()
    output_fn = VAR(0)
    register = FeedbackRegister(12345, cmpr)
    secret = [int(bit) for bit in register._state]
    access_fn, sim_fn, test_fn = access_fns(register, output_fn, tweakable, init_rounds=_INIT_ROUNDS)
    known_bits = {bit: secret[bit] for bit in tweakable}

    attack_data = cmpr_cube_attack_offline(
        cmpr, output_fn, sim_fn, tweakable, known_bits, max_degree=max_degree, time_limit=60,
    )
    assert not attack_data.equations

    with pytest.raises(ValueError, match="No cubes given"):
        cube_attack_online(cmpr, output_fn, access_fn, test_fn, attack_data)

def test_the_offline_phase_requires_a_value_for_every_tweakable_bit():
    """Non-cube tweakable bits are held at their values in both phases; without
    one, the superpolys the offline phase computes would not be the target's."""
    cmpr = _c17()
    register = FeedbackRegister(12345, cmpr)
    _access_fn, sim_fn, _test_fn = access_fns(register, VAR(0), [0, 1, 2], init_rounds=_INIT_ROUNDS)

    with pytest.raises(ValueError, match=r"tweakable bits \[2\] have no value"):
        cmpr_cube_attack_offline(cmpr, VAR(0), sim_fn, [0, 1, 2], {0: 0, 1: 1})
