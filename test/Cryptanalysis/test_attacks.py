"""End-to-end tests for the algebraic attacks: NAA, RAA and FAA.

Each attack runs in two phases. The offline phase inspects the register and the
output function and produces `attack_data` -- an equation system, index maps, and
the number of keystream bits it needs. The online phase takes a keystream
produced from a secret state and solves for that state.

The test that matters is therefore the round trip: seed a register, emit the
keystream the offline phase asked for, and require the online phase to return
the state we started from. Anything weaker -- that the phases run, that they
return the right shape -- would pass against an attack that recovers nothing.

Every attack is run across both of its offline paths and every online store and
solver `docs/architecture/Attack_Compatibilities.md` documents as compatible.
The register is deliberately small (a 5-bit and a 3-bit MPR, chained) so each
round trip takes a second or two; these attacks are exponential in the register
size. Wiring follows `.claude/commands/attack-setup.md`.
"""
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, VAR, XOR, BooleanANF
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template

from PyPR.FeedbackFunctions import CMPR, MPR
from PyPR.FeedbackRegister import FeedbackRegister

from PyPR.Cryptanalysis.Attacks.fast_algebraic_attack import FAA_offline, FAA_online
from PyPR.Cryptanalysis.Attacks.naive_algebraic_attack import NAA_offline, NAA_online
from PyPR.Cryptanalysis.Attacks.reduced_algebraic_attack import RAA_offline, RAA_online
from PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator import annihilators
from PyPR.Cryptanalysis.Components.EquationSolving.GaussElim import GaussElimSolver
from PyPR.Cryptanalysis.Components.EquationSolving.Grob_Solver import GrobnerSolver
from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.Cryptanalysis.Components.EquationSolving.SplitGrob_Solver import (
    SplitGrobnerSolver,
)
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore

_SECRET = 42
_TIME_LIMIT = 120
# Chaining generation draws from np.random, and FAA's dynamic path draws from
# random, so both are seeded and every run attacks the same register. Seed 1 is
# not arbitrary: its register leaves one variable undetermined after Groebner
# reduction, so the split solver genuinely splits rather than finding the system
# already solved -- some seeds (3 and 4, for instance) never enter that loop.
_REGISTER_SEED = 1


def _target():
    """A small chained CMPR, seeded, and a degree-2 output function f = x0 x1 + x2."""
    random.seed(_REGISTER_SEED)
    np.random.seed(_REGISTER_SEED)
    M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
    M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
    cmpr = CMPR([M5, M3])
    cmpr.generateChaining(template=arman_template(max_and=2))
    return cmpr, XOR(AND(VAR(0), VAR(1)), VAR(2))


def _keystream(cmpr, output_fn, num_bits):
    """The output bit of a register seeded with the secret, for `num_bits` clocks.

    A uint8 array: RAA and FAA XOR the keystream into a uint8 coefficient vector
    in place, which numpy refuses from a wider dtype such as int64.
    """
    reg = FeedbackRegister(_SECRET, cmpr)
    return np.array(
        [int(output_fn.eval(state._state)) for state in reg.run(num_bits, compiled=False)],
        dtype=np.uint8,
    )


def _secret_state(cmpr):
    return list(FeedbackRegister(_SECRET, cmpr)._state)


def _offline_path(cmpr, use_profiles):
    """Keyword arguments selecting an offline path.

    With monomial profiles the variable count is known ahead of time, equations
    come from `CubeEqGenerator`, and no rank tracker is built. Without them the
    attack discovers variables through `SubstitutionEqGenerator` and tracks rank
    to decide when it has enough equations -- and never compiles the feedback
    function, which is what exposed guess-and-solve's old assumption that it had.
    """
    if not use_profiles:
        return {}
    return {"monomial_profiles": cmpr.monomial_profiles(), "variable_blocks": cmpr.blocks}


_PATHS = pytest.mark.parametrize(
    "use_profiles", [True, False], ids=["profile-path", "dynamic-path"]
)


# Every documented online store + solver pairing for RAA and FAA, as factories so
# each case gets fresh instances.
_ONLINE_COMBINATIONS = [
    pytest.param(dict, id="defaults-LU"),
    pytest.param(lambda: {"solver": GaussElimSolver()}, id="GaussElim"),
    pytest.param(
        lambda: {"solver": GrobnerSolver(),
                 "online_store": GroebnerEqStore(simplify_mode=None)},
        id="Groebner-store-and-solver",
    ),
    pytest.param(
        lambda: {"solver": SplitGrobnerSolver(),
                 "online_store": GroebnerEqStore(simplify_mode=None)},
        id="Groebner-store-split-solver",
    ),
]


# ── the pairs the attacks are built on ───────────────────────────────────────
# Two (g, h = f*g) pairs are used below.
#
# The annihilator pair is the best one `annihilators` finds. For f = x0 x1 + x2
# no degree-1 annihilator exists, and every degree-1 g leaves h of degree 2, so
# the degree-minimising search settles on g = f + 1 at degrees (2, 0): a true
# annihilator, h = 0. That is exactly what RAA wants. For FAA it is the
# degenerate case -- with h = 0 there is no linear relation to exploit.
#
# The fast pair is what FAA is actually for: g = x0 of degree 1, and
# h = x0 x1 + x0 x2 nonzero. FAA's cost is set by deg g, because the linear
# relation of the sequence h(s_t) cancels h whatever its degree, so FAA prefers
# this pair over the degree-2 true annihilator.

def test_the_annihilator_pair_annihilates_the_output_function():
    """g is an annihilator of f exactly when f * g is the zero function."""
    _cmpr, output_fn = _target()
    _degrees, basis = annihilators(output_fn, verbose=False)
    annihilator = basis[0]
    multiple = AND(output_fn, annihilator).translate_ANF()

    for bits in range(2 ** 3):
        state = [(bits >> i) & 1 for i in range(3)]
        assert output_fn.eval(state) * annihilator.eval(state) == 0, (
            f"g does not annihilate f at {state}"
        )
        assert multiple.eval(state) == 0, f"h = f*g is nonzero at {state}"

def test_the_fast_pair_has_a_nonzero_multiple():
    """If h were zero, FAA's tests would exercise nothing RAA doesn't."""
    _cmpr, output_fn = _target()
    g = VAR(0)
    h = AND(output_fn, g).translate_ANF()

    assert BooleanANF.from_BooleanFunction(h).terms, "h = f*g is the zero function"
    for bits in range(2 ** 3):
        state = [(bits >> i) & 1 for i in range(3)]
        assert h.eval(state) == (output_fn.eval(state) & g.eval(state))


# ── NAA ──────────────────────────────────────────────────────────────────────

@pytest.mark.slow
@_PATHS
@pytest.mark.parametrize("make_solver", [LUSolver, GaussElimSolver],
                         ids=["LUSolver", "GaussElimSolver"])
def test_naa_recovers_the_secret_state(use_profiles, make_solver):
    """NAA layers keystream constants onto an LU-shaped store; both solvers take them."""
    cmpr, output_fn = _target()
    attack_data = NAA_offline(cmpr, output_fn, 0, _TIME_LIMIT, verbose=False,
                              **_offline_path(cmpr, use_profiles))
    keystream = _keystream(cmpr, output_fn, attack_data["keystream needed"])
    recovered = NAA_online(cmpr, output_fn, keystream, attack_data,
                           verbose=False, solver=make_solver())

    assert recovered is not None, "NAA returned no solution"
    assert list(np.asarray(recovered).ravel()) == _secret_state(cmpr)

@pytest.mark.slow
def test_naa_rejects_the_grobner_solver():
    """The one enforced incompatibility: NAA defers keystream constants, which
    Groebner reduction has no way to accept separately."""
    cmpr, output_fn = _target()
    attack_data = NAA_offline(cmpr, output_fn, 0, _TIME_LIMIT, verbose=False,
                              **_offline_path(cmpr, True))
    keystream = _keystream(cmpr, output_fn, attack_data["keystream needed"])

    with pytest.raises(ValueError, match="GrobnerSolver"):
        NAA_online(cmpr, output_fn, keystream, attack_data,
                   verbose=False, solver=GrobnerSolver())


# ── RAA ──────────────────────────────────────────────────────────────────────

@pytest.mark.slow
@_PATHS
@pytest.mark.parametrize("make_kwargs", _ONLINE_COMBINATIONS)
def test_raa_recovers_the_secret_state(use_profiles, make_kwargs):
    cmpr, output_fn = _target()
    _degrees, basis = annihilators(output_fn, verbose=False)
    annihilator = basis[0]
    multiple = AND(output_fn, annihilator).translate_ANF()
    attack_data = RAA_offline(cmpr, annihilator, multiple, 0, 4, _TIME_LIMIT,
                              verbose=False, **_offline_path(cmpr, use_profiles))
    keystream = _keystream(cmpr, output_fn, attack_data["keystream needed"])
    recovered = RAA_online(cmpr, output_fn, keystream, attack_data,
                           verbose=False, **make_kwargs())

    assert recovered is not None, "RAA returned no solution"
    assert list(np.asarray(recovered).ravel()) == _secret_state(cmpr)


# ── FAA ──────────────────────────────────────────────────────────────────────

@pytest.mark.slow
@_PATHS
@pytest.mark.parametrize("make_kwargs", _ONLINE_COMBINATIONS)
def test_faa_recovers_the_secret_state(use_profiles, make_kwargs):
    """FAA in the case it exists for: h nonzero, cancelled by its linear relation."""
    cmpr, output_fn = _target()
    g = VAR(0)
    h = AND(output_fn, g).translate_ANF()
    attack_data = FAA_offline(cmpr, g, h, 0, 4, _TIME_LIMIT,
                              verbose=False, **_offline_path(cmpr, use_profiles))
    keystream = _keystream(cmpr, output_fn, attack_data["keystream needed"])
    recovered = FAA_online(cmpr, output_fn, keystream, attack_data,
                           verbose=False, **make_kwargs())

    assert recovered is not None, "FAA returned no solution"
    assert list(np.asarray(recovered).ravel()) == _secret_state(cmpr)

@pytest.mark.slow
@_PATHS
def test_faa_recovers_with_a_true_annihilator(use_profiles):
    """The degenerate case, h = 0: FAA must still work when handed RAA's pair."""
    cmpr, output_fn = _target()
    _degrees, basis = annihilators(output_fn, verbose=False)
    annihilator = basis[0]
    multiple = AND(output_fn, annihilator).translate_ANF()
    attack_data = FAA_offline(cmpr, annihilator, multiple, 0, 4, _TIME_LIMIT,
                              verbose=False, **_offline_path(cmpr, use_profiles))
    keystream = _keystream(cmpr, output_fn, attack_data["keystream needed"])
    recovered = FAA_online(cmpr, output_fn, keystream, attack_data, verbose=False)

    assert recovered is not None, "FAA returned no solution"
    assert list(np.asarray(recovered).ravel()) == _secret_state(cmpr)
