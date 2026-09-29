"""Tests for the branch-and-prune Groebner solver.

The solver guesses one unsolved variable at a time, forking the store into a
var=0 branch and a var=1 branch and processing them in alternating batches.
Whichever branch reaches an algebraic contradiction first confirms the other
one's value, and splitting continues from there.

The end-to-end recovery is covered in `test_attacks.py`, which runs this solver
as RAA's and FAA's online solver. What is covered here is the machinery that
drives it: the fork and how it chooses a variable, its termination condition,
the independence of the two branches, the contradiction that turns a guess into
a confirmed value, and that a real attack actually enters the loop.
"""
import copy
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, CONST, VAR, XOR
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template

from PyPR.FeedbackFunctions import CMPR, MPR
from PyPR.FeedbackRegister import FeedbackRegister

import PyPR.Cryptanalysis.Components.EquationSolving.SplitGrob_Solver as split_grob
from PyPR.Cryptanalysis.Attacks.reduced_algebraic_attack import RAA_offline, RAA_online
from PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator import annihilators
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore

# ── the fork ─────────────────────────────────────────────────────────────────

def test_split_forks_the_store_on_an_unsolved_variable():
    """One branch asserts the variable is 0, the other that it is 1."""
    store = GroebnerEqStore(simplify_mode=None)
    store.enqueue_equation(XOR(VAR(0), VAR(1)))
    store.enqueue_equation(XOR(VAR(1), VAR(2), CONST(1)))
    store.process_pending()

    result = split_grob._split(store, False, "")
    assert result is not None, "a store with unsolved variables must split"
    var, store_0, store_1 = result

    assert var in store.unknown_vars
    # the var=0 branch is the original store, extended in place
    assert store_0 is store
    assert store_1 is not store

    store_0.process_pending()
    store_1.process_pending()
    assert store_0.solved_vars.get(var) == 0
    assert store_1.solved_vars.get(var) == 1

def test_split_returns_none_when_nothing_is_left_to_guess():
    """The loop's termination condition: no unsolved variable means no split."""
    store = GroebnerEqStore(simplify_mode=None)
    store.enqueue_equation(XOR(VAR(0), CONST(1)))
    store.enqueue_equation(VAR(1))
    store.process_pending()

    assert not (store.unknown_vars - set(store.solved_vars.keys()))
    assert split_grob._split(store, False, "") is None

def test_the_two_branches_are_independent():
    """The var=1 branch is a deepcopy; solving in one must not touch the other."""
    store = GroebnerEqStore(simplify_mode=None)
    store.enqueue_equation(XOR(VAR(0), VAR(1)))
    store.enqueue_equation(XOR(VAR(1), VAR(2), CONST(1)))
    store.process_pending()

    result = split_grob._split(store, False, "")
    assert result is not None
    _var, store_0, store_1 = result
    store_0.process_pending()
    before = copy.deepcopy(store_1.solved_vars)
    store_0.enqueue_equation(VAR(2))
    store_0.process_pending()

    assert store_1.solved_vars == before, "solving branch 0 leaked into branch 1"


# ── the step the loop is built from ─────────────────────────────────────────

def test_split_picks_the_highest_indexed_unsolved_variable():
    """The choice is deterministic: `unsolved[-1]`, which the source notes works
    better for CMPRs. Two identical stores therefore split on the same variable.
    """
    choices = set()
    for _ in range(3):
        store = GroebnerEqStore(simplify_mode=None)
        store.enqueue_equation(XOR(VAR(0), VAR(1)))
        store.enqueue_equation(XOR(VAR(1), VAR(2), CONST(1)))
        store.process_pending()
        result = split_grob._split(store, False, "")
        assert result is not None
        choices.add(result[0])

    assert choices == {2}

def test_splitting_a_free_variable_gives_two_consistent_branches():
    """A variable the equations leave free can take either value.

    x0 is pinned to 1, so it is not a split candidate; x1 = x2 leaves the pair
    free, and assuming either value for the split variable is consistent -- so
    neither branch may contradict. A spurious contradiction here would make the
    loop "confirm" a value the system never implied.
    """
    store = GroebnerEqStore(simplify_mode=None)
    store.enqueue_equation(XOR(VAR(0), CONST(1)))      # x0 = 1
    store.enqueue_equation(XOR(VAR(1), VAR(2)))        # x1 = x2
    store.process_pending()
    assert store.solved_vars.get(0) == 1

    result = split_grob._split(store, False, "")
    assert result is not None
    var, store_0, store_1 = result
    assert var != 0, "x0 was already solved, so it should not be the split variable"

    for value, branch in ((0, store_0), (1, store_1)):
        branch.process_pending()                       # raises if inconsistent
        assert branch.solved_vars.get(var) == value

def test_a_branch_against_a_solved_variable_contradicts():
    """Asserting the opposite of a pinned value must be detected, not tolerated --
    it is this contradiction that lets the loop confirm the other branch."""
    store = GroebnerEqStore(simplify_mode=None)
    store.enqueue_equation(XOR(VAR(0), CONST(1)))      # x0 = 1
    store.process_pending()

    store.enqueue_equation(VAR(0))                     # ...and now x0 = 0
    with pytest.raises(ValueError, match="Inconsistent"):
        store.process_pending()


# ── the loop that drives it ──────────────────────────────────────────────────

@pytest.mark.slow
def test_split_loop_runs_during_a_real_attack(monkeypatch):
    """Guards the coverage the end-to-end tests claim.

    If Groebner reduction determined the whole system up front, `_split` would
    return None on its first call and the loop body would never run, so the
    solver would be "covered" without its main loop executing. Whether a split
    is needed depends on the register; chaining is seeded here to one (seed 1)
    that leaves a variable undetermined, and the count makes it visible if a
    change to the solver or the store stops that from being true.
    """
    calls = []
    original = split_grob._split

    def counting_split(grob_store, verbose, indent):
        result = original(grob_store, verbose, indent)
        calls.append(result)
        return result

    monkeypatch.setattr(split_grob, "_split", counting_split)

    random.seed(1)
    np.random.seed(1)
    cmpr = CMPR([MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1]),
                 MPR(3, [1, 1, 0, 1], [1, 0, 1])])
    cmpr.generateChaining(template=arman_template(max_and=2))
    output_fn = XOR(AND(VAR(0), VAR(1)), VAR(2))
    _degrees, basis = annihilators(output_fn, verbose=False)
    annihilator = basis[0]
    multiple = AND(output_fn, annihilator).translate_ANF()

    attack_data = RAA_offline(
        cmpr, annihilator, multiple, 0, 4, 120, verbose=False,
        monomial_profiles=cmpr.monomial_profiles(), variable_blocks=cmpr.blocks,
    )
    keystream = np.array(
        [int(output_fn.eval(state._state))
         for state in FeedbackRegister(42, cmpr).run(attack_data["keystream needed"],
                                                     compiled=False)],
        dtype=np.uint8,
    )
    recovered = RAA_online(
        cmpr, output_fn, keystream, attack_data, verbose=False,
        solver=split_grob.SplitGrobnerSolver(),
        online_store=GroebnerEqStore(simplify_mode=None),
    )

    splits = [call for call in calls if call is not None]
    assert splits, "the split loop never ran, so the solver's main path went untested"
    assert calls[-1] is None, "the loop should end by running out of variables to guess"
    assert list(np.asarray(recovered).ravel()) == list(FeedbackRegister(42, cmpr)._state)
