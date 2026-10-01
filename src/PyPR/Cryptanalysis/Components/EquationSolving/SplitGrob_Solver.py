"""Split Groebner solver -- branch-and-prune over GF(2) polynomials.

Picks unsolved variables one at a time, guessing both values (0 and 1).
Both branches are processed in alternating batches. The wrong guess
should hit an algebraic inconsistency quickly (the correct branch tends
to churn through many non-informative syzygies). When a branch detects
a contradiction, the variable's value is confirmed and the surviving
branch's store is carried forward to split on the next variable.

Occasionally neither guess ever contradicts -- being algebraically
consistent doesn't prove a guess matches the real keystream, just that
it doesn't violate the (possibly incomplete) equations gathered so far.
So whenever a branch fully determines the system, it's checked directly
against the real keystream. If it matches, that's the answer. If not,
the guess was wrong despite being consistent -- the other,
not-yet-finished branch is adopted instead and splitting continues from
there.

At INFO the loop is one progress bar: variables confirmed, out of those
unsolved at the start. At DEBUG each split also shows its two branches as
live lines (processed-equation count, and basis, queue and solved sizes),
replaced by one line giving the verdict when the split resolves.
"""
import copy
import logging

import numpy as np

from PyPR.Reporting import get_logger

from PyPR.BooleanLogic.FunctionInputs import CONST, VAR
from PyPR.BooleanLogic.Gates import XOR

from PyPR.Cryptanalysis.Components.Adapters.store_repr import to_anf_list
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore

log = get_logger(__name__)

_BRANCH_BATCH_SIZE = 50


def _split(grob_store):
    """Pick the next unsolved variable in `grob_store` and fork it: the
    store becomes the var=0 branch in place (no copy), and a single
    deepcopy becomes the var=1 branch.

    :return: ``(var, store_0, store_1)``, or ``None`` if every variable
        in `grob_store` is already solved.
    :rtype: tuple[int, GroebnerEqStore, GroebnerEqStore] | None
    """
    unsolved = sorted(grob_store.unknown_vars - set(grob_store.solved_vars.keys()))
    if not unsolved:
        return None
    var = unsolved[-1] # works better for CMPRs
    store_1 = copy.deepcopy(grob_store)
    grob_store.enqueue_equation(VAR(var))
    store_1.enqueue_equation(XOR(VAR(var), CONST(1)))
    return var, grob_store, store_1


def _branch_status(store):
    return f"Basis: {store.num_eqs} -- Queue: {len(store.queue)} -- Solved: {len(store.solved_vars)}"


@log.stage("Split-Groebner solve")
def solve(
    equation_store,
    feedback_fn, output_fn, keystream,
    test_length=1000, verify=None, simplify_mode=None,
):
    """Solve a GF(2) system via batch branch-and-prune over Groebner bases.

    Loads equations into a GroebnerEqStore without initial reduction,
    then picks unsolved variables one at a time, guessing both values
    (0 and 1) and processing alternating batches. The wrong guess
    should hit an algebraic inconsistency quickly; when it does,
    the variable's value is confirmed and the surviving branch's store
    is carried forward to split on the next unsolved variable. This
    repeats until a branch fully determines the system.

    Full determination is checked directly against the real keystream
    before being trusted -- it only proves consistency with the equations
    gathered so far, not that the guess is actually correct. A mismatch
    means the guess was wrong despite being consistent; the other,
    not-yet-finished branch is adopted and splitting continues.

    Any variables left unresolved (not referenced by any equation) are
    delegated to :func:`GuessSolver.guess_and_solve` for keystream
    verification.

    :param equation_store: Any equation store (converted to
        GroebnerEqStore if needed).
    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param keystream: The observed keystream to verify against.
    :type keystream: np.ndarray[np.uint8]
    :param test_length: Number of keystream bits to use for verification.
    :type test_length: int
    :param verify: Decides whether a candidate initial state is correct, in
        place of comparing its keystream with `keystream` (which may then be
        None). Forwarded to :func:`GuessSolver.guess_and_solve`.
    :type verify: Callable[[np.ndarray[np.uint8]], bool] | None
    :param simplify_mode: Simplification strategy for the GroebnerEqStore.
    :type simplify_mode: str | None
    :return: The recovered state (or None), the number of guesses tried, and
        the number of independent guess dimensions after pruning.
    :rtype: SolveResult
    """
    from PyPR.Cryptanalysis.Components.EquationSolving.GuessSolver import (
        SolveResult,
        guess_and_solve,
        keystream_verifier,
    )

    # --- Load equations into GroebnerEqStore (no initial reduction) ---
    if isinstance(equation_store, GroebnerEqStore):
        grob_store = equation_store
    else:
        grob_store = GroebnerEqStore(simplify_mode=simplify_mode)
        for eq in to_anf_list(equation_store):
            grob_store.enqueue_equation(eq.to_BooleanFunction())

    n = feedback_fn.size
    unsolved = sorted(grob_store.unknown_vars - set(grob_store.solved_vars.keys()))
    if unsolved:
        log.info("%d variables unsolved", len(unsolved))
    else:
        log.info("All variables already solved -- no splits needed")

    # a branch that determines the system is checked against the keystream
    # before it is trusted; this is the same check guess_and_solve applies
    matches = verify if verify is not None else keystream_verifier(
        feedback_fn, output_fn, keystream, test_length
    )

    # --- Split-and-prune loop ---
    log.step("Splitting")
    confirmations = log.progress("Variables confirmed", total=len(unsolved))
    guesses_made = 0
    confirmed = 0

    split = _split(grob_store)
    if split is not None:
        guesses_made += 1
    processed_0 = processed_1 = 0
    branches = None

    # `split` is the loop's state: a guessed variable and the two branch stores
    # that assume it 0 and 1. It is None exactly when there is nothing left to
    # guess, so testing it is both the termination check and the guarantee that
    # the three names below are bound.
    while split is not None:
        var, store_0, store_1 = split
        if branches is None:
            branches = (
                log.progress(f"x_{var} = 0) processed", level=logging.DEBUG),
                log.progress(f"x_{var} = 1) processed", level=logging.DEBUG),
            )

        contradicted = 0
        try:
            processed_0 += store_0.process_pending(batch_size=_BRANCH_BATCH_SIZE)
            contradicted = 1
            processed_1 += store_1.process_pending(batch_size=_BRANCH_BATCH_SIZE)
        except ValueError:
            confirmed += 1
            confirmed_value = 1 - contradicted
            surviving = store_1 if contradicted == 0 else store_0
            surviving.solved_vars[var] = confirmed_value
            surviving.unknown_vars.discard(var)
            surviving._simplify()
            for branch in branches:
                branch.close(quiet=True)
            branches = None
            confirmations.update()
            log.debug(
                "x_%d = %d) found inconsistent => x_%d = %d) confirmed (solved: %d/%d)",
                var, contradicted, var, confirmed_value, len(surviving.solved_vars), n,
            )
            grob_store = surviving
            split = _split(grob_store)
            if split is not None:
                guesses_made += 1
            processed_0 = processed_1 = 0
            continue

        # A branch is only trustworthy once its queue is *also* empty --
        # is_determined alone just means every variable seen so far has a
        # value; a still-pending syzygy could yet reduce to a contradiction.
        if store_0.is_determined:
            winner, other, winner_idx = store_0, store_1, 0
        elif store_1.is_determined:
            winner, other, winner_idx = store_1, store_0, 1
        else:
            branches[0].update_to(processed_0)
            branches[0].set_status(_branch_status(store_0))
            branches[1].update_to(processed_1)
            branches[1].set_status(_branch_status(store_1))
            continue

        for branch in branches:
            branch.close(quiet=True)
        branches = None

        # Being algebraically consistent only proves `winner` doesn't
        # violate the equations gathered so far -- it doesn't prove the
        # guess is actually correct. Check directly against the real
        # keystream before trusting it.
        candidate = np.zeros(n, dtype=np.uint8)
        for i, val in winner.solved_vars.items():
            candidate[i] = val
        if matches(candidate):
            confirmations.close()
            log.info("x_%d = %d) matches the keystream -- done", var, winner_idx)
            return SolveResult(list(candidate), 0, 0)

        # Wrong guess despite consistency -- the other, not-yet-finished
        # branch must be the correct one; adopt it and keep splitting.
        other_idx = 1 - winner_idx
        other.solved_vars[var] = other_idx
        other.unknown_vars.discard(var)
        other._simplify()
        confirmed += 1
        confirmations.update()
        log.debug(
            "x_%d = %d) consistent but wrong keystream => x_%d = %d) confirmed (solved: %d/%d)",
            var, winner_idx, var, other_idx, len(other.solved_vars), n,
        )
        grob_store = other
        split = _split(grob_store)
        if split is not None:
            guesses_made += 1
        processed_0 = processed_1 = 0

    confirmations.close()

    # --- Build base solution and effect vectors ---
    base_solution = np.zeros(n, dtype=np.uint8)
    effect_vectors = []

    for i in range(n):
        if i in grob_store.solved_vars:
            base_solution[i] = grob_store.solved_vars[i]
        else:
            effect = np.zeros(n, dtype=np.uint8)
            effect[i] = 1
            effect_vectors.append(effect)

    log.info(
        "Variables solved: %d/%d, free: %d, guesses made: %d (%d confirmed)",
        n - len(effect_vectors), n, len(effect_vectors), guesses_made, confirmed,
    )

    log.step("Guessing")
    return guess_and_solve(
        feedback_fn, output_fn, base_solution, effect_vectors,
        keystream, test_length=test_length, verify=verify,
    )


class SplitGrobnerSolver:
    """Batch branch-and-prune Groebner solver.

    After initial Groebner reduction, guesses unsolved variables by
    forking the store into two branches (var=0, var=1) and processing
    them in alternating batches. When one branch detects an algebraic
    inconsistency, the variable's value is confirmed and the surviving
    branch is carried forward.

    :param simplify_mode: Simplification strategy for the underlying
        GroebnerEqStore. Defaults to ``None``.
    :type simplify_mode: str | None
    """

    def __init__(self, *, simplify_mode=None):
        self.simplify_mode = simplify_mode

    def solve(
        self, equation_store, feedback_fn, output_fn, keystream, *,
        test_length=1000, verify=None,
    ):
        return solve(
            equation_store, feedback_fn, output_fn, keystream,
            test_length=test_length, verify=verify, simplify_mode=self.simplify_mode,
        )
