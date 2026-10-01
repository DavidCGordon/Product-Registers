"""Guess-and-prune solver for incomplete equation systems.

After a solver produces a partial solution, remaining unsolved variables
must be guessed. This module provides effect-vector pruning (linear
independence) and exhaustive search with keystream verification.
"""
import logging
from collections.abc import Callable
from itertools import product
from typing import Any, NamedTuple

import numpy as np

from PyPR import FeedbackRegister
from PyPR.Reporting import get_logger

log = get_logger(__name__)


class SolveResult(NamedTuple):
    """What every solver returns.

    A named tuple, so existing `state, guesses, bits = solver.solve(...)`
    unpacking keeps working alongside the named fields.

    :ivar state: The recovered initial state as a list of bits (numpy uint8
        values, as the solvers' arrays hold them), or None if no candidate
        passed verification.
    :ivar guesses: How many candidates were tried after the base solution.
    :ivar guess_bits: The number of independent guess dimensions left after
        pruning, so the search space was `2 ** guess_bits`.
    """
    state: list[Any] | None
    guesses: int
    guess_bits: int


def keystream_verifier(
    feedback_fn,
    output_fn,
    keystream,
    test_length: int = 1000,
) -> Callable[[np.ndarray], bool]:
    """Build the default check: does a candidate state reproduce the keystream?

    The candidate is loaded into a register over `feedback_fn` and its first
    `test_length` output bits (through `output_fn`) are compared with the
    observed keystream.

    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param keystream: The observed keystream.
    :type keystream: np.ndarray[np.uint8]
    :param test_length: Number of keystream bits to compare. Defaults to 1000.
    :type test_length: int
    :raises ValueError: If `keystream` is None.
    :return: A function accepting a candidate state and returning whether it
        reproduces the keystream.
    :rtype: Callable[[np.ndarray[np.uint8]], bool]
    """
    if keystream is None:
        raise ValueError("guess_and_solve needs a keystream or a verify function")

    F = FeedbackRegister(0, feedback_fn)
    # Use the compiled clocking when the caller has compiled the function,
    # and the uncompiled path otherwise -- they produce the same states.
    # Defaulting to compiled made this raise for any caller that hadn't,
    # which is every attack's dynamic (no monomial profile) path except
    # FAA's, whose offline phase happens to compile the function for its
    # own use.
    compiled = getattr(feedback_fn, "_compiled", None) is not None
    test_length = min(test_length, len(keystream))
    test_keystream = keystream[:test_length]

    def keystream_matches(candidate):
        F.set_state(candidate)
        for t, state in enumerate(F.run(test_length, compiled=compiled)):
            if output_fn.eval(state) != test_keystream[t]:
                return False
        return True

    return keystream_matches


@log.stage("Guessing remaining state")
def guess_and_solve(
    feedback_fn,
    output_fn,
    base_solution,
    effect_vectors,
    keystream,
    test_length=1000,
    verify=None,
) -> SolveResult:
    """Prune effect vectors and exhaustively search for a valid initial state.

    Given a base solution and a list of effect vectors (how each guess
    changes the state), this function:

    1. Tests whether the base solution already passes verification
    2. Prunes linearly dependent effect vectors
    3. Exhaustively tries all 2^k independent guess combinations

    A candidate passes when its keystream matches the observed one, or, if
    `verify` is given, when `verify` accepts it. The hook lets an attack check
    candidates against a target it reaches only through an interface -- a
    chosen-IV oracle, say, whose initialization rounds and keystream stay
    inside it -- instead of against a keystream handed to the solver.

    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param base_solution: Initial state vector (unknowns defaulted).
    :type base_solution: np.ndarray[np.uint8]
    :param effect_vectors: List of effect vectors, each showing how a
        single guess bit changes the base solution.
    :type effect_vectors: list[np.ndarray[np.uint8]]
    :param keystream: The observed keystream to verify against. Unused, and
        may be None, when `verify` is given.
    :type keystream: np.ndarray[np.uint8] | None
    :param test_length: Number of keystream bits to use for verification.
    :type test_length: int
    :param verify: Decides whether a candidate initial state is correct, in
        place of comparing its keystream with `keystream` (which may then be
        None).
    :type verify: Callable[[np.ndarray[np.uint8]], bool] | None
    :return: The recovered state (or None), the number of guesses tried, and
        the number of independent guess dimensions after pruning.
    :rtype: SolveResult
    :raises ValueError: If neither a keystream nor `verify` is given.
    """
    if verify is None:
        verify = keystream_verifier(feedback_fn, output_fn, keystream, test_length)

    # test if base solution is already correct:
    if verify(base_solution):
        log.info("Base solution is already correct")
        return SolveResult(list(base_solution), 0, 0)

    log.debug("Base solution failed verification")

    # prune linearly dependent effect vectors:
    log.step("Pruning effect vectors", level=logging.DEBUG)
    pruning = log.progress("Vectors pruned", total=len(effect_vectors))
    pruned_guesses = []
    already_solved = set()
    reduced_matrix = np.zeros([feedback_fn.size, feedback_fn.size], dtype=np.uint8)

    for effect_vector in effect_vectors:
        effect_vector_copy = effect_vector.copy()
        for idx in range(len(effect_vector)):
            if effect_vector[idx] == 1:
                if idx in already_solved:
                    effect_vector ^= reduced_matrix[idx]
                else:
                    already_solved.add(idx)
                    pruned_guesses.append(effect_vector_copy)
                    reduced_matrix[idx] = effect_vector
                    break
        pruning.update()
    pruning.close()
    log.info("Guess space: 2^%d, pruned to 2^%d", len(effect_vectors), len(pruned_guesses))

    # exhaustively test pruned guesses:
    log.step("Guessing")
    guessing = log.progress("Guesses", total=2 ** len(pruned_guesses), unit="guesses")
    for guess_assignment in product((0, 1), repeat=len(pruned_guesses)):
        guessing.update()

        candidate = base_solution.copy()
        for idx, assigned in enumerate(guess_assignment):
            if assigned:
                candidate ^= pruned_guesses[idx]

        if verify(candidate):
            guessing.close()
            log.info("Solution found after %d guesses", guessing.count)
            return SolveResult(list(candidate), guessing.count, len(pruned_guesses))

    guessing.close()
    log.warning("Guessing exhausted -- no solution found")
    return SolveResult(None, guessing.count, len(pruned_guesses))
