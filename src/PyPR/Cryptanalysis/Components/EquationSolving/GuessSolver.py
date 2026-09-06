"""Guess-and-prune solver for incomplete equation systems.

After a solver produces a partial solution, remaining unsolved variables
must be guessed. This module provides effect-vector pruning (linear
independence) and exhaustive search with keystream verification.
"""
from PyPR import FeedbackRegister

from itertools import product
import numpy as np
import time


def guess_and_solve(
    feedback_fn,
    output_fn,
    base_solution,
    effect_vectors,
    keystream,
    test_length=1000,
    verbose=False,
    _print_depth=0,
):
    """Prune effect vectors and exhaustively search for a valid initial state.

    Given a base solution and a list of effect vectors (how each guess
    changes the state), this function:

    1. Tests whether the base solution already reproduces the keystream
    2. Prunes linearly dependent effect vectors
    3. Exhaustively tries all 2^k independent guess combinations

    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param base_solution: Initial state vector (unknowns defaulted).
    :type base_solution: np.ndarray[np.uint8]
    :param effect_vectors: List of effect vectors, each showing how a
        single guess bit changes the base solution.
    :type effect_vectors: list[np.ndarray[np.uint8]]
    :param keystream: The observed keystream to verify against.
    :type keystream: np.ndarray[np.uint8]
    :param test_length: Number of keystream bits to use for verification.
    :type test_length: int
    :param verbose: Whether to print progress.
    :type verbose: bool
    :param _print_depth: Indentation level for verbose output.
    :type _print_depth: int
    :return: ``(initial_state, guesses_tried, pruned_guess_bits)`` —
        the recovered state (or None), total guesses tested, and
        number of independent guess dimensions after pruning.
    :rtype: tuple[list[int] | None, int, int]
    """
    start_time = time.time()

    F = FeedbackRegister(0, feedback_fn)
    test_length = min(test_length, len(keystream))
    test_keystream = keystream[:test_length]

    _indent_1 = '|   ' * (_print_depth+1)
    _indent_2 = '|   ' * (_print_depth+2)
    _indent_3 = '|   ' * (_print_depth+3)

    # test if base solution is already correct:
    F.set_state(base_solution)
    test_seq = [output_fn.eval(state) for state in F.run(test_length)]
    if np.all(test_seq == test_keystream):
        if verbose:
            print(f"{_indent_1}Solve complete -- correct base solution")
            print(f"{_indent_1}Time: {time.time() - start_time} s")
        return (list(base_solution), 0, 0)

    if verbose:
        print(_indent_1)
        print(f"{_indent_1}Initial solve failed, guessing remaining information:")
        print(f"{_indent_2}Starting effect pruning:")

    # prune linearly dependent effect vectors:
    effect_pruning_time = time.time()
    pruned_guesses = []
    already_solved = set()
    reduced_matrix = np.zeros([feedback_fn.size, feedback_fn.size], dtype=np.uint8)

    for i, effect_vector in enumerate(effect_vectors):
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

        if verbose:
            print(f"\r{_indent_3}Vectors Pruned: {i+1}/{len(effect_vectors)}", end='')

    if verbose:
        print(f"\n{_indent_2}Pruning finished:")
        print(f"{_indent_2}Max number of guesses (original): 2^{len(effect_vectors)}")
        print(f"{_indent_2}Max number of guesses (pruned): 2^{len(pruned_guesses)}")
        print(f"{_indent_2}Time: {time.time() - effect_pruning_time} s")
        print(_indent_2)
        print(f"{_indent_2}Starting to Guess:")

    # exhaustively test pruned guesses:
    guess_count = 0
    guess_start_time = time.time()
    for guess_assignment in product((0, 1), repeat=len(pruned_guesses)):
        guess_count += 1

        if verbose:
            print(f"\r{_indent_3}Guess count: {guess_count}", end='')

        F.set_state(base_solution)
        for idx, assigned in enumerate(guess_assignment):
            if assigned:
                F._state ^= pruned_guesses[idx]
        F.seed(F._state.copy())

        mismatch = False
        for t, state in enumerate(F.run(test_length)):
            if output_fn.eval(state) != test_keystream[t]:
                mismatch = True
                break

        if not mismatch:
            if verbose:
                print(f"\n{_indent_2}Guessing Finished:")
                print(f"{_indent_2}Time: {time.time() - guess_start_time} s")
                print(f"{_indent_1}Solution Found!")
            F.reset()
            return (list(F), guess_count, len(pruned_guesses))

    if verbose:
        print(f"\n{_indent_2}Guessing exhausted — no solution found.")
    return (None, guess_count, len(pruned_guesses))
