"""Batch Gröbner basis solver with Gebauer-Möller optimizations.

Drop-in replacement for :mod:`Grob_Solver` backed by
:class:`GroebnerEqStore2`, which uses deferred S-polynomial
computation, Gebauer-Möller B/F/M pair pruning, interreduction,
and non-constant linear-lead propagation.

Follows the same function + thin-wrapper-class convention as
the other solvers (see ``EquationSolving/__init__.py``).
"""
import numpy as np
import time

from PyPR.Cryptanalysis.Components.Adapters.store_repr import to_anf_list
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore2 import (
    GroebnerEqStore2,
)


def reduce(equation_store, simplify_mode=None, linear_sub_threshold=4,
           verbose=False, _print_depth=0):
    """Reduce a GF(2) polynomial system via Gröbner basis computation.

    Produces a :class:`GroebnerEqStore2` with as many variables
    determined as Buchberger reduction (with GM pair management)
    + unit propagation + linear-lead propagation can reach.

    Accepts any equation store.  A ``GroebnerEqStore2`` is returned
    as-is; other stores are converted via ``to_anf_list``.

    :param equation_store: Any equation store.
    :param simplify_mode: Simplification strategy for the
        GroebnerEqStore2.
    :type simplify_mode: str | None
    :param linear_sub_threshold: Maximum monomial count in a
        linear polynomial's tail for non-constant propagation.
    :type linear_sub_threshold: int
    :param verbose: Print progress.
    :type verbose: bool
    :param _print_depth: Indentation level for verbose output.
    :type _print_depth: int
    :return: A GroebnerEqStore2 containing the reduced basis.
    :rtype: GroebnerEqStore2
    """
    if isinstance(equation_store, GroebnerEqStore2):
        return equation_store

    anf_list = to_anf_list(equation_store)
    gb = GroebnerEqStore2(
        simplify_mode=simplify_mode,
        linear_sub_threshold=linear_sub_threshold,
    )
    for eq in anf_list:
        gb.enqueue_equation(eq.to_BooleanFunction())

    gb.process_pending(verbose=verbose, _print_depth=_print_depth)
    return gb


def solve(
    equation_store,
    feedback_fn, output_fn, keystream,
    test_length=1000, simplify_mode=None, linear_sub_threshold=4,
    verbose=False, _print_depth=0,
):
    """Solve a GF(2) system via Gröbner (GM) reduction + exhaustive guess.

    Runs the GM-enhanced Gröbner reducer, builds a base solution from
    ``solved_vars``, and delegates to
    :func:`GuessSolver.guess_and_solve` for exhaustive search over
    the remaining free variables.

    :param equation_store: Any equation store (converted to
        GroebnerEqStore2 if needed).
    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param keystream: The observed keystream to verify against.
    :type keystream: np.ndarray[np.uint8]
    :param test_length: Number of keystream bits for verification.
    :type test_length: int
    :param simplify_mode: Simplification strategy for GroebnerEqStore2.
    :type simplify_mode: str | None
    :param linear_sub_threshold: Maximum monomial count in a linear
        polynomial's tail for non-constant propagation.
    :type linear_sub_threshold: int
    :param verbose: Print progress.
    :type verbose: bool
    :param _print_depth: Indentation level for verbose output.
    :type _print_depth: int
    :return: ``(initial_state, guesses_tried, pruned_guess_bits)``
    :rtype: tuple[list[int] | None, int, int]
    """
    from PyPR.Cryptanalysis.Components.EquationSolving.GuessSolver import (
        guess_and_solve,
    )

    reduction_start = time.time()
    grob_store = reduce(
        equation_store,
        simplify_mode=simplify_mode,
        linear_sub_threshold=linear_sub_threshold,
        verbose=verbose,
        _print_depth=_print_depth,
    )
    reduction_time = time.time() - reduction_start

    n = feedback_fn.size
    base_solution = np.zeros(n, dtype=np.uint8)
    effect_vectors = []

    for i in range(n):
        if i in grob_store.solved_vars:
            base_solution[i] = grob_store.solved_vars[i]
        else:
            effect = np.zeros(n, dtype=np.uint8)
            effect[i] = 1
            effect_vectors.append(effect)

    if verbose:
        print(f"{'|   ' * (_print_depth+1)}Groebner (GM) solve complete:")
        print(f"{'|   ' * (_print_depth+2)}Variables solved:"
              f" {n - len(effect_vectors)}/{n}")
        print(f"{'|   ' * (_print_depth+2)}Free variables:"
              f" {len(effect_vectors)}")
        print(f"{'|   ' * (_print_depth+2)}Time: {reduction_time} s")

    return guess_and_solve(
        feedback_fn, output_fn, base_solution, effect_vectors,
        keystream, test_length=test_length,
        verbose=verbose, _print_depth=_print_depth,
    )


class GrobnerSolver2:
    """Gröbner solver with Gebauer-Möller pair management.

    :param simplify_mode: Simplification strategy for the underlying
        GroebnerEqStore2.  Defaults to ``None``.
    :type simplify_mode: str | None
    :param linear_sub_threshold: Maximum monomial count in a linear
        polynomial's tail before non-constant propagation is skipped.
    :type linear_sub_threshold: int
    """

    def __init__(self, *, simplify_mode=None, linear_sub_threshold=4):
        self.simplify_mode = simplify_mode
        self.linear_sub_threshold = linear_sub_threshold

    def reduce(self, equation_store, verbose=False, _print_depth=0):
        return reduce(
            equation_store,
            simplify_mode=self.simplify_mode,
            linear_sub_threshold=self.linear_sub_threshold,
            verbose=verbose,
            _print_depth=_print_depth,
        )

    def solve(
        self, equation_store, feedback_fn, output_fn, keystream, *,
        test_length=1000, verbose=False, _print_depth=0,
    ):
        return solve(
            equation_store, feedback_fn, output_fn, keystream,
            test_length=test_length,
            simplify_mode=self.simplify_mode,
            linear_sub_threshold=self.linear_sub_threshold,
            verbose=verbose,
            _print_depth=_print_depth,
        )
