"""Batch Gröbner basis solver.

Accepts equations from any store (via ``to_anf_list``), feeds them
into a fresh :class:`GroebnerEqStore`, and returns the solved
variable assignments from Buchberger reduction + unit propagation.
"""
import numpy as np

from PyPR.Reporting import get_logger

from PyPR.Cryptanalysis.Components.Adapters.store_repr import to_anf_list
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore

log = get_logger(__name__)


def reduce(equation_store, simplify_mode=None):
    """Reduce a system of GF(2) polynomial equations via Gröbner basis computation.

    Produces a GroebnerEqStore with as many variables determined as
    Buchberger reduction + unit propagation can reach. Undetermined
    variables remain in ``unknown_vars``.

    Accepts any equation store. Equations are extracted as BooleanANF
    (via ``to_anf_list``), then fed into a fresh GroebnerEqStore which
    performs Buchberger reduction with unit propagation.

    :param equation_store: Any equation store. GroebnerEqStore is used
        natively; others are converted via ``to_anf_list``.
    :param simplify_mode: Simplification strategy for the underlying
        GroebnerEqStore. Defaults to None.
    :return: A GroebnerEqStore containing the reduced Gröbner basis.
        Access ``.solved_vars`` for determined variable values.
    :rtype: GroebnerEqStore
    """
    if isinstance(equation_store, GroebnerEqStore):
        return equation_store

    anf_list = to_anf_list(equation_store)
    gb = GroebnerEqStore(simplify_mode=simplify_mode)
    for eq in anf_list:
        bf = eq.to_BooleanFunction()
        gb.enqueue_equation(bf)

    gb.process_pending()
    return gb


@log.stage("Groebner solve")
def solve(
    equation_store,
    feedback_fn, output_fn, keystream,
    test_length=1000, verify=None, simplify_mode=None,
):
    """Solve a system of GF(2) equations via Gröbner reduction + exhaustive guess.

    Runs the Gröbner reducer to extract determined variable assignments,
    builds a base solution from ``solved_vars``, and delegates to
    :func:`GuessSolver.guess_and_solve` for exhaustive search over
    the remaining free variables. Guarantees a solution if one
    exists and the guess space is feasible.

    In a fully reduced Gröbner basis, unsolved variables are genuinely
    free — they have no polynomial relationships left in the basis.
    Each free state bit contributes a unit-vector effect, so the
    guess space is exactly 2^k where k is the number of undetermined
    state bits.

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
        guess_and_solve,
    )

    log.step("Groebner reduction")
    grob_store = reduce(equation_store, simplify_mode=simplify_mode)

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

    log.info("Variables solved: %d/%d, free: %d", n - len(effect_vectors), n, len(effect_vectors))

    log.step("Guessing")
    return guess_and_solve(
        feedback_fn, output_fn, base_solution, effect_vectors,
        keystream, test_length=test_length, verify=verify,
    )


class GrobnerSolver:
    """Object-oriented wrapper for Gröbner basis solving.

    :param simplify_mode: Simplification strategy for the underlying
        GroebnerEqStore. Defaults to ``None``.
    :type simplify_mode: str | None
    """

    def __init__(self, *, simplify_mode=None):
        self.simplify_mode = simplify_mode

    def reduce(self, equation_store):
        return reduce(equation_store, simplify_mode=self.simplify_mode)

    def solve(
        self, equation_store, feedback_fn, output_fn, keystream, *,
        test_length=1000, verify=None,
    ):
        return solve(
            equation_store, feedback_fn, output_fn, keystream,
            test_length=test_length, verify=verify, simplify_mode=self.simplify_mode,
        )
