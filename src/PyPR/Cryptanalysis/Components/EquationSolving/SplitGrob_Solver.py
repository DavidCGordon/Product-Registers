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
against the real keystream via :func:`GuessSolver.guess_and_solve`. If it
matches, that's the answer. If not, the guess was wrong despite being
consistent -- the other, not-yet-finished branch is adopted instead and
splitting continues from there.
"""
import copy
import time
import numpy as np

from PyPR.Cryptanalysis.Components.Adapters.store_repr import to_anf_list
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore
from PyPR.BooleanLogic.FunctionInputs import VAR, CONST
from PyPR.BooleanLogic.Gates import XOR


_BRANCH_BATCH_SIZE = 50


def _split(grob_store, verbose, indent):
    """Pick the next unsolved variable in `grob_store` and fork it: the
    store becomes the var=0 branch in place (no copy), and a single
    deepcopy becomes the var=1 branch.

    Prints (and \\033[s-saves the position of) the "Splitting on x_i..."
    header -- everything printed by :func:`_print_footer`/:func:`_replace_header`
    afterward is anchored to that single saved position, so this must be
    the first of the three to run for a given split.

    :return: ``(var, store_0, store_1)``, or ``None`` if every variable
        in `grob_store` is already solved.
    :rtype: tuple[int, GroebnerEqStore, GroebnerEqStore] | None
    """
    unsolved = sorted(grob_store.unknown_vars - set(grob_store.solved_vars.keys()))
    if not unsolved:
        return None
    var = unsolved[-1] # works better for CMPRs
    if verbose:
        print("\033[s", end='')
        print(f"{indent}Splitting on x_{var}...")
    store_1 = copy.deepcopy(grob_store)
    grob_store.enqueue_equation(VAR(var))
    store_1.enqueue_equation(XOR(VAR(var), CONST(1)))
    return var, grob_store, store_1


def _print_footer(line0, line1):
    """Redraw the live 2-line batch-progress footer in place, directly
    below the header `_split` printed (and saved the position of).
    """
    # \033[u jumps back to exactly where \033[s saved it (the start of
    # the header line), regardless of how many rows anything since has
    # wrapped into. \033[1B\r then steps down past that one header line.
    print("\033[u\033[1B\r", end='')
    print(f"\033[K{line0}")
    print(f"\033[K{line1}", end='', flush=True)


def _replace_header(indent, text):
    """Replace the header `_split` printed -- and clear any footer drawn
    below it -- with a single permanent summary line.
    """
    # \033[J erases from the restored position to the end of the screen,
    # so this works whether or not a footer was ever drawn below the
    # header, and regardless of how many rows it wrapped into.
    print("\033[u\033[J", end='')
    print(f"{indent}{text}")


def solve(
    equation_store,
    feedback_fn, output_fn, keystream,
    test_length=1000, simplify_mode=None,
    verbose=False, _print_depth=0,
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
    :param simplify_mode: Simplification strategy for the GroebnerEqStore.
    :type simplify_mode: str | None
    :param verbose: Whether to print progress.
    :type verbose: bool
    :param _print_depth: Indentation level for verbose output.
    :type _print_depth: int
    :return: ``(initial_state, guesses_tried, pruned_guess_bits)`` --
        the recovered state (or None), total guesses tested, and
        independent guess dimensions after pruning.
    :rtype: tuple[list[int] | None, int, int]
    """
    from PyPR.Cryptanalysis.Components.EquationSolving.GuessSolver import guess_and_solve

    _indent_1 = '|   ' * (_print_depth + 1)
    _indent_2 = '|   ' * (_print_depth + 2)
    _indent_3 = '|   ' * (_print_depth + 3)

    # --- Load equations into GroebnerEqStore (no initial reduction) ---
    if isinstance(equation_store, GroebnerEqStore):
        grob_store = equation_store
    else:
        grob_store = GroebnerEqStore(simplify_mode=simplify_mode)
        for eq in to_anf_list(equation_store):
            grob_store.enqueue_equation(eq.to_BooleanFunction())

    n = feedback_fn.size
    unsolved = sorted(grob_store.unknown_vars - set(grob_store.solved_vars.keys()))

    if verbose:
        if unsolved:
            print(f"{_indent_1}Starting split-reduce ({len(unsolved)} unsolved):")
        else:
            print(f"{_indent_1}All variables already solved -- no splits needed.")

    # --- Split-and-prune loop ---
    split_start = time.time()
    guesses_made = 0
    confirmed = 0

    split = _split(grob_store, verbose, _indent_2)
    if split is None:
        finished = True
    else:
        var, store_0, store_1 = split
        guesses_made += 1
        finished = False
    processed_0 = processed_1 = 0

    while not finished:
        try:
            contradicted = 0
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
            if verbose:
                _replace_header(
                    _indent_2,
                    f"x_{var} = {contradicted}) found inconsistent"
                    f" => x_{var} = {confirmed_value}) confirmed"
                    f" (solved: {len(surviving.solved_vars)}/{n})",
                )
            grob_store = surviving
            split = _split(grob_store, verbose, _indent_2)
            if split is None:
                finished = True
            else:
                var, store_0, store_1 = split
                guesses_made += 1
                finished = False
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
            if verbose:
                _print_footer(
                    f"{_indent_3}x_{var} = 0) Processed: {processed_0}"
                    f"  --  Basis: {store_0.num_eqs}"
                    f"  --  Queue: {len(store_0.queue)}"
                    f"  --  Solved: {len(store_0.solved_vars)}",
                    f"{_indent_3}x_{var} = 1) Processed: {processed_1}"
                    f"  --  Basis: {store_1.num_eqs}"
                    f"  --  Queue: {len(store_1.queue)}"
                    f"  --  Solved: {len(store_1.solved_vars)}",
                )
            continue

        # Being algebraically consistent only proves `winner` doesn't
        # violate the equations gathered so far -- it doesn't prove the
        # guess is actually correct. Check directly against the real
        # keystream before trusting it.
        candidate = np.zeros(n, dtype=np.uint8)
        for i, val in winner.solved_vars.items():
            candidate[i] = val
        result = guess_and_solve(
            feedback_fn, output_fn, candidate, [], keystream,
            test_length=test_length, verbose=False, _print_depth=_print_depth,
        )
        if result[0] is not None:
            if verbose:
                _replace_header(_indent_2, f"x_{var} = {winner_idx}) matches the keystream -- done")
            return result

        # Wrong guess despite consistency -- the other, not-yet-finished
        # branch must be the correct one; adopt it and keep splitting.
        other_idx = 1 - winner_idx
        other.solved_vars[var] = other_idx
        other.unknown_vars.discard(var)
        other._simplify()
        confirmed += 1
        if verbose:
            _replace_header(
                _indent_2,
                f"x_{var} = {winner_idx}) consistent but wrong keystream"
                f" => x_{var} = {other_idx}) confirmed"
                f" (solved: {len(other.solved_vars)}/{n})",
            )
        grob_store = other
        split = _split(grob_store, verbose, _indent_2)
        if split is None:
            finished = True
        else:
            var, store_0, store_1 = split
            guesses_made += 1
            finished = False
        processed_0 = processed_1 = 0

    split_time = time.time() - split_start

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

    if verbose:
        print(f"{_indent_1}Split-Groebner solve complete:")
        print(f"{_indent_2}Variables solved: {n - len(effect_vectors)}/{n}")
        print(f"{_indent_2}Free variables: {len(effect_vectors)}")
        print(f"{_indent_2}Guesses made: {guesses_made} ({confirmed} confirmed)")
        print(f"{_indent_2}Time: {split_time:.3f} s")

    return guess_and_solve(
        feedback_fn, output_fn, base_solution, effect_vectors,
        keystream, test_length=test_length,
        verbose=verbose, _print_depth=_print_depth,
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
        test_length=1000, verbose=False, _print_depth=0,
    ):
        return solve(
            equation_store, feedback_fn, output_fn, keystream,
            test_length=test_length, simplify_mode=self.simplify_mode,
            verbose=verbose, _print_depth=_print_depth,
        )
