import numpy as np
import numba
import time

from PyPR.Cryptanalysis.Components.Adapters.store_repr import to_coef_matrix

u8 = numba.types.uint8
@numba.njit(u8[:](u8[:,:],u8[:,:],u8[:]))
def lu_solve(L,U,z):

    # backsolve L: z  =>  L^{-1}z
    for i in range(len(z)-1):
        for j in range(i+1,len(z)):
            z[j] ^= L[j,i] * z[i]

    # backsolve U: L^{-1}z => U^{-1}L^{-1}z
    # equivalently: (LU)^{-1}(z + c)
    for i in range(len(z)-1,0,-1):
        for j in range(i):
            z[j] ^= U[j,i] * z[i]

    return z


def _as_lu_store(equation_store):
    """Return an LUEqStore, converting from other store types if needed."""
    if hasattr(equation_store, 'upper_matrix') and hasattr(equation_store, 'lower_matrix'):
        return equation_store

    from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

    matrix, c2i, _ = to_coef_matrix(equation_store)
    lu = LUEqStore(comb_to_idx=c2i)
    for row in range(matrix.shape[0]):
        lu.insert_equation(matrix[row])
    return lu


def reduce(equation_store, additional_constants=None):
    """Reduce a system of GF(2) equations via LU back-substitution.

    Produces a solution vector from the LU decomposition alone, without
    guessing free variables. If the system is underdetermined, unsolved
    variables default to zero.

    Accepts any equation store. Non-LU stores are converted via
    ``to_coef_matrix`` → LUEqStore (expensive but valid).

    :param equation_store: Any equation store (LUEqStore is native;
        others are converted automatically).
    :param additional_constants: Optional constant vector of length
        ``num_vars``. Defaults to zero vector.
    :type additional_constants: np.ndarray, optional
    :return: Solution vector of length ``num_vars``.
    :rtype: np.ndarray[np.uint8]
    """
    equation_store = _as_lu_store(equation_store)

    if (
        additional_constants is not None and
        len(additional_constants) != (equation_store.num_vars)
    ):
        raise ValueError(
            f"additional_constants has length {len(additional_constants)}, "
            f"but equation store has {equation_store.num_vars} variables."
        )

    if additional_constants is None:
        constant_vector = np.zeros([equation_store.num_vars], dtype=np.uint8)
    else:
        constant_vector = additional_constants.copy()

    if equation_store.consistent and tuple() in equation_store.comb_to_idx:
        constant_vector[equation_store.comb_to_idx[tuple()]] = 1

    n = equation_store.num_vars
    solution = lu_solve(
        equation_store.lower_matrix[:n,:n],
        equation_store.upper_matrix[:n,:n],
        constant_vector,
    )

    return solution


def solve(
    equation_store, feedback_fn, output_fn, keystream, *,
    test_length=1000, additional_constants=None,
    verbose=False, _print_depth=0,
):
    """Solve a system of GF(2) equations via LU reduction + exhaustive guess.

    Computes a base solution from the LU decomposition, collects
    per-guess effect vectors, prunes impossible monomials (those
    containing a stable zero), and delegates to
    :func:`GuessSolver.guess_and_solve` for linear-independence
    pruning and exhaustive search. Guarantees a solution if one
    exists and the guess space is feasible.

    Unsolved variables and register-state column indices are derived
    from the store's monomial index (``comb_to_idx`` / ``idx_to_comb``)
    and ``solved_for`` vector.

    :param equation_store: An LUEqStore (or any store; converted if needed).
    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param keystream: The observed keystream to verify against.
    :type keystream: np.ndarray[np.uint8]
    :param test_length: Number of keystream bits to use for verification.
    :type test_length: int
    :param additional_constants: Base constant vector for the LU reduce.
        For NAA this contains keystream values at equation positions.
        For FAA/RAA (which use consistency mode) this is None.
    :type additional_constants: np.ndarray[np.uint8] | None
    :param verbose: Whether to print progress.
    :type verbose: bool
    :param _print_depth: Indentation level for verbose output.
    :type _print_depth: int
    :return: ``(initial_state, guesses_tried, pruned_guess_bits)`` —
        the recovered state (or None), total guesses tested, and
        independent guess dimensions after pruning.
    :rtype: tuple[list[int] | None, int, int]
    """
    from PyPR.Cryptanalysis.Components.EquationSolving.GuessSolver import guess_and_solve

    equation_store = _as_lu_store(equation_store)
    num_vars = equation_store.num_vars

    idx_to_comb = equation_store.idx_to_comb
    variable_indices = [equation_store.comb_to_idx[(v,)] for v in range(feedback_fn.size)]
    guess_bits = [(v, idx_to_comb[v]) for v in range(num_vars) if not equation_store.solved_for[v]]

    base_solution = reduce(
        equation_store, additional_constants
    )[variable_indices].copy()

    if verbose:
        print(f"{'|   ' * (_print_depth+1)}Collecting guess effect vectors:")

    effect_collection_start = time.time()
    effects_with_monomials = []
    unstable_bits = np.zeros_like(base_solution)
    for t in range(len(guess_bits)):
        if verbose:
            print(f"\r{'|   ' * (_print_depth+2)}Matrix Solves: {t+1}/{len(guess_bits)}", end='')

        v, comb = guess_bits[t]
        if additional_constants is not None:
            guess_vector = additional_constants.copy()
        else:
            guess_vector = np.zeros(num_vars, dtype=np.uint8)
        guess_vector[v] = 1

        solution = reduce(equation_store, guess_vector)[variable_indices]
        difference = solution ^ base_solution
        effects_with_monomials.append((difference, comb))
        unstable_bits |= difference

    if verbose:
        print(f"\n{'|   ' * (_print_depth+1)}Finished collecting guess effect vectors:")
        print(f"{'|   ' * (_print_depth+1)}Time: {time.time() - effect_collection_start} s")

    # prune impossible monomials: a monomial can't be 1 if it contains a stable zero
    pruned_effects = []
    for effect_vector, comb in effects_with_monomials:
        impossible = False
        for var in comb:
            if (not unstable_bits[var]) and (base_solution[var] == 0):
                impossible = True
                break
        if not impossible:
            pruned_effects.append(effect_vector)

    return guess_and_solve(
        feedback_fn, output_fn, base_solution, pruned_effects,
        keystream, test_length=test_length,
        verbose=verbose, _print_depth=_print_depth,
    )


class LUSolver:
    """Object-oriented wrapper for LU back-substitution solving.

    :param additional_constants: Base constant vector for the LU reduce.
        For NAA this contains keystream values at equation positions.
        For FAA/RAA (which use consistency mode) this is ``None``.
    :type additional_constants: np.ndarray[np.uint8] | None
    """

    def __init__(self, *, additional_constants=None):
        self.additional_constants = additional_constants

    def reduce(self, equation_store):
        return reduce(equation_store, self.additional_constants)

    def solve(
        self, equation_store, feedback_fn, output_fn, keystream, *,
        test_length=1000, verbose=False, _print_depth=0,
    ):
        return solve(
            equation_store, feedback_fn, output_fn, keystream,
            test_length=test_length,
            additional_constants=self.additional_constants,
            verbose=verbose, _print_depth=_print_depth,
        )
