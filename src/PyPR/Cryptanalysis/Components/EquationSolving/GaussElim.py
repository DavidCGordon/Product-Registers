import numpy as np
import numba
import time

from PyPR.Cryptanalysis.Components.Adapters.store_repr import to_coef_matrix

u8 = numba.types.uint8
@numba.njit(numba.types.Tuple((u8[:,:],u8[:]))(u8[:,:]))
def reduce_matrix(matrix):
    """Compute the reduced row echelon form (RREF) over GF(2).

    :param matrix: Binary coefficient matrix (uint8, entries 0 or 1).
        Modified in place.
    :type matrix: np.ndarray
    :return: ``(rref, free_vars)`` — the RREF (truncated to pivot rows)
        and a binary vector marking free (non-pivot) columns.
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    rows, cols = matrix.shape
    p_row, p_col = 0, 0
    free_vars = np.zeros((cols,),dtype='uint8')

    while p_row < rows and p_col < cols:
        # Find the pivot element/swap rows
        for i in range(p_row + 1, rows):
            if matrix[i,p_col] > matrix[p_row,p_col]:
                matrix[np.array([p_row, i])] = matrix[np.array([i, p_row])]
                break

        # Identify free vars
        if matrix[p_row,p_col] == 0:
            free_vars[p_col] = 1
            p_col += 1
            continue

        # No normalization needed due to 0/1 only

        # Eliminate other rows
        for i in range(rows):
            if i != p_row and matrix[i,p_col]:
                matrix[i] ^= matrix[p_row]

        p_row += 1
        p_col += 1

    # make sure the last columns are counted as free:
    for i in range(p_col,cols):
        free_vars[i] = 1

    return matrix[:p_row], free_vars


def reduce(equation_store, constants=None):
    """Reduce a system of GF(2) equations via Gaussian elimination (RREF).

    Produces a solution vector from the RREF alone, without guessing
    free variables. If the system is underdetermined, free variables
    default to zero.

    Accepts any equation store. Non-EqStore types are converted via
    ``to_coef_matrix`` (expensive but valid).

    :param equation_store: Any equation store.
    :param constants: Unused (reserved for future augmented-column support).
    :return: Solution vector of length ``num_vars``.
    :rtype: np.ndarray
    """
    if (
        hasattr(equation_store, 'equations') and
        isinstance(equation_store.equations, np.ndarray) and
        hasattr(equation_store, 'comb_to_idx')
    ):
        matrix = equation_store.equations[
            :equation_store.num_eqs,
            :equation_store.num_vars
        ].copy()
        num_vars = equation_store.num_vars
    else:
        matrix, _, _ = to_coef_matrix(equation_store)
        num_vars = matrix.shape[1]

    rref, free_vars = reduce_matrix(matrix)
    num_pivot_rows = rref.shape[0]

    solution = np.zeros(num_vars, dtype=np.uint8)
    row = num_pivot_rows - 1
    col = num_vars - 1
    while row >= 0 and col >= 0:
        if free_vars[col]:
            col -= 1
            continue

        if rref[row, col]:
            for i in range(row - 1, -1, -1):
                solution[i] ^= rref[i, col]

        row -= 1
        col -= 1

    # For consistent stores, CONST is always 1
    if (
        hasattr(equation_store, 'consistent') and
        equation_store.consistent and
        hasattr(equation_store, 'comb_to_idx') and
        tuple() in equation_store.comb_to_idx
    ):
        solution[equation_store.comb_to_idx[tuple()]] = 1

    return solution


def solve(
    equation_store, feedback_fn, output_fn, keystream, *,
    test_length=1000, additional_constants=None,
    verbose=False, _print_depth=0,
):
    """Solve a system of GF(2) equations via RREF + exhaustive guess.

    Computes the RREF once, extracts a base solution and per-guess
    effect vectors directly from the reduced matrix, prunes impossible
    monomials, and delegates to :func:`GuessSolver.guess_and_solve`
    for linear-independence pruning and exhaustive search.

    Unsolved variables are identified from the RREF's free-variable
    vector; register-state column indices are derived from the store's
    monomial index (``comb_to_idx``).

    .. note::

        This is dramatically less efficient than
        :func:`LU_Solver.solve`, which reuses a cached LU
        decomposition for O(n^2) back-substitution per guess bit.
        Here the RREF itself costs O(n^3), and while effect vectors
        are read from the RREF in O(n), the initial reduction
        dominates. Prefer ``LU_Solver`` for attack pipelines;
        this function exists for completeness and for stores that
        cannot be converted to LUEqStore cheaply.

    :param equation_store: Any equation store.
    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function for keystream generation.
    :type output_fn: BooleanFunction
    :param keystream: The observed keystream to verify against.
    :type keystream: np.ndarray[np.uint8]
    :param test_length: Number of keystream bits to use for verification.
    :type test_length: int
    :param additional_constants: Base constant vector of length
        ``num_vars``. For LUEqStore inputs, entries align by variable
        index with the upper_matrix rows.
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

    # Extract augmented matrix [A | b] and monomial index maps
    if (
        additional_constants is not None and
        hasattr(equation_store, 'upper_matrix') and
        hasattr(equation_store, 'solved_for')
    ):
        n = equation_store.num_vars
        idx_to_comb = equation_store.idx_to_comb
        comb_to_idx = equation_store.comb_to_idx

        # Forward-solve additional_constants through L to get the
        # correct RHS for each reduced (upper-matrix) equation.
        L = equation_store.lower_matrix[:n, :n]
        transformed = additional_constants[:n].copy()
        for i in range(n - 1):
            for j in range(i + 1, n):
                transformed[j] ^= L[j, i] * transformed[i]

        # Extract only pivot rows; skip the const axiom row for consistent stores
        const_idx = comb_to_idx.get(tuple())
        pivot_rows = [i for i in range(n) if equation_store.solved_for[i]]
        if equation_store.consistent and const_idx is not None:
            pivot_rows = [i for i in pivot_rows if i != const_idx]

        if pivot_rows:
            matrix = np.stack([equation_store.upper_matrix[i, :n] for i in pivot_rows])
            const_col = np.array(
                [transformed[i] for i in pivot_rows], dtype=np.uint8
            ).reshape(-1, 1)
            # Restore zeroed-out constants for post-const_idx rows
            if equation_store.consistent and const_idx is not None:
                for row_idx, i in enumerate(pivot_rows):
                    if i > const_idx:
                        matrix[row_idx, const_idx] = equation_store._extra_constants[i]
        else:
            matrix = np.zeros((0, n), dtype=np.uint8)
            const_col = np.zeros((0, 1), dtype=np.uint8)

        augmented = np.hstack([matrix, const_col])
        num_vars = n
    elif (
        hasattr(equation_store, 'equations') and
        isinstance(equation_store.equations, np.ndarray) and
        hasattr(equation_store, 'comb_to_idx')
    ):
        matrix = equation_store.equations[
            :equation_store.num_eqs,
            :equation_store.num_vars
        ].copy()
        num_vars = equation_store.num_vars
        idx_to_comb = equation_store.idx_to_comb
        comb_to_idx = equation_store.comb_to_idx
        const_col = np.zeros((matrix.shape[0], 1), dtype=np.uint8)
        augmented = np.hstack([matrix, const_col])
    else:
        matrix, comb_to_idx, idx_to_comb = to_coef_matrix(equation_store)
        num_vars = matrix.shape[1]
        const_col = np.zeros((matrix.shape[0], 1), dtype=np.uint8)
        augmented = np.hstack([matrix, const_col])

    variable_indices = [comb_to_idx[(v,)] for v in range(feedback_fn.size)]

    rref, free_vars_aug = reduce_matrix(augmented)
    num_pivot_rows = rref.shape[0]

    # separate coefficient columns from the constant column
    rref_coeffs = rref[:, :num_vars]
    rref_constants = rref[:, -1]
    coeff_free = free_vars_aug[:num_vars]

    # build pivot map: pivot_cols[r] = pivot column of row r
    pivot_cols = []
    c = 0
    for r in range(num_pivot_rows):
        while c < num_vars and coeff_free[c]:
            c += 1
        if c < num_vars:
            pivot_cols.append(c)
            c += 1
        else:
            break

    # base solution: x[p] = constant[r] for each pivot (r, p), free vars = 0
    base_full = np.zeros(num_vars, dtype=np.uint8)
    for r, p in enumerate(pivot_cols):
        base_full[p] = rref_constants[r]
    base_solution = base_full[variable_indices].copy()

    # guess_bits: free columns with their monomial tuples
    guess_bits = [(v, idx_to_comb[v]) for v in range(num_vars) if coeff_free[v]]

    if verbose:
        print(f"{'|   ' * (_print_depth+1)}Computing effect vectors from RREF:")

    # effect vectors: for free variable v, flipping it changes each
    # pivot variable p by rref_coeffs[r, v] (read directly from RREF)
    effect_collection_start = time.time()
    effects_with_monomials = []
    unstable_bits = np.zeros_like(base_solution)
    for t, (v, comb) in enumerate(guess_bits):
        if verbose:
            print(f"\r{'|   ' * (_print_depth+2)}Effect Vectors: {t+1}/{len(guess_bits)}", end='')

        effect_full = np.zeros(num_vars, dtype=np.uint8)
        effect_full[v] = 1
        for r, p in enumerate(pivot_cols):
            effect_full[p] ^= rref_coeffs[r, v]

        effect = effect_full[variable_indices]
        effects_with_monomials.append((effect, comb))
        unstable_bits |= effect

    if verbose:
        print(f"\n{'|   ' * (_print_depth+1)}Finished computing effect vectors:")
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


class GaussElimSolver:
    """Object-oriented wrapper for Gaussian elimination (RREF) solving.

    :param additional_constants: Base constant vector for the RREF.
        For NAA this contains keystream values at equation positions.
        For FAA/RAA (which use consistency mode) this is ``None``.
    :type additional_constants: np.ndarray[np.uint8] | None
    """

    def __init__(self, *, additional_constants=None):
        self.additional_constants = additional_constants

    def reduce(self, equation_store):
        return reduce(equation_store)

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
