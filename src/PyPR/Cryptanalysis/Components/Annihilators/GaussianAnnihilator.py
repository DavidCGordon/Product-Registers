from PyPR.BooleanLogic import BooleanFunction
from PyPR.BooleanLogic.Gates import AND, XOR
from PyPR.BooleanLogic.FunctionInputs import VAR, CONST

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

from itertools import product,combinations,chain
import numpy as np
import numba 

import time

# internal use only:
def _generate_monomials(bits, degree=None, verbose=False):
    if degree == None:
        degree = len(bits)

    combs = chain.from_iterable(
        combinations(bits, r) 
        for r in range(degree+1)
    )

    count = 0
    output: list["BooleanFunction"] = [CONST(1)]
    for comb in combs:
        if comb != tuple():
            if verbose:
                print(f"\rBuilding Monomials: {count}",end='')
                count += 1

            output.append(AND(*(VAR(v) for v in comb)))
    print("\n")
    return output

u8 = numba.types.uint8
@numba.njit(numba.types.Tuple((u8[:,:],u8[:]))(u8[:,:]))
def gaussian_elim(
    matrix: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
) -> tuple[
    np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    np.ndarray[tuple[int],np.dtype[np.uint8]]
]:
    """Generic GF(2) gaussian elimination method which row reduces and identifies the free variables.

    reduces in place and returns a view of the reduced matrics, so copy the matrix beforehand if you
    dont want it changed by this method. The vector returned for free variables is 1 if the variable
    is free, and 0 otherwise.

    :param matrix: an NxM matrix of uint8's
    :type matrix: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :return: A pair containing a view of the reduced matrix, an an array indicating which variables are free
    :rtype: tuple[np.ndarray[tuple[int,int],np.dtype[np.uint8]],np.ndarray[tuple[int],np.dtype[np.uint8]]]
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

        # Normalize the pivot row
        if matrix[p_row,p_col] == 0:
            free_vars[p_col] = 1
            p_col += 1
            continue

        # Eliminate other rows
        for i in range(rows):
            if i != p_row and matrix[i,p_col]:
                matrix[i] ^= matrix[p_row]

        p_row += 1
        p_col += 1

    # make sure the last columns are counted as free:
    for i in range(p_col,cols):
        free_vars[i] = 1
 
    return matrix[:p_row], free_vars, 

def build_constraint_data(
    input_fn: BooleanFunction, 
    candidate_anns: list[BooleanFunction], 
    verbose: bool = False
) -> tuple[int, tuple[
    np.ndarray[tuple[int,int],np.dtype[np.uint8]], EqStore,
    np.ndarray[tuple[int,int],np.dtype[np.uint8]], EqStore
]]:
    """Build constraints which specify an annihilator/multiple combination of a certain degree.
     
    The constraints are composed by a 4-tuple, which contains both the matrix and eq store for both
    the annihilator and multiple. The end result the type expression is rather nasty,
    but you can safely ignore it mostly, since the entire 4-tuple is passed around as just 
    "contraints". Here is the further information if you need more than that for debugging purposes:

    ### Further Explanation:

    for each candidate, if it is linearly independent with the previously collected candidates,
    it is inserted in the annihilator equations, and its corresponding multiple is inserted in the
    equation store for multiples.

    these (dynamic) equation stores then have each row corresponding to a candidate, with columns
    corresponding the the (dynamically collected) monomials that comprise them. The full constraint
    matrices are formed by taking the relevant portion of the backing matrix for each eq store, and
    transposing it. This way, each column corresponds to the selection of a candiate, and each row
    corresponds to a monomial.

    In `ann_solve` we subselect the rows which correspond to monomials of a prohibitively large degree,
    and then by asserting the matrix product equals the zero vector, we can solve for combinations of
    candidates which zero out all of the monomials with a degree too large.

    :param input_fn: The function we wish to annihilate
    :type input_fn: BooleanFunction
    :param candidate_anns: A basis for the space of annihilators we wish to search
    :type candidate_anns: list[BooleanFunction]
    :param verbose: Whether or not to print updates as the function runs, defaults to False
    :type verbose: bool, optional
    :return: A tuple containing the rank of the annihilator subspace, and the contraints.
    :rtype: tuple[int, tuple[
        np.ndarray[tuple[int,int],np.dtype[np.uint8]], DynamicEqStore,
        np.ndarray[tuple[int,int],np.dtype[np.uint8]], DynamicEqStore
    ]]
    """
    dependence_check = LUEqStore()
    anns = EqStore()
    mults = EqStore()

    # build a matching basis for ann/mult constraints:
    for i,candidate in enumerate(candidate_anns):
        if verbose: 
            print(f"\rBuilding Constraints: {i+1}/{len(candidate_anns)}",end='')

        candidate_anf = candidate.translate_ANF()
        linearly_independent = dependence_check.insert_equation(
            candidate_anf, translate_ANF=False, identifier=i
        )

        if linearly_independent:
            anns.insert_equation(candidate_anf,translate_ANF=False,identifier=i)
            mults.insert_equation(
                AND(input_fn,candidate_anf),
                translate_ANF=True,
                identifier=i
            )
        else:
            print("DEPENDENT!!! BAD!!! ", i,candidate)
    degree_constraints = anns.equations[:anns.num_eqs,:anns.num_vars].T
    zero_constraints = mults.equations[:mults.num_eqs,:mults.num_vars].T
    num_candidates = dependence_check.rank

    return num_candidates, (# all 4 of these just bundled as constraints
        degree_constraints, anns,
        zero_constraints, mults
    )

def ann_solve(
    ann_degree: int, 
    mult_degree: int,
    constraints: tuple[
        np.ndarray[tuple[int,int],np.dtype[np.uint8]], EqStore,
        np.ndarray[tuple[int,int],np.dtype[np.uint8]], EqStore
    ],
) -> tuple[
    list, list,
    np.ndarray[tuple[int,int],np.dtype[np.uint8]]
]:
    """Given the constraints from `build_constrain_data` find the subspace of annihilators
    which meet given degree targets (i.e. annihilator degree at most ann_degree and multiple 
    degree at most mult_degree).

    In `build_constraint_data` the (dynamic) equation stores then have each row corresponding to a 
    candidate, with columns corresponding the the (dynamically collected) monomials that comprise them.
    The full constraint matrices are formed by taking the relevant portion of the backing matrix for 
    each eq store, and transposing it. This way, each column corresponds to the selection of a candiate,
    and each row corresponds to a monomial.

    In this method, we use the comb_to_index maps from the equation stores to subselect the rows which
    correspond to monomials of a prohibitively large degree, and then by concatenating all these rows and
    asserting the matrix product equals the zero vector, we can solve for combinations of candidates which
    zero out all of the monomials with a degree too large in either the annihilator or the multiple.

    :param ann_degree: The maximum degree allowed for the annihilator
    :type ann_degree: int
    :param mult_degree: The maximum degree allowed for the multiple
    :type mult_degree: int
    :param constraints: part of the output of `build_constraint_data`
    :type constraints: tuple[
        np.ndarray[tuple[int,int],np.dtype[np.uint8]], DynamicEqStore, 
        np.ndarray[tuple[int,int],np.dtype[np.uint8]], DynamicEqStore 
    ]
    :return: A tuple containing the pivots, free variables, and reduced matrix, from which a solution
        can be constructed and the annihilators extracted. 
    :rtype: tuple[ list, list, np.ndarray[tuple[int,int],np.dtype[np.uint8]] ]
    """
    ann_constraints,anns_eq_store,mult_constraints,mults_eq_store = constraints
    ann_idxs = anns_eq_store.comb_to_idx
    mult_idxs = mults_eq_store.comb_to_idx

    ann_rows = [i for c,i in ann_idxs.items() if len(c) > ann_degree]
    mult_rows = [i for c,i in mult_idxs.items() if len(c) > mult_degree]
    
    reduced_matrix,free_vars = gaussian_elim(np.concatenate((
        mult_constraints[mult_rows],
        ann_constraints[ann_rows]
    ),axis=0))

    # convert to list of idx, rather than indicator vector
    pivots = [i for i in range(len(free_vars)) if not free_vars[i]]
    free_vars = [i for i in range(len(free_vars)) if free_vars[i]]
    return pivots, free_vars, reduced_matrix

def annihilators(
    input_fn: BooleanFunction,
    subspace: list[BooleanFunction] | None = None,
    annihilator_only: bool = False,
    verbose: bool = False
) -> tuple[
        tuple[int,int],
        list[BooleanFunction]
    ]:
    """Find a basis for the space of optimal annihilators given a set of candidates.

    If no list is passed for the candidates, the full list of monomials is used, and
    thus any function is a viable annihilator candidate.

    ### Degree Walk:
    This function performs a walk on the space of degree pairs `(ann_degree, mult_degree)`.
    for any functions `F`, `F+1` is an annihilator, so `(F.degree,0)` is a known feasible pair,
    and where we begin our walk, and we progress the walk according to simple rules.
    
     - For any feasible pair of degrees, attempt to decrease the annihilator degree by 1
     - For any infeasible pair of degrees, we increase the multiple degree by 1.
    
    We end the walk when `ann_degree` goes below 0 (exiting the grid on the bottom), or 
    `mult_degree` goes above `F.degree()` (as this is worse than just using the base function).
    If `annhilator_only = True`, than we end the walk early, as soon as `mult_degree != 0`, to avoid
    considering non-strict annihilators.
    
    ### Explanation / Proof:
    This is because we have a form of monotonicity in the degrees: as ann degree decreases, 
    the minimum `mult_degree` must be increase or stay the same. The easy proof of this is that our 
    `ann_solve` function allows any degree up to `ann_degree`. Thus all of the solutions possible 
    in the solve for `(ann_degree - 1)` were also viable solutions in the solve for `ann_degree`;
    Since we are looking at a strictly smaller space of solutions, any minimum over the set can only
    increase, including the minimum multiple degree.

    Given that fact, this strategy ensures that any time an encountered pair is feasible, it is the
    lowest multiple degree which is possible at that annihilator degree. Any time we have a feasible
    set of solutions, we also decrease the annihilator degree as much as possible to find the minimum
    annihilator degree for which that multiple degree is still possible. This process allows us to find
    minimal annihilator/multiple degree pairs in a linear number of solves (as we walk along the optimal
    frontier, rather than grid-searching the entire space).

    At the end of the walk, we select the best degree pair, with "best" measured by minimizing
    the maximum degree over the tuple, as this is what will typically most dramatically affect the
    times for algebraic cryptanalysis.

    :param input_fn: The function we are attempting to annihilate.
    :type input_fn: BooleanFunction
    :param subspace: A basis for the subspace of candidate annihilators, or None, 
        if you wish to use the entire function space. defaults to None
    :type subspace: list[BooleanFunction] | None, optional
    :param annihilator_only: Whether to consider all low degree pairs, or only solve
        for strict annihilators, defaults to False
    :type annihilator_only: bool, optional
    :param verbose: Whether to print output as the function runs, defaults to False
    :type verbose: bool, optional
    :return: A tuple containting the degree pair (ann_degree, mult_degree) and basis for the best 
        (i.e. minimal maximum degree) annihilator space (if `annihilator_only = True` then this only
        consider strict annihilators)
    :rtype: list[tuple[int,int], list[BooleanFunction]]
    """
    print("Starting!")
    points = {}

    # allow users to pass no subspace to use any ann up to the degree of F
    if subspace == None:
        subspace = _generate_monomials(
            input_fn.idxs_used(),
            input_fn.degree(), 
            verbose
        )

    # build contraints:
    num_candidates, constraints = build_constraint_data(
        input_fn, subspace, verbose
    )

    # pull out equation ids to recontruct later
    _, ann_eq_store, _, _ = constraints
    ann_eq_ids = ann_eq_store.equation_ids

    print("\n\nSolving Constraints:\n\n")

    mult_degree = 0
    ann_degree = input_fn.degree()
    while (
        mult_degree <= input_fn.degree() and 
        ann_degree >= 0
    ):
        if verbose: print(
            "\r\x1B[2A" + 
            f"|   Annihilator Degree: {ann_degree}\n" + 
            f"|   Multiple Degree: {mult_degree}\n"
            ,end=''
        )

        # solve system using contraints   
        pivots, free_vars, reduced_matrix = ann_solve(
            ann_degree, mult_degree, constraints
        )

        # update the degrees
        if free_vars:
            points[(ann_degree,mult_degree)] = pivots, free_vars, reduced_matrix
            ann_degree -= 1
        else:
            mult_degree += 1
            if annihilator_only:
                break
 
    if verbose:
        print("\nPOINTS: ", points.keys())

    selected = min(points.items(), key = lambda x: sorted(x[0],reverse=True))
    pivots, free_vars, reduced_matrix = selected[1]
    degrees = selected[0]

    outputs = []
    for v in free_vars:
        # Free var at index v => use eq_ids to get the original candidate
        # Row i in reduced => row pivots[i] in original
        ann_components = [subspace[ann_eq_ids[v]]]
        for row, value in enumerate(reduced_matrix[:,v]):
            if value:
                adjusted_idx = pivots[row]
                subspace_idx = ann_eq_ids[adjusted_idx]
                ann_components.append(subspace[subspace_idx])
        outputs.append(XOR(*(f for f in ann_components)))
    return degrees, outputs

def ann_iterator(
    input_fn: BooleanFunction, 
    subspace: list[BooleanFunction] | None = None,
    annihilator_only: bool = False,
    verbose: bool = True,
    yield_rate: int = 1
):
    """Like `annihilators`, this function finds a basis for the space of optimal annihilators
    given a set of candidates, but it also yields during the degree walk, allowing  more insight
    into the algorithm as is runs. It also constructs the BooleanFunction basis from the matrix
    for all degree pairs, not just the selected on at the end.

    If no list is passed for the candidates, the full list of monomials is used, and
    thus any function is a viable annihilator candidate.

    ### Degree Walk:
    This function performs a walk on the space of degree pairs `(ann_degree, mult_degree)`.
    for any functions `F`, `F+1` is an annihilator, so `(F.degree,0)` is a known feasible pair,
    and where we begin our walk, and we progress the walk according to simple rules.
    
     - For any feasible pair of degrees, attempt to decrease the annihilator degree by 1
     - For any infeasible pair of degrees, we increase the multiple degree by 1.
    
    We end the walk when `ann_degree` goes below 0 (exiting the grid on the bottom), or 
    `mult_degree` goes above `F.degree()` (as this is worse than just using the base function).
    If `annhilator_only = True`, than we end the walk early, as soon as `mult_degree != 0`, to avoid
    considering non-strict annihilators.
    
    ### Explanation / Proof:
    This is because we have a form of monotonicity in the degrees: as ann degree decreases, 
    the minimum `mult_degree` must be increase or stay the same. The easy proof of this is that our 
    `ann_solve` function allows any degree up to `ann_degree`. Thus all of the solutions possible 
    in the solve for `(ann_degree - 1)` were also viable solutions in the solve for `ann_degree`;
    Since we are looking at a strictly smaller space of solutions, any minimum over the set can only
    increase, including the minimum multiple degree.

    Given that fact, this strategy ensures that any time an encountered pair is feasible, it is the
    lowest multiple degree which is possible at that annihilator degree. Any time we have a feasible
    set of solutions, we also decrease the annihilator degree as much as possible to find the minimum
    annihilator degree for which that multiple degree is still possible. This process allows us to find
    minimal annihilator/multiple degree pairs in a linear number of solves (as we walk along the optimal
    frontier, rather than grid-searching the entire space).

    At the end of the walk, we select the best degree pair, with "best" measured by minimizing
    the maximum degree over the tuple, as this is what will typically most dramatically affect the
    times for algebraic cryptanalysis.

    :param input_fn: The function we are attempting to annihilate.
    :type input_fn: BooleanFunction
    :param subspace: A basis for the subspace of candidate annihilators, or None, 
        if you wish to use the entire function space. defaults to None
    :type subspace: list[BooleanFunction] | None, optional
    :param annihilator_only: Whether to consider all low degree pairs, or only solve
        for strict annihilators, defaults to False
    :type annihilator_only: bool, optional
    :param verbose: Whether to print output as the function runs, defaults to False
    :type verbose: bool, optional
    :return: A tuple containting the degree pair (ann_degree, mult_degree) and basis for the best 
        (i.e. minimal maximum degree) annihilator space (if `annihilator_only = True` then this only
        consider strict annihilators)
    :rtype: list[tuple[int,int], list[BooleanFunction]]
    """
    print("Starting!")
    points = {}

    # allow users to pass no subspace to use any ann up to the degree of F
    if subspace == None:
        subspace = _generate_monomials(
            input_fn.idxs_used(),
            input_fn.degree(), 
            verbose
        )

    # build contraints:
    num_candidates, constraints = build_constraint_data(
        input_fn, subspace, verbose
    )

    # pull out equation ids to recontruct later
    _, ann_eq_store, _, _ = constraints
    ann_eq_ids = ann_eq_store.equation_ids


    print("\n\nSolving Constraints:\n\n\n")

    mult_degree = 0
    ann_degree = input_fn.degree()
    count = 0
    while (
        mult_degree <= input_fn.degree() and 
        ann_degree >= 0
    ):
        count += 1
        if verbose: print(
            "\r\x1B[3A" + 
            f"|   Iteration: {count}\n" + 
            f"|   Annihilator Degree: {ann_degree}\n" + 
            f"|   Multiple Degree: {mult_degree}\n",
            end=''
        )

        # solve system using contraints   
        pivots, free_vars, reduced_matrix = ann_solve(
            ann_degree, mult_degree, constraints
        )

        # construct functions from free vars
        # (This is different from main alg, where these values are only constructed at the end)
        ann_list = []
        for v in free_vars:
            ann_components = [subspace[ann_eq_ids[v]]]
            for row, value in enumerate(reduced_matrix[:,v]):
                if value:
                    adjusted_idx = pivots[row]
                    subspace_idx = ann_eq_ids[adjusted_idx]
                    ann_components.append(subspace[subspace_idx])
            ann_list.append(XOR(*(f for f in ann_components)))

        # update the degrees
        if free_vars:
            points[(ann_degree,mult_degree)] = ann_list
            ann_degree -= 1
        else:
            mult_degree += 1
            if annihilator_only:
                break

        # for convenience, calculate the best item. 
        # (This is different from main alg, where these values are only constructed at the end)
        selected = min(points.items(), key = lambda x: sorted(x[0],reverse=True))
        
        # yield
        if count % yield_rate == 0:
            yield (
                (ann_degree,mult_degree),
                selected,
                points,
            )
 
    if verbose:
        print("\nPOINTS: ", points.keys())

    selected = min(points.items(), key = lambda x: sorted(x[0],reverse=True))
    return selected
