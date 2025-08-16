from PyPR.BooleanLogic import BooleanFunction
from PyPR.BooleanLogic.Gates import AND, XOR
from PyPR.BooleanLogic.FunctionInputs import VAR, CONST

from PyPR.Cryptanalysis.Components.EquationStores.DynamicEqStore import DynamicEqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUDynamicEqStore import LUDynamicEqStore

from itertools import product,combinations,chain
import numpy as np
import numba 

import time

# internal use only:
def _generate_monomials(bits,degree=None, verbose=False):
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
def gaussian_elim(matrix):
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

def build_constraint_data(f, candidate_anns, verbose=False):
    dependence_check = LUDynamicEqStore(integrated_constants=True)
    anns = DynamicEqStore(integrated_constants=True)
    mults = DynamicEqStore(integrated_constants=True)

    # build a matching basis for ann/mult constraints:
    for i,candidate in enumerate(candidate_anns):
        if verbose: print(f"\rBuilding Constraints: {i+1}/{len(candidate_anns)}",end='')
        candidate_anf = candidate.translate_ANF()
        linearly_independent = dependence_check.insert_equation(
            candidate_anf,translate_ANF=False,identifier=i
        )

        if linearly_independent:
            anns.insert_equation(candidate_anf,translate_ANF=False,identifier=i)
            mults.insert_equation(
                AND(f,candidate_anf),
                translate_ANF=True,
                identifier=i
            )
    degree_constraints = anns.equations[:anns.num_eqs,:anns.num_vars].T
    zero_constraints = mults.equations[:mults.num_eqs,:mults.num_vars].T
    num_candidates = dependence_check.rank

    return num_candidates, (# all 4 of these just bundled as constraints
        degree_constraints, anns,
        zero_constraints, mults
    )

def ann_solve(constraints, ann_degree, mult_degree):
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
    f, subspace = None,
    annihilator_only=False,
    verbose=True
):
    print("Starting!")
    points = {}

    # allow users to pass no subspace to 
    # use any ann up to the degree of F
    if subspace == None:
        subspace = _generate_monomials(f.idxs_used(),f.degree(), verbose)

    # build contraints:
    num_candidates, constraints = build_constraint_data(
        f, subspace, verbose
    )

    print("\n\nSolving Constraints:\n\n")

    mult_degree = 0
    ann_degree = f.degree()
    while (
        mult_degree <= f.degree() and 
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
            constraints, ann_degree, mult_degree
        )

        # pull out equation ids to recontruct later
        _, ann_eq_store, _, _ = constraints
        ann_eq_ids = ann_eq_store.equation_ids

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
        ann_components = [subspace[v]]
        for row, value in enumerate(reduced_matrix[:,v]):
            if value:
                ann_components.append(subspace[pivots[row]])
        outputs.append(XOR(*(f for f in ann_components)))
    return degrees, outputs

def ann_iterator(
    f, subspace = None,
    annihilator_only=False,
    verbose=True,
    yield_rate = 1
):
    print("Starting!")
    points = {}

    # allow users to pass no subspace to 
    # use any ann up to the degree of F
    if subspace == None:
        subspace = _generate_monomials(f.idxs_used(),f.degree(), verbose)

    # build contraints:
    num_candidates, constraints = build_constraint_data(
        f, subspace, verbose
    )

    print("\n\nSolving Constraints:\n\n\n")

    mult_degree = 0
    ann_degree = f.degree()
    count = 0
    while (
        mult_degree <= f.degree() and 
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
            constraints, ann_degree, mult_degree
        )

        # construct functions from free vars
        # (This is different from main alg, where these values are only constructed at the end)
        ann_list = []
        for v in free_vars:
            ann_components = [subspace[v]]
            for row, value in enumerate(reduced_matrix[:,v]):
                if value:
                    ann_components.append(subspace[pivots[row]])
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
