from PyPR.BooleanLogic import BooleanFunction

from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import CubeEqGenerator, get_var_map
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import SubstitutionEqGenerator

from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.Cryptanalysis.Components.EquationSolving.Grob_Solver import GrobnerSolver

import numpy as np
import time

# small helper function to help pretty-print:
def indent(n):
    return ("|   " * n)

def NAA_offline(
    feedback_fn, output_fn, init_rounds,
    time_limit, verbose = False, print_depth=0,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles = None,
    variable_blocks = None
):
    if verbose:
        print(f"{indent(print_depth)}Starting offline phase (Naive Algebraic Attack):")

    # initialize variables
    start_time = time.time()
    if monomial_profiles != None and variable_blocks != None:
        if verbose:
            print(f"{indent(print_depth+1)}using monomial profile optimization: True")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating larger monomial profile:")
            mp_time = time.time()

        # A map with all subsets filled in, to sum over cubes
        output_mp = output_fn.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)

        if verbose:
            print(f"{indent(print_depth+1)}Monomial profile computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - mp_time} s")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating variable_map:")
            var_map_time = time.time()

        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, output_mp, variable_blocks, complete_subsets = True, include_constant=True
        )

        if verbose:
            print(f"{indent(print_depth+1)}Variable map computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - var_map_time} s")
            #print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating linear relation:")

        # use precomputed maps for faster eq generation and storage
        eqs = LUEqStore(variable_indices)

        check_ranks = False

        eq_gen = CubeEqGenerator(
            feedback_fn, output_fn, 2**feedback_fn.size, 
            variable_indices
        )

    else:
        if verbose:
            print(f"{indent(print_depth+1)}Using monomial profile optimization: False")

        eqs = LUEqStore()
        # ensure all variables are in the eq store:
        for v in range(len(feedback_fn)):
            eqs._update_known_monomials(tuple([v]))

        eq_gen = SubstitutionEqGenerator(feedback_fn, output_fn, 2**feedback_fn.size)

    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Generating Equations:")
        eq_time = time.time()

    # main loop:
    for t, equation in enumerate(eq_gen): #type: ignore  (to narrow types correctly)
        equation: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]

        if t < init_rounds: 
            continue

        linearly_independent = eqs.insert_equation(
            equation,
            identifier = t,
            translate_ANF = False
        )

        if verbose: 
            print(f'\r{indent(print_depth+2)}Equations Found: {eqs.num_eqs} / {eqs.num_vars}',end='')

        if not linearly_independent:
            # all equations from this point are linearly dependent.
            if verbose:
                print(f"\n{indent(print_depth+2)}\n{indent(print_depth+2)}Linear complexity reached!",end='')
            break

        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print(f"\n{indent(print_depth+2)}\n{indent(print_depth+2)}Time limit reached!",end='')
            break

    not_solved = [(x,eqs.idx_to_comb[x]) for x in range(eqs.num_vars) if x not in eqs.equation_ids]

    if verbose:
        #print(f'\r{indent(print_depth+2)}Equations Found: {eqs.num_eqs} / {eqs.num_vars}')
        print(f"\n{indent(print_depth+1)}Finished equation generation: ")
        print(f"{indent(print_depth+1)}Time: {time.time() - eq_time} s")
        print(f"Offline phase complete -- Total time: ", time.time() - start_time)
    
    output = {}
    output['guess vars'] = not_solved
    output['equation times'] = eqs.equation_ids
    output['idx to comb map'] = eqs.idx_to_comb
    output['comb to idx map'] = eqs.comb_to_idx
    output['upper matrix'] = eqs.upper_matrix[:eqs.num_vars,:eqs.num_vars]
    output['lower matrix'] = eqs.lower_matrix[:eqs.num_vars,:eqs.num_vars]
    output['keystream needed'] = max(eqs.equation_ids.values()) + 1

    return output



def NAA_online(
    feedback_fn, output_fn, keystream, attack_data,
    test_length=1000, time_limit=None, verbose=False, print_depth=0,
    solver=None,
):
    if isinstance(solver, GrobnerSolver):
        raise ValueError(
            "NAA with GrobnerSolver is not supported: NAA's offline phase produces "
            "LU matrices with separate constants, which Gröbner solving cannot use. "
            "Use LUSolver or GaussElimSolver instead."
        )

    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Naive Algebraic Attack):")
    start_time = time.time()

    # unpack attack_data
    var_map = attack_data['equation times']
    upper_matrix = attack_data['upper matrix']
    lower_matrix = attack_data['lower matrix']
    num_vars = len(upper_matrix)
    comb_to_idx = attack_data['comb to idx map']

    # reconstruct an LUEqStore from the stored matrices
    solved_store = LUEqStore(comb_to_idx)
    solved_store.upper_matrix[:num_vars,:num_vars] = upper_matrix
    solved_store.lower_matrix[:num_vars,:num_vars] = lower_matrix
    for v in range(num_vars):
        if v in var_map:
            solved_store.solved_for[v] = 1
    solved_store.num_eqs = sum(1 for v in range(num_vars) if v in var_map)
    solved_store.equation_ids = {v: var_map[v] for v in var_map}

    # build constant vector from keystream
    additional_constants = np.zeros([num_vars], dtype=np.uint8)
    for v in range(num_vars):
        if v in var_map:
            additional_constants[v] = keystream[var_map[v]]

    if solver is None:
        solver = LUSolver(additional_constants=additional_constants)
    else:
        solver.additional_constants = additional_constants

    if time_limit and (time.time() - start_time >= time_limit):
        if verbose:
            print(f"{indent(print_depth)}Online phase timed out -- Total time: {time.time() - start_time} s")
        return None

    if verbose:
        print(f"{indent(print_depth+1)}Starting solve:")

    initial_state, _, _ = solver.solve(
        solved_store, feedback_fn, output_fn, keystream,
        test_length=test_length,
        verbose=verbose, _print_depth=print_depth+1,
    )

    if verbose:
        print(f"{indent(print_depth)}Online phase complete -- Total time: {time.time() - start_time} s")

    return initial_state
 