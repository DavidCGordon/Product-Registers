from PyPR.BooleanLogic import BooleanFunction

from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import CubeEqGenerator, get_var_map
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import SubstitutionEqGenerator

from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver

import numpy as np
import time

# small helper function to help pretty-print:
def indent(n):
    return ("|   " * n)


def RAA_offline(
    feedback_fn, annihilator, multiple, 
    init_rounds, margin,
    time_limit, verbose = False, print_depth=0,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles = None,
    variable_blocks = None
    ):

    if verbose:
        print(f"{indent(print_depth)}Starting offline phase (Reduced Algebraic Attack):")
    start_time = time.time()

    if monomial_profiles != None and variable_blocks != None:
        if verbose:
            print(f"{indent(print_depth+1)}using monomial profile optimization: True")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating larger monomial profile:")
            mp_time = time.time()

        selected = max((annihilator, multiple), key = lambda x: x.degree())
        selected_mp = selected.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)
        max_LC = selected_mp.upper()


        if verbose:
            print(f"{indent(print_depth+1)}Monomial profile computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - mp_time} s")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating variable_map:")
            var_map_time = time.time()
        
        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, selected_mp, variable_blocks, complete_subsets = True
        )

        if verbose:
            print(f"{indent(print_depth+1)}Variable map computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - var_map_time} s")
            #print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating linear relation:")

        # use precomputed maps for faster eq generation and storage
        annihilator_eqs = EqStore(variable_indices)
        multiple_eqs = EqStore(variable_indices)

        check_ranks = False

        eq_gen = CubeEqGenerator(
            feedback_fn, [annihilator, multiple], (len(variable_indices) + margin), 
            variable_indices, #output_map = variable_indices
        )

    else:
        if verbose:
            print(f"{indent(print_depth+1)}Using monomial profile optimization: False")

        # create dynamic storage and generation
        annihilator_eqs = EqStore()
        annihilator_LU = LUEqStore()
        multiple_eqs = EqStore()
        multiple_LU = LUEqStore()

        # link all eq stores
        multiple_eqs.link(annihilator_eqs)
        annihilator_eqs.link(multiple_eqs)
        annihilator_LU.link(annihilator_eqs)
        annihilator_eqs.link(annihilator_LU)
        multiple_LU.link(multiple_eqs)
        multiple_eqs.link(multiple_LU)

        # ensure all variables are in the eq store:
        for v in range(len(feedback_fn)):
            annihilator_eqs._update_known_monomials(tuple([v]))

        eq_gen = SubstitutionEqGenerator(
            feedback_fn, [annihilator, multiple], 2**feedback_fn.size
        )

        # have to check ranks, since number of
        # variables isnt known ahead of time
        check_ranks = True
        count_into_margin = 0


    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Generating Equations:")
        eq_time = time.time()

    # main equation loop
    for t, (ann_eq, mult_eq) in enumerate(eq_gen): # type: ignore (to narrow types correctly)
        ann_eq:  BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]
        mult_eq: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]

        # don't generate equations for initializatipon rounds
        if t < init_rounds: continue

        annihilator_eqs.insert_equation(ann_eq, identifier = t)
        multiple_eqs.insert_equation(mult_eq, identifier = t)

        if verbose: 
            print(f'\r{indent(print_depth+2)}Equations Found: {multiple_eqs.num_eqs} / {multiple_eqs.num_vars + margin}',end='')

        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print(f"\n{indent(print_depth+2)}\n{indent(print_depth+2)}Time limit reached!",end='')
            break

        # break step only necessary for dynamic stores
        if check_ranks:
            ann_indep = annihilator_LU.insert_equation(ann_eq, identifier = t)
            mult_indep = multiple_LU.insert_equation(mult_eq, identifier = t)

            # continue for margin more steps after both have hit their
            # linear recurrence phase (not perfect but better than nothing)
            if not (ann_indep or mult_indep):
                count_into_margin += 1
                if count_into_margin == margin:
                    break


    if verbose:
        #print(f'\r{indent(print_depth+2)}Equations Found: {annihilator_eqs.num_eqs} / {annihilator_eqs.num_vars + margin}',end='\n')
        print(f"\n{indent(print_depth+1)}Finished equation generation: ")
        print(f"{indent(print_depth+1)}Time: {time.time() - eq_time} s")
        print(f"Offline phase complete -- Total time: ", time.time() - start_time)

    output = {}
    output['idx to comb map'] = multiple_eqs.idx_to_comb
    output['comb to idx map'] = multiple_eqs.comb_to_idx
    output['annihilator equations'] = annihilator_eqs.equations[:annihilator_eqs.num_eqs,:annihilator_eqs.num_vars]
    #output['annihilator consts'] = annihilator_eqs.constants[:annihilator_eqs.num_eqs]
    output['multiple equations'] = multiple_eqs.equations[:multiple_eqs.num_eqs,:multiple_eqs.num_vars]
    #output['multiple consts'] = multiple_eqs.constants[:multiple_eqs.num_eqs]
    output['num variables'] =  multiple_eqs.num_vars
    output['keystream needed'] = max(multiple_eqs.equation_ids.values()) + 1
    output['margin'] = margin

    return output




# Dont need known bits: this is because each equation is cheap (relative to cube attacks)
# and the known bits doesnt /really/ help with the monomials (without a big loop), so it
# doesnt shrink the system that much, but does introduce a lot of overhead.
def RAA_online(
    feedback_fn, output_fn, keystream, attack_data,
    test_length=1000, time_limit=None, verbose=False, print_depth=0,
    solver=None, online_store=None,
):
    if solver is None:
        solver = LUSolver()

    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Reduced Algebraic Attack):")
    start_time = time.time()

    # unpack attack_data
    num_vars = attack_data['num variables']
    num_eqs = attack_data['keystream needed']

    annihilator_eqs = attack_data['annihilator equations']
    multiple_eqs = attack_data['multiple equations']

    comb_to_idx = attack_data['comb to idx map']
    idx_to_comb = attack_data['idx to comb map']

    if online_store is None:
        online_store = LUEqStore(comb_to_idx, consistent=(tuple() in comb_to_idx))

    if verbose:
        print(f"{indent(print_depth+1)}Starting Equation Substitution:")

    from PyPR.Cryptanalysis.Components.Adapters.online_insertion import make_online_inserter

    insert_eq, finalize = make_online_inserter(
        online_store, idx_to_comb,
        total_eqs=num_eqs, num_vars=num_vars,
        verbose=verbose, print_depth=print_depth+2,
    )

    for eq_idx in range(num_eqs):
        coef_vector = np.zeros([num_vars], dtype="uint8")
        coef_vector ^= multiple_eqs[eq_idx]
        coef_vector ^= keystream[eq_idx] * annihilator_eqs[eq_idx]

        if insert_eq(coef_vector, eq_idx):
            break
        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print(f"\n{indent(print_depth+1)}Time limit reached during substitution.")
            break

    finalize()

    if verbose:
        if hasattr(online_store, 'solved_vars') and isinstance(online_store.solved_vars, dict):
            solved_count = len(online_store.solved_vars)
        else:
            solved_count = online_store.num_eqs
        print(f"\n{indent(print_depth+1)}Finished substituting key stream:")
        print(f"{indent(print_depth+1)}Variables Solved: {solved_count}/{num_vars}")
        print(f"{indent(print_depth+1)}Time: {time.time() - start_time} s")

    if time_limit and (time.time() - start_time >= time_limit):
        if verbose:
            print(f"{indent(print_depth)}Online phase timed out -- Total time: {time.time() - start_time} s")
        return None

    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Starting solve:")

    initial_state, _, _ = solver.solve(
        online_store, feedback_fn, output_fn, keystream,
        test_length=test_length,
        verbose=verbose, _print_depth=print_depth+1,
    )

    if verbose:
        print(f"{indent(print_depth)}Online phase complete -- Total time: {time.time() - start_time} s")

    return initial_state
 