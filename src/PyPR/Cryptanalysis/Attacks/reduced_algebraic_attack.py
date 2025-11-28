from PyPR import FeedbackRegister
from PyPR.BooleanLogic import BooleanFunction

from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import CubeEqGenerator, get_var_map
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import SubstitutionEqGenerator

import PyPR.Cryptanalysis.Components.EquationSolvers.LU_Solver as LU_Solver

from itertools import product
import numpy as np
import numba
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
def RAA_online(feedback_fn, output_fn, keystream, attack_data, test_length = 1000, verbose = False, print_depth = 0):
    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Reduced Algebraic Attack):")
    start_time = time.time()

    # unpack attack_data
    num_vars = attack_data['num variables']
    num_eqs = attack_data['keystream needed']
    margin = attack_data['margin']

    annihilator_eqs = attack_data['annihilator equations']
    multiple_eqs = attack_data['multiple equations']
    
    comb_to_idx = attack_data['comb to idx map']
    idx_to_comb = attack_data['idx to comb map']
    variable_indices = [
        attack_data['comb to idx map'][(v,)] 
        for v in range(feedback_fn.size)
    ]

    if verbose:
        print(f"{indent(print_depth+1)}Starting Equation Substitution:")

    # main loop:
    combined_eqs = LUEqStore(comb_to_idx, consistent=(tuple() in comb_to_idx))
    for eq_idx in range(num_eqs):
        # initialize vector/const:
        coef_vector = np.zeros([num_vars], dtype="uint8")

        # use keystream to construct final equation and const:
        coef_vector ^= multiple_eqs[eq_idx]
        coef_vector ^= keystream[eq_idx] * annihilator_eqs[eq_idx]
        combined_eqs.insert_equation(coef_vector, identifier=eq_idx)

        if verbose:
            print(
                f"\r{indent(print_depth+2)}Equations Substituted: {eq_idx+1} / {num_eqs}" +
                f"  --  Current Rank: {combined_eqs.rank} / {num_vars}",
                end = ''
            )

        if combined_eqs.rank == combined_eqs.num_vars:
            if verbose:
                print(f"\n{indent(print_depth+2)}\n{indent(print_depth+2)}Substitution finished early!")
                print(f"{indent(print_depth+2)}Equations processed: {eq_idx+1}/{num_eqs}", end='')
            break

    if verbose:
        print(f"\n{indent(print_depth+1)}Finished substituting key stream:")
        print(f"{indent(print_depth+1)}Variables Solved: {combined_eqs.num_eqs}/{num_vars}")
        print(f"{indent(print_depth+1)}Time: {time.time() - start_time} s")    
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Starting initial matrix solve:")

    # initialize new data for guessing:
    initial_guess_start = time.time()
    guess_count = 0
    guess_bits = [
        (v, idx_to_comb[v]) for v in range(num_vars)
        if not combined_eqs.solved_for[v]
    ]

    # determine base solution:
    base_solution = LU_Solver.solve(combined_eqs)[variable_indices].copy()

    # data / buffers for testing an candidate initial state
    F = FeedbackRegister(0,feedback_fn)
    F.set_state(base_solution)

    test_length = min(test_length,len(keystream))
    test_keystream = keystream[:test_length]
    test_seq = [output_fn.eval(state) for state in F.run(test_length)]
    
    # test if this was the correct initial_state
    if np.all(test_seq == test_keystream):
        if verbose:
            print(f"{indent(print_depth+1)}Initial matrix solve complete -- correct base solution")
            print(f"{indent(print_depth+1)}Time: {time.time() - initial_guess_start} s")
            print(f"{indent(print_depth)}Online phase complete -- Total time: ", time.time() - start_time)
        return base_solution
    
    # otherwise we need to try different guesses
    if verbose:
        print(f"{indent(print_depth+1)}Initial solution failed, guessing remaining information:")
        print(f"{indent(print_depth+2)}Collecting guess effect vectors:")

    # first, collect the effects of every guessed bit independently
    effect_collection_start = time.time()
    guess_vector = np.zeros(combined_eqs.num_vars, dtype = np.uint8)
    unstable_bits = np.zeros_like(base_solution)
    guess_effect_map = {}
    for t in range(len(guess_bits)):
        if verbose:
            print(f"\r{indent(print_depth+3)}Matrix Solves: {t+1}/{len(guess_bits)}",end='')

        # fill in guess vector
        guess_assignment = [0]*len(guess_bits)
        guess_assignment[t] = 1
        for i,(v,comb) in enumerate(guess_bits):
            guess_vector[v] = guess_assignment[i]

        # solve the equation.
        solution = LU_Solver.solve(combined_eqs, guess_vector)[variable_indices]

        difference = (solution ^ base_solution)
        guess_effect_map[guess_bits[t]] = difference
        unstable_bits |= difference

    if verbose:
        print(f"\n{indent(print_depth+2)}Finished collecting guess effect vectors:")
        print(f"{indent(print_depth+2)}Time: {time.time() - effect_collection_start} s")
        print(f"{indent(print_depth+2)}\n{indent(print_depth+2)}Starting effect pruning:")

    # prune guesses by removing impossible and dependent guesses:
    effect_pruning_time = time.time()
    pruned_guesses = []
    already_solved = set()
    reduced_matrix = np.zeros([feedback_fn.size,feedback_fn.size], dtype = np.uint8)

    for (v,comb), effect_vector in guess_effect_map.items():
        # don't guess a monomial which contains a known 0
        impossible_comb = False
        for var in comb:
            if (not unstable_bits[var]) and (base_solution[var] == 0):
                impossible_comb = True
        
        if impossible_comb:
            continue

        # check that this effect vector is linearly independent:
        # note that this is a slightly simplified LU build-up, with
        # reduced_mat = upper_matrix, and already_solved = var_map
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
            print(f"\r{indent(print_depth+3)}Vectors Pruned: {i+1}/{len(guess_effect_map)}",end='')

    if verbose:
        print(f"\n{indent(print_depth+2)}Pruning finished:")
        print(f"{indent(print_depth+2)}Max number of guesses (original): 2^{len(guess_bits)}")
        print(f"{indent(print_depth+2)}Max number of guesses (pruned): 2^{len(pruned_guesses)}")
        print(f"{indent(print_depth+2)}Time: {time.time() - effect_pruning_time} s")
        #print(f"{indent(print_depth+2)}Time for total pruning process: {time.time() - effect_collection_start} s")
        print(f"{indent(print_depth+2)}\n{indent(print_depth+2)}Starting to Guess:")
    
    # Now test using the pruned guesses:
    guess_count = 0
    guess_start_time = time.time()
    for guess_assignment in product((0,1), repeat = len(pruned_guesses)):
        guess_count += 1

        if verbose:
            print(f"\r{indent(print_depth+3)}Guess count: {guess_count}",end='')

        F.set_state(base_solution)
        for idx, assigned in enumerate(guess_assignment):
            if assigned:
                F._state ^= pruned_guesses[idx]
        F.seed(F._state.copy())

        mismatch = False
        for t,state in enumerate(F.run(test_length)):
            if output_fn.eval(state) != test_keystream[t]:
                mismatch = True
                break
        
        if not mismatch:
            if verbose:
                print(f"\n{indent(print_depth+2)}Guessing Finished:")
                print(f"{indent(print_depth+2)}Time: {time.time() - guess_start_time} s")
                print(f"{indent(print_depth+1)}Solution Found!")
                print(f"{indent(print_depth)}Online phase complete -- Total time: ", time.time() - start_time)
            F.reset()
            return list(F)
    return None
 