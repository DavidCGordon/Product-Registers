from ProductRegisters import FeedbackRegister
from ProductRegisters.BooleanLogic import BooleanFunction, AND, XOR, CONST

from ProductRegisters.Tools.RootCounting.MonomialProfile import MonomialProfile
from ProductRegisters.Tools.RootCounting.RootExpression import RootExpression
from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet

from ProductRegisters.Cryptanalysis.Components.EquationStores import *
from ProductRegisters.Cryptanalysis.Components.EquationGenerators import *
from ProductRegisters.Cryptanalysis.utils import *

from itertools import product
import numpy as np
import numba
import time


# small helper function to help pretty-print:
def indent(n):
    return ("|   " * n)



def NAA_offline(
    feedback_fn, output_fn,
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
        output_mp = output_fn.remap_constants({
            0: MonomialProfile.logical_zero(),
            1: MonomialProfile.logical_one()
        }).eval_ANF(monomial_profiles)

        if verbose:
            print(f"{indent(print_depth+1)}Monomial profile computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - mp_time} s")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating variable_map:")
            var_map_time = time.time()

        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, output_mp, variable_blocks, complete_subsets = True
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

        eqs = LUDynamicEqStore()
        eq_gen = EqGenerator(feedback_fn, output_fn, 2**feedback_fn.size)

    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Generating Equations:")
        eq_time = time.time()

    # main loop:
    for t, equation, extra_const in eq_gen:
        linearly_independent = eqs.insert_equation(
            equation, extra_const,
            identifier = t,
            # equations are generated in ANF,
            # don't need to translate again
            translate_ANF = False
        )

        if verbose: 
            print(f'\r{indent(print_depth+2)}Equations Found: {eqs.num_eqs} / {eqs.num_vars}',end='')

        # all equations from this point are linearly dependent.
        if not linearly_independent:
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
    output['constant vector'] = eqs.constants[:eqs.num_vars]
    output['keystream needed'] = max(eqs.equation_ids.values()) + 1


    return output



u8 = numba.types.uint8
@numba.njit(u8[:](u8[:,:],u8[:,:],u8[:]))
def lu_solve(L,U,b):
    c = b.copy()

    # backsolve L
    for i in range(len(b)-1):
        for j in range(i+1,len(b)):
            c[j] ^= L[j,i] * c[i]

    # backsolve U
    for i in range(len(b)-1,0,-1):
        for j in range(i):
            c[j] ^= U[j,i] * c[i]

    return c




# Dont need known bits: this is because each equation is cheap (relative to cube attacks)
# and the known bits doesnt /really/ help with the monomials (without a big loop), so it
# doesnt shrink the system that much, but does introduce a lot of overhead.

def NAA_online(feedback_fn, output_fn, keystream, attack_data, test_length = 1000, verbose = False,print_depth=0):
    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Naive Algebraic Attack):")
    start_time = time.time()

    # unpack attack_data
    guess_bits = attack_data['guess vars']
    var_map = attack_data['equation times']
    upper_matrix = attack_data['upper matrix']
    lower_matrix = attack_data['lower matrix']
    const_vector = attack_data['constant vector']
    num_vars = len(upper_matrix)



    if verbose:    
        print(f"{indent(print_depth+1)}Starting initial matrix solve:")

    # initialize new data:
    initial_guess_start = time.time()
    guess_count = 0
    online_vector = np.zeros([num_vars],dtype=np.uint8)

    # determine base solution:
    for v in range(num_vars):
        if v in var_map:
            online_vector[v] = keystream[var_map[v]] ^ const_vector[v]
    
    base_solution = lu_solve(
        lower_matrix,
        upper_matrix,
        online_vector
    )[[attack_data['comb to idx map'][(v,)] for v in range(feedback_fn.size)]].copy()

    # data / buffers for testing an candidate initial state
    F = FeedbackRegister(0,feedback_fn)
    test_length = min(test_length,len(keystream))
    test_keystream = keystream[:test_length]

    # test if this was the correct initial_state
    F._state = base_solution.copy()
    test_seq = [output_fn.eval(state) for state in F.run(test_length)]
    if np.all(test_seq == test_keystream):
        if verbose:
            print(f"{indent(print_depth+1)}Initial matrix solve complete -- correct base solution")
            print(f"{indent(print_depth+1)}Time: {time.time() - initial_guess_start} s")
            print(f"{indent(print_depth)}Online phase complete -- Total time: ", time.time() - start_time)
        return list(base_solution)

    # otherwise we need to try different guesses
    if verbose:
        print(f"{indent(print_depth+1)}Initial solution failed, guessing remaining information:")
        print(f"{indent(print_depth+2)}Collecting guess effect vectors:")

    # first, collect the effects of every guessed bit independently
    effect_collection_start = time.time()
    guess_effect_map = {}
    unstable_bits = np.zeros_like(base_solution)
    for t in range(len(guess_bits)):
        if verbose:
            print(f"\r{indent(print_depth+3)}Matrix Solves: {t+1}/{len(guess_bits)}",end='')

        guess_assignment = [0]*len(guess_bits)
        guess_assignment[t] = 1

        # fill in the vector with known equations + guesses
        for v in range(num_vars):
            if v in var_map:
                online_vector[v] = keystream[var_map[v]] ^ const_vector[v]
        for i, (v,c) in enumerate(guess_bits):
            online_vector[v] = guess_assignment[i]

        # solve the equation.
        solution = lu_solve(
            lower_matrix,
            upper_matrix,
            online_vector
        )[:feedback_fn.size]

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

        F._state = base_solution.copy()
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
 