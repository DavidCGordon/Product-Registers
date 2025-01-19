from ProductRegisters import FeedbackRegister
from ProductRegisters.BooleanLogic import BooleanFunction, AND, XOR, CONST

from ProductRegisters.Tools.RootCounting.MonomialProfile import MonomialProfile
from ProductRegisters.Tools.RootCounting.RootExpression import RootExpression
from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet

from ProductRegisters.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey, berlekamp_massey_iterator

from ProductRegisters.Cryptanalysis.utils import *
from ProductRegisters.Cryptanalysis.Components.EquationGenerators import *
from ProductRegisters.Cryptanalysis.Components.EquationStores import *

from itertools import product
import numpy as np
import numba
import time
import random

from matplotlib import pyplot as plt


# this file is attacks on filter generators - combination generators do not work here. 

# def access_fns(lfsr_fn, filter_fn):
#     # given a full state, simulate that state to get the bit
#     def sim_fn(state):
#         register = FeedbackRegister(state,lfsr_fn)
#         for i in range(init_rounds):
#             register.clock_compiled()
#             output = register._state[0]

#         register.reset()
#         return output


# small helper function to help pretty-print:
def indent(n):
    return ("|   " * n)



def FAA_offline(
    feedback_fn, annihilator, multiple, 
    init_rounds, margin,
    time_limit, 
    verbose = False,
    print_depth = 0,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles = None,
    variable_blocks = None
    ):

    if verbose:
        print(f"{indent(print_depth)}Starting offline phase (Fast Algebraic Attack):")
    start_time = time.time()

    # compute equations for the annihilator:
    if monomial_profiles != None and variable_blocks != None:

        if verbose:
            print(f"{indent(print_depth+1)}using monomial profile optimization: True")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating monomial profile for annihilator:")
            mp_a_time = time.time()

        annihilator_mp = annihilator.remap_constants({
            0: MonomialProfile.logical_zero(),
            1: MonomialProfile.logical_one()
        }).eval_ANF(monomial_profiles)

        if verbose:
            print(f"{indent(print_depth+1)}Monomial profile computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - mp_a_time} s")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating monomial profile for low degree multiple:")
            mp_m_time = time.time()

        # Precompute LC for low degree multiple:
        multiple_mp = multiple.remap_constants({
            0: MonomialProfile.logical_zero(),
            1: MonomialProfile.logical_one()
        }).eval_ANF(monomial_profiles)
        max_LC = multiple_mp.upper()

        if verbose:
            print(f"{indent(print_depth+1)}Monomial profile computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - mp_m_time} s")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating variable_map:")
            var_map_time = time.time()
        
        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, annihilator_mp, variable_blocks, complete_subsets = True
        )

        if verbose:
            print(f"{indent(print_depth+1)}Variable map computed:")
            print(f"{indent(print_depth+1)}Time: {time.time() - var_map_time} s")
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating linear relation:")
            lin_rel_time = time.time()
       
        # use berlekamp_massey to get the exact relation
        feedback_fn.compile()
        test_register = FeedbackRegister(random.random(), feedback_fn)
        seq = [multiple.eval(state) for state in test_register.run_compiled(2*max_LC+32)]
        linear_complexity, linear_relation = berlekamp_massey(seq)

        # flip linear relation, due to dot product vs convolution
        linear_relation = linear_relation[::-1]
        margin += linear_complexity

        if verbose:
            print(f"{indent(print_depth+1)}Linear relation found:")
            print(f"{indent(print_depth+1)}Linear complexity: {linear_complexity}")
            print(f"{indent(print_depth+1)}Time: {time.time()-lin_rel_time} s")

        # use precomputed maps for faster eq generation and storage
        annihilator_eqs = EqStore(variable_indices)
        eq_gen = CubeEqGenerator(
            feedback_fn, annihilator, (len(variable_indices) + margin), 
            variable_indices, verbose=True, print_depth=print_depth+2
        )

        # additional flags
        check_ranks = False
        max_LC = len(variable_indices)

    else:
        if verbose:
            print(f"{indent(print_depth+1)}Using monomial profile optimization: False")

        # create dynamic storage and generation
        annihilator_eqs = DynamicEqStore()
        annihilator_LU = LUDynamicEqStore()
        annihilator_LU.link(annihilator_eqs)

        eq_gen = EqGenerator(
            feedback_fn, annihilator, 2**feedback_fn.size
        )

        # have to check ranks, since number of
        # variables isnt known ahead of time
        check_ranks = True
        count_into_margin = 0

        # Precompute LC for low degree multiple:
        # because max_LC isnt known, test until there are no changes:
        feedback_fn.compile()
        test_register = FeedbackRegister(random.random(), feedback_fn)
        

        if verbose:
            print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating linear complexity dynamically:")
            lin_rel_time = time.time()

        count = 0
        curr_LC = 0
        curr_relation = []
        for linear_complexity, linear_relation in berlekamp_massey_iterator(
            seq = (multiple.eval(state) for state in test_register.run_compiled(2**(feedback_fn.size))),
            yield_rate=1000
        ):
            if verbose:
                print(f"\r{indent(print_depth+2)}Bits processed (thousands): {count} -- Linear Complexity: {curr_LC}", end='')

            # check lengths first for more efficient short circuit:
            if (linear_complexity == curr_LC) and np.all(linear_relation == curr_relation):
                break

            count += 1
            curr_LC = linear_complexity
            curr_relation = linear_relation
        
        # flip linear relation, due to dot product vs convolution
        linear_relation = linear_relation[::-1]
        margin += linear_complexity

        if verbose:
            print(f"\n{indent(print_depth+1)}Linear relation found:")
            print(f"{indent(print_depth+1)}Linear complexity: {curr_LC}")
            print(f"{indent(print_depth+1)}Time: {time.time()-lin_rel_time} s")

    
    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Generating Equations:")
        eq_time = time.time()
   
    # main equation loop
    for t, ann_eq, ann_extra_const in eq_gen:

        # don't generate equations for initialization rounds
        if t < init_rounds: continue
        
        annihilator_eqs.insert_equation(ann_eq, ann_extra_const, identifier = t)

        if verbose: 
            print(f'\r{indent(print_depth+2)}Equations Found: {annihilator_eqs.num_eqs} / {annihilator_eqs.num_vars + margin}',end='')

        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print(f"\n{indent(print_depth)}Time limit reached!")
            break

        # break step only necessary for dynamic stores:
        # reduces speed a fair bit, due to extra insert
        if check_ranks:
            ann_independent = annihilator_LU.insert_equation(ann_eq, ann_extra_const, identifier = t)
        
            # continue for margin more steps after both have hit their
            # linear recurrence phase (not perfect but better than nothing)
            if not (ann_independent):
                count_into_margin += 1
                if count_into_margin == margin:
                    break

    if verbose:
        print(f'\r{indent(print_depth+2)}Equations Found: {annihilator_eqs.num_eqs} / {annihilator_eqs.num_vars + margin}',end='\n')
        print(f"{indent(print_depth+1)}Finished equation generation: ")
        print(f"{indent(print_depth+1)}Time: {time.time() - eq_time} s")
        print(f"Offline phase complete -- Total time: ", time.time() - start_time)

    output = {}
    output['idx to comb map'] = annihilator_eqs.idx_to_comb
    output['comb to idx map'] = annihilator_eqs.comb_to_idx
    output['annihilator equations'] = annihilator_eqs.equations[:annihilator_eqs.num_eqs,:annihilator_eqs.num_vars]
    output['annihilator consts'] = annihilator_eqs.constants[:annihilator_eqs.num_eqs]
    output['linear relation'] = linear_relation
    output['num variables'] = annihilator_eqs.num_vars
    output['keystream needed'] = max(annihilator_eqs.equation_ids.values()) + 1
    output['margin'] = margin - (linear_complexity)

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






@numba.njit(numba.types.Tuple((u8[:],u8))(u64,u8[:],u8[:,:],u8[:],u8[:]))
def sum_over_linear_relationship(start_idx, keystream, equations, constants, linear_relation):
    coef_vector = np.zeros((equations.shape[1],), dtype="uint8")
    const_val = 0

    for i in range(len(linear_relation)):
        if (keystream[start_idx + i])==1 and (linear_relation[i]==1):
            for j in range(equations.shape[1]):
                coef_vector[j] ^= equations[start_idx + i, j]
            const_val ^= constants[start_idx + i]

    return coef_vector,const_val





# Dont need known bits: this is because each equation is cheap (relative to cube attacks)
# and the known bits doesnt /really/ help with the monomials (without a big loop), so it
# doesnt shrink the system that much, but does introduce a lot of overhead.
def FAA_online(
    feedback_fn, output_fn, keystream, attack_data, 
    test_length = 1000, verbose = False, print_depth = 0
):
    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Fast Algebraic Attack):")
    start_time = time.time()

    if type(keystream) != np.ndarray:
        keystream = np.array(keystream, dtype = 'uint8')

    # unpack attack_data
    num_vars = attack_data['num variables']
    num_eqs = attack_data['keystream needed']
    margin = attack_data['margin']

    annihilator_eqs = attack_data['annihilator equations']
    annihilator_consts = attack_data['annihilator consts']
    linear_relation = attack_data['linear relation']
    
    comb_to_idx = attack_data['comb to idx map']
    idx_to_comb = attack_data['idx to comb map']
    variable_indices = [
        attack_data['comb to idx map'][(v,)] 
        for v in range(feedback_fn.size)
    ]

    if verbose:
        print(f"{indent(print_depth+1)}Starting Equation Substitution:")

    # main loop:
    combined_eqs = LUEqStore(comb_to_idx)
    for eq_idx in range(num_eqs - len(linear_relation)):
        coef_vector,const_val = sum_over_linear_relationship(
            eq_idx, keystream, annihilator_eqs, annihilator_consts, linear_relation
        )

        combined_eqs.insert_equation(
            coef_vector, const_val, identifier=eq_idx
        )

        if verbose:
            print(
                f"\r{indent(print_depth+2)}Equations Substituted: {eq_idx+1} / {num_eqs - len(linear_relation)}" +
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
        if v not in combined_eqs.equation_ids
    ]

    # determine base solution:
    base_solution = lu_solve(
        combined_eqs.lower_matrix,
        combined_eqs.upper_matrix,
        combined_eqs.constants
    )[variable_indices].copy()

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
    online_vector = combined_eqs.constants.copy()
    guess_effect_map = {}
    unstable_bits = np.zeros_like(base_solution)
    for t in range(len(guess_bits)):
        if verbose:
            print(f"\r{indent(print_depth+3)}Matrix Solves: {t+1}/{len(guess_bits)}",end='')

        guess_assignment = [0]*len(guess_bits)
        guess_assignment[t] = 1

        # fill in the vector with known equations + guesses
        for v in range(num_vars):
            if v in combined_eqs.equation_ids:
                online_vector[v] = combined_eqs.constants[v]
        for i,(v,comb) in enumerate(guess_bits):
            online_vector[v] = guess_assignment[i]

        # solve the equation.
        solution = lu_solve(
            combined_eqs.lower_matrix,
            combined_eqs.upper_matrix,
            online_vector
        )[variable_indices]

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
    
    
    for i, ((v,comb), effect_vector) in enumerate(guess_effect_map.items()):
        # ignore any guessed_monomial which contains a known 0
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

    # This should never be reached:
    return None
 