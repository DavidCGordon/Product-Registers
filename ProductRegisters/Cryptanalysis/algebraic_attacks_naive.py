from ProductRegisters import FeedbackRegister
from ProductRegisters.BooleanLogic import BooleanFunction, AND, XOR, CONST

from itertools import product
import numpy as np
import numba
import time

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

def alg_attack_offline(feedback_fn, output_fn, time_limit, verbose = False):
    start_time = time.time()
    idx_to_comb = {}    # variable idx -> monomial
    comb_to_idx = {}    # monomial -> variable idx
    var_map = {}        # variable idx -> clock cycle which gives eq responsible for it

    # initialize index maps:
    for v in range(feedback_fn.size):
        idx_to_comb[v] = (v,)
        comb_to_idx[(v,)] = v

    # initialize var counts:
    num_vars = feedback_fn.size
    max_vars = 1
    while max_vars <= num_vars:
        max_vars *= 2

    # initialize matrices:
    upper_matrix = np.eye(max_vars, dtype=np.uint8)
    lower_matrix = np.eye(max_vars, dtype=np.uint8)
    const_vec = np.zeros([max_vars,1],dtype=np.uint8)
    
    # main loop
    for t,anfs in enumerate(feedback_fn.anf_iterator(2**feedback_fn.size,bits = output_fn.idxs_used())):
        if verbose: 
            print(f'\rEquations Found: {t} / {num_vars}',end='')
        
        equation_anf = output_fn.compose(anfs).translate_ANF()

        # find the equation vector from anf:
        const_val = 0
        coef_vector = np.zeros([max_vars], dtype=np.uint8)
        for term in equation_anf.args:

            # handle constant values
            if type(term) == CONST:
                const_val = term.value
                continue

            comb = tuple(sorted([var.index for var in term.args]))

            # expand matrices as necessary:
            if comb not in comb_to_idx:
                if num_vars == max_vars:
                    max_vars *= 2

                    # initialize
                    new_upper_matrix = np.eye(max_vars, dtype=np.uint8)
                    new_lower_matrix = np.eye(max_vars, dtype=np.uint8)
                    new_const_vec = np.zeros([max_vars,1],dtype=np.uint8)
                    new_coef_vector = np.zeros([max_vars], dtype=np.uint8)

                    # copy
                    new_upper_matrix[:num_vars, :num_vars] = upper_matrix
                    new_lower_matrix[:num_vars, :num_vars] = lower_matrix
                    new_const_vec[:num_vars] = const_vec
                    new_coef_vector[:num_vars] = coef_vector

                    # replace
                    upper_matrix = new_upper_matrix
                    lower_matrix = new_lower_matrix
                    const_vec = new_const_vec
                    coef_vector = new_coef_vector

                idx_to_comb[num_vars] = comb
                comb_to_idx[comb] = num_vars
                num_vars += 1

            # after expanding as necessary, still set the appropriate var:
            coef_vector[comb_to_idx[comb]] = 1

        linearly_independent = False
        modification_vector  = np.zeros_like(coef_vector)
        for idx in range(len(coef_vector)):
            if coef_vector[idx] == 1:
                modification_vector[idx] = 1
                if idx in var_map:
                    coef_vector ^= upper_matrix[idx]
                else:
                    linearly_independent = True

                    var_map[idx] = t
                    lower_matrix[idx] = modification_vector
                    upper_matrix[idx] = coef_vector
                    const_vec[idx] = const_val
                    break

        # the first time you get a nonlinear equation you have hit LC
        # all equations from this point are NL.
        if not linearly_independent:
            if verbose:
                print("\nLinear complexity reached! ")
            break

        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print("\nTime limit reached!")
            break

    not_solved = [(x,idx_to_comb[x]) for x in range(num_vars) if x not in var_map]

    if verbose:
        print(f"Final number of variables: {num_vars}")
        print(f"Final number of equations: {len(var_map)}")
        print(len(not_solved))
        print(f"Keystream Required: {max(var_map.values()) + 1} bits")

    
    output = {}
    output['guess vars'] = not_solved
    output['equation times'] = var_map
    output['idx to comb map'] = idx_to_comb
    output['comb to idx map'] = comb_to_idx
    output['upper matrix'] = upper_matrix[:num_vars,:num_vars]
    output['lower matrix'] = lower_matrix[:num_vars,:num_vars]
    output['constant vector'] = const_vec[:num_vars]
    output['keystream needed'] = max(var_map.values()) + 1

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

def alg_attack_online(feedback_fn, output_fn, keystream, attack_data, test_length = 1000, verbose = False):
    start_time = time.time()

    # unpack attack_data
    guess_bits = attack_data['guess vars']
    var_map = attack_data['equation times']
    upper_matrix = attack_data['upper matrix']
    lower_matrix = attack_data['lower matrix']
    const_vector = attack_data['constant vector']
    num_vars = len(upper_matrix)

    # initialize new data:
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
    )[:feedback_fn.size].copy()

    # data / buffers for testing an candidate initial state
    F = FeedbackRegister(0,feedback_fn)
    test_length = min(test_length,len(keystream))
    test_keystream = keystream[:test_length]

    # test if this was the correct initial_state
    F._state = base_solution.copy()
    test_seq = [output_fn.eval(state) for state in F.run(test_length)]
    if np.all(test_seq == test_keystream):
        if verbose:
            print(f'Correct Base Solution  --  total time:', time.time() - start_time)
        return list(base_solution)
    
    if verbose:
        print("Initial solution failed, guessing remaining information:")

    # otherwise we need to try different guesses
    # first, collect the effects of every guessed bit independently
    effect_collection_start = time.time()
    guess_effect_map = {}
    unstable_bits = np.zeros_like(base_solution)
    for t in range(len(guess_bits)):
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
        print(f"Finished collecting guess effect vectors:")
        print('Matrix solves:', len(guess_bits), '\tTotal time: ', time.time() - effect_collection_start)

    # prune guesses by removing impossible and dependent guesses:
    effect_pruning_time = time.time()
    pruned_guesses = []
    already_solved = set()
    reduced_matrix = np.zeros([feedback_fn.size,feedback_fn.size], dtype = np.uint8)
    for (v,comb), effect_vector in guess_effect_map.items():
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
        print(f"Pruning finished: {time.time() - effect_pruning_time} s")
        print(f"Max number of guesses (original): 2^{len(guess_bits)}")
        print(f"Max number of guesses (pruned): 2^{len(pruned_guesses)}")

    # Now test using the pruned guesses:
    guess_count = 0
    guess_start_time = time.time()
    for guess_assignment in product((0,1), repeat = len(pruned_guesses)):
        guess_count += 1

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
                print('Guess count:', guess_count, '\tTotal guess time: ', time.time() - guess_start_time)
            F.reset()
            return list(F)
    return None
 