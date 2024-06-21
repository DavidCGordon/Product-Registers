from ProductRegisters.BooleanLogic import BooleanFunction, AND, XOR, CONST

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

# establish matrix indexing

def alg_attack_offline(feedback_fn, output_fn, time_limit):
    count = 0
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
    time_limit = None
    
    for t,anfs in enumerate(feedback_fn.anf_iterator(2**feedback_fn.size,bits = output_fn.idxs_used())):
        print(t, end = ", ", flush=True)
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

            # if monomial already has a variable, simply set it:
            if comb in comb_to_idx:
                coef_vector[comb_to_idx[comb]] = 1 

            # otherwise, add new variable entries, then set it
            else:
                #expand matrices as necessary:
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
            break

        # when LC is slightly degenerate, we have slightly fewer variables than equations
        # this leaves some bits to guess, but they seem not to affect the outcome

        # each time a new variable is discovered, nothing before it can depend on it? creating 

        if time_limit and time.time() - start_time > time_limit:
            break

    not_solved = [(x,idx_to_comb[x]) for x in range(num_vars) if x not in var_map]
    for v,c in not_solved:
        print(v,c," - ")
        print("Lower: ", [i for i in range(num_vars) if lower_matrix[i,v] == 1])
        print("Upper: ", [i for i in range(num_vars) if upper_matrix[i,v] == 1])

    output = {}
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
















def alg_attack_online(keystream, state_size, known_bits, attack_data):
    var_map = attack_data['equation times']
    upper_matrix = attack_data['upper matrix']
    lower_matrix = attack_data['lower matrix']
    const_vector = attack_data['constant vector']
    comb_to_idx = attack_data['comb to idx map']
    idx_to_comb = attack_data['idx to comb map']
    num_vars = len(upper_matrix)

    secret_bits = [i for i in range(state_size) if i not in known_bits]
    guess_bits = [i for i in secret_bits if i not in var_map]
    #print("Guessing Bits: ", guess_bits)


    online_vector = np.zeros([num_vars],dtype=np.uint8)
    for v in range(num_vars):
        if v in var_map:
            online_vector[v] = keystream[var_map[v]] ^ const_vector[v]

    solution = lu_solve(
        lower_matrix,
        upper_matrix,
        online_vector
    )

    initial_state = []
    for bit in range(state_size):
        if (bit,) in comb_to_idx:
            initial_state.append(solution[comb_to_idx[(bit,)]])
        else:
            initial_state.append(0)

    return initial_state

# # Online phase:
# total_time = max(var_map.values())

# guess_vars = []
# for i in range(num_vars):
#     if i not in var_map:
#         guess_vars.append(i)

# print("Guessing Bits: ", guess_vars)

# F = FeedbackRegister(random.random(),C)
# print(F)

# seq = np.array([output_filter.eval(state) for state in F.run(total_time+1)],dtype='uint8')
# outs = np.zeros([num_vars],dtype=np.uint8)
# for v in range(num_vars):
#     if v in var_map:
#         outs[v] = (seq[var_map[v]])
    
# for assignment in product((0,1), repeat = len(guess_vars)):
#     for i,v in enumerate(guess_vars):
#         outs[v] = (assignment[i])

#     gf_mat = galois.GF2(lower_matrix) @ galois.GF2(upper_matrix)
#     gf_vec = galois.GF2(outs)

#     sol = [int(x) for x in np.linalg.solve(gf_mat,gf_vec)[:8]]
#     print(C.convert_state(sol))