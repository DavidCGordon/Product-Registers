import numpy as np
import time

import numba
u8 = numba.types.uint8
i64 = numba.types.int64
u64 = numba.types.uint64
feedback_function_type = numba.types.FunctionType(u8[:](u8[:]))
output_function_type = numba.types.FunctionType(u8(u8[:]))

def get_var_map(
    feedback_fn,
    monomial_profile,
    variable_blocks,
    include_variables = True, # include all base variables
    complete_subsets = False, # ensure variable map is closed under subsets
    lexicographic = True      # sort monomials lexicographically
):
    mons_by_len = {}

    for selectors in monomial_profile.get_monomials(complete_subsets=complete_subsets):
        # construct monomial
        mon = []
        for block_idx in range(len(selectors)):
            for bit_idx in selectors[block_idx]:
                mon.append(variable_blocks[block_idx][bit_idx])
        
        # append to the appropriate list
        mon = tuple(sorted(mon))
        if len(mon) in mons_by_len:
            mons_by_len[len(mon)].append(mon)
        else:
            mons_by_len[len(mon)] = [mon]

    # replace length 1 segment if necessary:
    if include_variables:
        mons_by_len[1] = [(i,) for i in range(feedback_fn.size)] 

    # cosmetic changes to order:
    list_segments = []
    for length, mon_list in mons_by_len.items():
        if lexicographic:
            sorted_mons = [m for m in mon_list]
            for i in range(length):
                sorted_mons = sorted(sorted_mons, key=lambda x: x[i])
            list_segments.append((length,sorted_mons))
        else:
            list_segments.append((length,mon_list))
    list_segments = sorted(list_segments, key = lambda x: x[0])

    # merging list segments into output maps
    var_idx = 0
    comb_to_idx = {}
    idx_to_comb = {}
    for segment in list_segments:
        for monomial in segment[1]:
            comb_to_idx[monomial] = var_idx
            idx_to_comb[var_idx] = monomial
            var_idx += 1

    return comb_to_idx #, idx_to_comb

# small helper function to help pretty-print:
def indent(n):
    return ("|   " * n)

# main method
def CubeEqGenerator(
    feedback_fn, output_fn, limit, var_map, verbose=False, print_depth = 0
):
    # set flags to match outputs shape to input shape:
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    else:
        return_list = False
        output_fn_list = [output_fn]


    # move constants to the end of the maps:
    var_map = {k:v for k,v in var_map.items() if k != tuple()}
    var_map[tuple()] = len(var_map)

    # variable inits
    num_bits = len(feedback_fn)
    feedback_fn = feedback_fn.compile()
    output_fn_list = [fn.compile() for fn in output_fn_list]

    # create arrays for the state values
    prev_states = np.zeros([len(var_map),num_bits], dtype='uint8')
    curr_states = np.zeros([len(var_map),num_bits], dtype='uint8')
    for comb,idx in var_map.items():
        for v in comb:
            curr_states[idx,v] = 1

    # create array to hold evaluations for reuse:
    evaluations = np.zeros([len(var_map),len(output_fn_list)], dtype='uint8')
    for fn_idx in range(len(output_fn_list)):
        for comb,idx in var_map.items():
            evaluations[idx,fn_idx] = output_fn_list[fn_idx](curr_states[idx])

    if verbose:
        print(f"{indent(print_depth)}Precomputing splits for CubeEqGenerator:")
        precomp_time = time.time()
    
    subcomb_precomputed, subcomb_evals, subcomb_bounds = compute_splits(var_map)

    if verbose:
        print(f"{indent(print_depth)}Finished computing splits:")
        print(f"{indent(print_depth)}Time: {time.time()-precomp_time}")
        print(f"{indent(print_depth)}\n{indent(print_depth)}Main Loop:")

    eq_vec =  np.zeros([len(output_fn_list),len(var_map)], dtype='uint8')
    for t in range(limit):

        # compiled loop with some DP to reduce number of things to be summed:
        eq_vec = combine_vecs(
            eq_vec, evaluations, 
            subcomb_precomputed, subcomb_evals, subcomb_bounds
        )

        # yield, matching the input format:
        if return_list:
            yield [
                (t, eq_vec[i,:-1], eq_vec[i, -1])
                for i in range(len(output_fn_list))
            ]
        else:
            yield (t, eq_vec[0,:-1], eq_vec[0, -1])
           
            
        # update the current states and evaluations:
        prev_states,curr_states = update_states(
            feedback_fn,prev_states,curr_states
        )

        for fn_idx, output_fn in enumerate(output_fn_list):
            evaluations = update_evals(
                fn_idx, output_fn, curr_states, evaluations
            )

# update prev and current state arrays using the feedback fn:
@numba.njit(numba.types.Tuple((u8[:,:],u8[:,:]))(feedback_function_type,u8[:,:],u8[:,:]))
def update_states(feedback_fn,prev_states,curr_states):
    prev_states, curr_states = curr_states, prev_states
    for i in range(len(curr_states)):
        curr_states[i] = feedback_fn(prev_states[i])

    return prev_states, curr_states

# update the evaluations array using a specifc output fn:
@numba.njit((u8[:,:])(i64,output_function_type,u8[:,:],u8[:,:]))
def update_evals(fn_idx, output_fn, curr_states, evals):
    for idx in range(len(curr_states)):
        evals[idx,fn_idx] = output_fn(curr_states[idx])

    return evals

# use the subcomb data to combine evaluations:
@numba.njit(u8[:,:](u8[:,:],u8[:,:],u64[:],u64[:],u64[:,:]))
def combine_vecs(
    eq_vec, # equation vector (to read precomputed sums and to write to)
    evals,  # matrix of function evaluations to compute the sum over
    subcomb_precomputed, # indices to use for summing precomputed sum
    subcomb_evals, # indices to sum over for summing new evals
    subcomb_bounds # start/stop indices to interpret the above arrays properly
    ):

    # Note that for any sets of evaluations, we have a variant of the principle of inclusion-exclusion: 
    #   (Sum over Set A) + (Sum over Set B) = (Sum over A union B) + (Sum over A intersect B)

    # Consider for some set of variables, S, consider the cubes with each variable missing: 
    #   -  i.e. the cube over (V - s) for each s in S
    #
    # the union of all of these sets is every evaluation except those which contain all of S
    # the intersection of any two cubes (e.g. over V_1 = V - S_1, V2 = V - S_2, with S_i being
    # arbitrary subsets of S) is the cube over V_1 intersect V_2 = V - (S_1 union S_2). All of
    # these intersections are cubes which we have already computed, and can be used directly.
    #
    # one branch of the split is the P.I.E sum over every intersection,
    # the other is the remaining set of evaluations (having all of S)

    for term_idx in range(eq_vec.shape[1]):
        for fn_idx in range(eq_vec.shape[0]):
            eq_vec[fn_idx, term_idx] = 0

        # first half of the split uses precomputed sums (P.I.E sum)
        for i in range(subcomb_bounds[term_idx,0],subcomb_bounds[term_idx+1,0]):
            for fn_idx in range(eq_vec.shape[0]):
                eq_vec[fn_idx, term_idx] ^= eq_vec[fn_idx,subcomb_precomputed[i]]
        
        # second half of split add the remaining subcomb evaluations
        for i in range(subcomb_bounds[term_idx,1],subcomb_bounds[term_idx+1,1]):
            for fn_idx in range(eq_vec.shape[0]):
                eq_vec[fn_idx, term_idx] ^= evals[subcomb_evals[i],fn_idx]

    return eq_vec

def compute_splits(var_map):
    subcomb_precomputed = []
    subcomb_evals = []

    # use a random split:
    for comb in var_map.keys():
        split_1 = tuple(sorted(np.random.choice(comb,len(comb)//2, replace=False)))
        split_2 = tuple(sorted([x for x in comb if x not in split_1]))

        subcomb_precomputed.append(split_1)
        subcomb_evals.append(split_2)

    # convert splits into index data:
    eval_indices = np.zeros([sum(2**len(x) for x in subcomb_evals)], dtype='uint64')
    precomputed_indices = np.zeros([sum(2**len(x) for x in subcomb_precomputed)], dtype='uint64')
    output_bounds = np.zeros([len(var_map)+1,2],dtype='uint64')
    for term_idx,(precompute_vars, eval_vars) in enumerate(
        zip(subcomb_precomputed, subcomb_evals)
    ):
        
        # precompute indices formed by holding eval variables fixed and
        # summing over a cube of the precompute variables (except all variables)
        for i in range(2**len(precompute_vars)-1):
            subcomb = eval_vars
            subcomb += tuple([
                precompute_vars[idx] for idx in range(len(precompute_vars))
                if (i & (1 << idx))]
            )
        
            subcomb = tuple(sorted(subcomb))
            precomputed_indices[output_bounds[term_idx,0] + np.uint64(i)] = var_map[subcomb]
        output_bounds[term_idx+1,0] = output_bounds[term_idx,0] + (2**len(precompute_vars)-1)

        # similarly, eval indices formed by holding precompute variables fixed and
        # summing over a cube of the eval variables (full cube this time)
        for i in range(2**len(eval_vars)):
            subcomb = precompute_vars
            subcomb += tuple([
                eval_vars[idx] for idx in range(len(eval_vars))
                if (i & (1 << idx))]
            )
        
            subcomb = tuple(sorted(subcomb))
            eval_indices[output_bounds[term_idx,1] + np.uint64(i)] = var_map[subcomb]
        output_bounds[term_idx+1,1] = output_bounds[term_idx,1] + (2**len(eval_vars))

    return precomputed_indices, eval_indices, output_bounds
