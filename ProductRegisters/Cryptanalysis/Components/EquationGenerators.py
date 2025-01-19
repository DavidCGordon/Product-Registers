import numba.typed
import numpy as np
import time

import numba
u8 = numba.types.uint8
i64 = numba.types.int64
u64 = numba.types.uint64
feedback_function_type = numba.types.FunctionType(u8[:](u8[:]))
output_function_type = numba.types.FunctionType(u8(u8[:]))

# small helper function to help pretty-print:
def indent(n):
    return ("|   " * n)


def EqGenerator(feedback_fn, output_fn, limit):
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    else:
        return_list = False
        output_fn_list = [output_fn]

    for t, equations in enumerate(
        feedback_fn.anf_iterator(
            limit, 
            bits = list(sorted(set().union(*(
                fn.idxs_used() for fn in output_fn_list
            ))))
        )
    ):
        if not return_list:
            # because constant is included in ANF, no extra const needed
            yield (t, output_fn.compose(equations).translate_ANF(), 0)
        else:
            yield [
                (t, output_fn.compose(equations).translate_ANF(), 0)
                for output_fn in output_fn_list
            ]


















# reasonable default:
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


    # Note that for any sets of evaluationsd, we have a variant of the principle of inclusion-exclusion: 
    #   (Sum over Set A) + (Sum over Set B) = (Sum over A union B) + (Sum over A intersect B)

    # Consider for some set of variables, S, the sets with each variable missing: (V - s) for each s in S
    # the union of all of these sets is every evaluation except those which contain all of S
    # one branch of the split is the P.I.E sum over every intersection,
    # the other is the last set of evaluation (having all of S)

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
    #subcomb_precomputed, subcomb_evals = split_fn(var_map, output_map)
    subcomb_precomputed = []
    subcomb_evals = []

    # use a random split:
    for comb,idx in var_map.items():
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























def block_splits(blocks):
    def inner_function(var_map,output_map):
        subcomb_precomputed = []
        subcomb_evals = []

        # if blocks are known we can guarantee a good split
        if blocks != None:
            for comb, idx in output_map.items():
                seen_blocks = []
                data_1 = []
                data_2 = []

                # first pass collect minimum set in data 1
                for bit in comb:
                    for block_idx in range(len(blocks)):
                        if (block_idx not in seen_blocks) and (bit in blocks[block_idx]):
                            seen_blocks.append(block_idx)
                            data_1.append(bit)
                
                # second pass, assign bits to data 2 until they are even
                for bit in comb:
                    if bit not in data_1:
                        if len(data_2) < len(comb)//2:
                            data_2.append(bit)
                        else:
                            data_1.append(bit)

                # insert into lists
                subcomb_precomputed.append(tuple(sorted(data_1)))
                subcomb_evals.append(tuple(sorted(data_2)))
        return subcomb_precomputed, subcomb_evals
    return inner_function

def random_split(num_trials):
    def inner_function(var_map,output_map):
        s = 0
        subcomb_precomputed = []
        subcomb_evals = []

        for comb,idx in output_map.items():
            # try to randomly guess a good split:
            success = False
            for i in range(num_trials):
                data_1 = tuple(sorted(np.random.choice(comb,(len(comb)+1)//2,replace=False)))
                if data_1 in output_map and output_map[data_1] < idx:
                    success = True
                    break

            if success:
                s += 1
                # if you got a good split add all remaining bits to data 2
                data_2 = tuple(sorted([x for x in comb if x not in data_1]))
            else:
                # otherwise the whole comb as to be data 1
                data_1 = comb
                data_2 = tuple()

            subcomb_precomputed.append(data_1)
            subcomb_evals.append(data_2)
        print(s / len(output_map))
        return subcomb_precomputed, subcomb_evals
    return inner_function

def compute_splits_restricted(var_map,output_map, split_fn = random_split(1000)):
    subcomb_precomputed, subcomb_evals = split_fn(var_map, output_map)

    output_1_data = np.zeros([sum(2**len(x) for x in subcomb_evals)], dtype='uint64')
    output_2_data = np.zeros([sum(2**len(x) for x in subcomb_precomputed)], dtype='uint64')
    output_bounds = np.zeros([len(output_map)+1,2],dtype='uint64')
    for term_idx,(data_1,data_2) in enumerate(zip(subcomb_precomputed, subcomb_evals)):
        
        # output 1 formed by holding data 1 fixed:
        for i in range(2**len(data_2)-1):
            subcomb = tuple([data_2[idx] for idx in range(len(data_2)) if (i & (1 << idx))])
            subcomb += data_1 # <= data 1 is fixed
            subcomb = tuple(sorted(subcomb))
            output_1_data[output_bounds[term_idx,0] + np.uint64(i)] = output_map[subcomb]
        output_bounds[term_idx+1,0] = output_bounds[term_idx,0] + (2**len(data_2)-1)

        # output 2 formed by holding data 2 fixed:
        for i in range(2**len(data_1)):
            subcomb = tuple([data_1[idx] for idx in range(len(data_1)) if (i & (1 << idx))])
            subcomb += data_2 # <= data 2 is fixed
            subcomb = tuple(sorted(subcomb))
            output_2_data[output_bounds[term_idx,1] + np.uint64(i)] = var_map[subcomb]
        output_bounds[term_idx+1,1] = output_bounds[term_idx,1] + 2**len(data_1)

    return output_1_data, output_2_data, output_bounds




def CubeEqGenerator_restricted(
    feedback_fn, output_fn, limit, var_map, output_map,
    split_function = random_split(num_trials = 1000)
):
    # set flags to match outputs shape to input shape:
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    else:
        return_list = False
        output_fn_list = [output_fn]


    # move constants to the end of the maps:
    output_map = {k:v for k,v in output_map.items() if k != tuple()}
    var_map = {k:v for k,v in var_map.items() if k != tuple()}
    output_map[tuple()] = len(output_map)
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


    subcomb_precomputed, subcomb_evals, subcomb_bounds = compute_splits_restricted(
        var_map, output_map, split_fn=split_function
    )

    eq_vec =  np.zeros([len(output_fn_list),len(output_map)], dtype='uint8')
    for t in range(limit):

        eq_vec = combine_vecs(eq_vec, evaluations, subcomb_precomputed, subcomb_evals, subcomb_bounds)
        #eq_vec = combine_vecs(eq_vec, evaluations, subcomb_data, subcomb_bounds)

        if not return_list:
            yield (t, eq_vec[0,:-1], eq_vec[0, -1])
        else:
            yield [
                (t, eq_vec[i,:-1], eq_vec[i, -1])
                for i in range(len(output_fn_list))
            ]
            
        # update the current states and evaluations:
        prev_states,curr_states = update_states(
            feedback_fn,prev_states,curr_states
        )

        for fn_idx, output_fn in enumerate(output_fn_list):
            evaluations = update_evals(
                fn_idx, output_fn, curr_states, evaluations
            )


    
def compute_splits_restricted(var_map,output_map,split_fn=random_split(1000)):
    subcomb_precomputed, subcomb_evals = split_fn(var_map, output_map)

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
        output_bounds[term_idx+1,1] = output_bounds[term_idx,1] + (2**len(eval_vars)-1)

    return precomputed_indices, eval_indices, output_bounds