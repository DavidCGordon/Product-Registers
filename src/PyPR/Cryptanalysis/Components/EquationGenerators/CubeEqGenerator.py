from typing import Iterator, Any

from PyPR.BooleanLogic import BooleanFunction
from PyPR.FeedbackFunctions import FeedbackFunction

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
    include_constant = True,  # whether or not to include an empty tuple
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
    for segment in list_segments:
        # skip the constant segment if needed
        if segment[0] == 0 and not include_constant:
            continue

        for monomial in segment[1]:
            comb_to_idx[monomial] = var_idx
            var_idx += 1

    return comb_to_idx

# small helper function to help pretty-print:
def indent(n)->str:
    """Small helper for pretty printing. Should not be used externally"""
    return ("|   " * n)

# main method
def CubeEqGenerator(
    feedback_fn: FeedbackFunction,
    output_fn: BooleanFunction | list[BooleanFunction], 
    limit: int, 
    var_map: dict[tuple[int,...],int],
    verbose: bool = False, 
    _print_depth: int = 0
) -> Iterator[
    np.ndarray[tuple[int],np.dtype[np.uint8]] |
    list[np.ndarray[tuple[int],np.dtype[np.uint8]]]
]:
    """Generates the Equations for a given feedback function and output function. Where possible,
    this is by far the fastest of the equation generators.

    Like all equation generators, the return is an iterator which iterates through tuples of the
    form `(time, equation, extra_constant)`. Time is an integer which denotes which clock cycle the
    equation represents. The equation can be different types depending on the equation generator,
    but for the CubeEqGenerator, equations are directly generated as np.ndarrays (with `dtype=uint8`,
    and a layout determined by the `var_map` passed in). The extra constant is either 0 or 1. 
    
    This function also allows you to pass in a list of output functions. If you do this,
    the, output which is yielded on each cycle will be a list of tuples, with each corresponding
    to one of the output functions passed in. These tuples each have the same (time, equation, 
    extra_constant) format as if only one output function is passed.\n

    ## How it works:
    The mechanism through which the CubeEqGenerator works is the most complicated of the eq generators,
    and it is the fastest, in places where it is applicable. The core idea is that the coefficient of
    each monomial is determinable by summing over a cube of evaluations with all variable set to zero 
    except those in the cube. However many of these evaluations are reused across cubes. Thus, we keep
    an large array of states which we use to determine the evaluations (simulating the register in parallel,
    with no symbolics). We precompute the indices needed to sum over this array, and then loop through that
    precomputed order to combine evaluations into the coeeficients needed to output, reusing each 
    evaluation many times. If the monomial order is increasing (i.e. each subset of a monomial is computed
    before it), we can reuse previous coefficent to cover large parts of this sum in a principle of
    inclusion/exclusion style sum. This reduces the work needed to compute each monomial  coefficient 
    from 2^(size of monomial) to 2^(size of monomial / 2), effectively square rooting the number of 
    combinations needed for each coefficent (delaying the exponential blowup from cube sums).

    ## Performance:
    The main speed gains come from a couple of things:
    1. There is huge amount of reuse, both in reusing evalutions of the states in different cubes, and
        reusing older coefficients to cover large amounts of the sum. 
    2. The summation indices are the same everytime; after precomputation, this looks like a very tight loop
        of summing values in an array. The storage of everything in dense arrays also helps with cache locality
    3. because the algorithm is so computer-friendly (especially relative to handling symbolics), we can
        numba compile it to speed it up even further. 
    altogether this is probably 2-3 orders of magnitude faster than the other methods

    :param feedback_fn: The function used to update the state on each clock cycle
    :type feedback_fn: FeedbackFunction
    :param output_fn: The function used to derive output from the state
    :type output_fn: BooleanFunction
    :param limit: The number of cycles to generate output for
    :type limit: int
    :param var_map: A map which describes to assign monomials (as sorted tuples) to indices
    :type var_map: dict[tuple[int, ...], int]
    :param verbose: whether or not to print output, defaults to False
    :type verbose: bool, optional
    :param _print_depth: The indentation level to print at (advise not to touch this, it's 
        mostly internal to make printing prettier, and doesnt change much), defaults to 0
    :type _print_depth: int, optional
    :return: An Iterator which yields (time,equation,extra_const) tuples. If a list of output 
        functions are passed as input, the iterator will yield a list of such Tuples on each iteration.
        Otherwise, only one tuple will be yielded each iteration. 
    :rtype: Iterator[
        tuple[int, np.ndarray[tuple[int],np.dtype[np.uint8]], int] | 
        list[tuple[int, np.ndarray[tuple[int],np.dtype[np.uint8]], int]] 
    ]
    """
    # set flags to match outputs shape to input shape:
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    elif isinstance(output_fn,BooleanFunction):
        return_list = False
        output_fn_list = [output_fn]
    else:
        raise TypeError(
            f"output_fn must be a BooleanFunction or list of Boolean functions. " +
            f"Got {type(output_fn)} instead."
        )

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
        print(f"{indent(_print_depth)}Precomputing splits for CubeEqGenerator:")
        precomp_time = time.time()
    
    subcomb_precomputed, subcomb_evals, subcomb_bounds = compute_splits(var_map)

    if verbose:
        print(f"{indent(_print_depth)}Finished computing splits:")
        print(f"{indent(_print_depth)}Time: {time.time()-precomp_time}")
        print(f"{indent(_print_depth)}\n{indent(_print_depth)}Main Loop:")

    eq_vec =  np.zeros([len(output_fn_list),len(var_map)], dtype='uint8')
    for t in range(limit):

        # compiled loop with some DP to reduce number of things to be summed:
        eq_vec = combine_vecs(
            eq_vec, evaluations, 
            subcomb_precomputed, subcomb_evals, subcomb_bounds
        )

        # yield, matching the input format:
        if return_list:
            yield [eq_vec[i] for i in range(len(output_fn_list))]
        else:
            yield eq_vec[0]
           
            
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
def update_states(
    feedback_fn: Any,
    prev_states: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    curr_states: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
) -> tuple[
    np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    np.ndarray[tuple[int,int],np.dtype[np.uint8]]
]:
    """Update all of the states using the feedback function

    :param feedback_fn: The Feedback Function used to update the states
    :type feedback_fn: FeedbackFunction (compiled)
    :param prev_states: An array holding the previous state values
    :type prev_states: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param curr_states: An array to write the updated states to
    :type curr_states: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :return: The updated previous and current state arrays
    :rtype: tuple[
        np.ndarray[tuple[int,int],np.dtype[np.uint8]],
        np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    ]
    """
    prev_states, curr_states = curr_states, prev_states
    for i in range(len(curr_states)):
        curr_states[i] = feedback_fn(prev_states[i])

    return prev_states, curr_states

# update the evaluations array using a specifc output fn:
@numba.njit((u8[:,:])(i64,output_function_type,u8[:,:],u8[:,:]))
def update_evals(
    fn_idx: int, 
    output_fn: Any,
    curr_states: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    evals: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
) -> np.ndarray[tuple[int,int],np.dtype[np.uint8]]:
    """Update the evaluations array from the state array.

    :param fn_idx: the index representing which output function we are using
    :type fn_idx: int
    :param output_fn: the output function used to get the evaluation from the state
    :type output_fn: BooleanFunction (compiled)
    :param curr_states: the array of states to read from
    :type curr_states: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param evals: the array of evaluations to write to
    :type evals: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :return: the updated array of evaluations
    :rtype: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    """
    for idx in range(len(curr_states)):
        evals[idx,fn_idx] = output_fn(curr_states[idx])

    return evals

# use the subcomb data to combine evaluations:
@numba.njit(u8[:,:](u8[:,:],u8[:,:],u64[:],u64[:],u64[:,:]))
def combine_vecs(
    eq_vec: np.ndarray[tuple[int,int],np.dtype[np.uint8]], # equation vector (to read precomputed sums and to write to)
    evals: np.ndarray[tuple[int,int],np.dtype[np.uint8]],  # matrix of function evaluations to compute the sum over
    subcomb_precomputed: np.ndarray[tuple[int],np.dtype[np.int64]], # indices to use for summing precomputed sum
    subcomb_evals: np.ndarray[tuple[int],np.dtype[np.int64]], # indices to sum over for summing new evals
    subcomb_bounds: np.ndarray[tuple[int,int],np.dtype[np.int64]] # start/stop indices to interpret the above arrays properly
    ) -> np.ndarray[tuple[int,int],np.dtype[np.uint8]]:
    """Sum over the evaluations to get the coefficients

    **Proof sketch:** For any sets of evaluations, we have a variant of the principle of
    inclusion-exclusion:
    ```
    (Sum over Set A) xor (Sum over Set B) = (Sum over A union B) xor (Sum over A intersect B)
    ```
     
    For a monomial `M` with variables from `V`, consider the cubes each missing one variable from some 
    set, `S` (i.e. for each `s` in `S`, consider the cube over `(V - s)`) The union of these cubes contains
    every evaluation needed to compute the coefficient of M (except those which contain all of S). 
    Using the identity above, we can calculate the sum over the union as:
    ```
    (Sum over A union B) = (Sum over A) xor (Sum over B) xor (Sum over A intersect B)
    ```

    And the intersection of any two cubes with missing variables (e.g. over `V_1 = V - S_1`, `V2 = V - S_2`, with 
    `S_i` being arbitrary subsets of `S`) is easy to compute: it's the cube over `V_1 intersect V_2` = `V - (S_1 union S_2)`. 
    All of these intersections are cubes which we have already computed, and can be used directly. This gives a P.I.E style
    decomposition of the precomputed sums, and the new evaluations, which is optimized when `|S|` is close to `|V|/2`. 
    
    One branch of the split (subcombs_precomputed) is the P.I.E sum over the already computed coefficients,
    while the the other (subcombs_evals) is the remaining set of evaluations (the combinations having all of S)
    needed to complete the sum

    :param eq_vec: The current equation vector to write to (and to read previously computed sums from)
    :type eq_vec: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param evals: The array of evaluations to read from when summing
    :type evals: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param subcomb_precomputed: The subcombinations for precomputed sums (reading from the equation vector)
    :type subcomb_precomputed: np.ndarray[tuple[int],np.dtype[np.int64]]
    :param subcomb_evals: The subcombinations for new evaluations (reading from the evals array)
    :type subcomb_evals: np.ndarray[tuple[int],np.dtype[np.int64]]
    :param subcomb_bounds: An array which contains the bounds for the subcomb arrays, which helps to
        know which indices correspond to each coefficient
    :type subcomb_bounds: np.ndarray[tuple[int,int],np.dtype[np.int64]]
    :return: The updated / filled out equation vector
    :rtype: np.ndarray[tuple[int, int],np.dtype[np.uint8]]
    """
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

def compute_splits(
    var_map:dict[tuple[int,...],int]
) -> tuple[
    np.ndarray[tuple[int],np.dtype[np.int64]],
    np.ndarray[tuple[int],np.dtype[np.int64]],
    np.ndarray[tuple[int,int],np.dtype[np.int64]]
]:
    """precompute the subcomb splits and iteration order

    :param var_map: the mapping from combinations to indices
    :type var_map: dict[tuple[int,...], int]
    :return: the two subcombination splits (precomputed and evals)
        and the bounds to interpret them
    :rtype: tuple[ndarray, ndarray, ndarray]
    """
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
