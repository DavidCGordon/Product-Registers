import time

import numba
import numpy as np
from itertools import combinations,chain,cycle,product,tee

from PyPR import FeedbackRegister
from PyPR.BooleanLogic import BooleanFunction
from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile, TermSet

def indent(n):
    return ("|   " * n)

# Cube attacks need to tweak/query the actual register:
# because of this, we need pass functions to the attack
# which allow it to interface with the target system
def access_fns(register, output_fn, tweakable_bits, init_rounds=100, keystream_len=None):
    # default keystream len:
    if keystream_len == None:
        keystream_len = max(100,2*register.size)

    # compile as needed:
    if not hasattr(register.fn,'_compiled'):
        register.fn.compile()
    if not hasattr(output_fn,'_compiled'):
        output_fn.compile()

    for i in range(init_rounds):
        register.clock()

    keystream = [
        output_fn._compiled(state._state) 
        for state in register.run(keystream_len)
    ]

    register.reset()

    # given a full state, simulate that state to get the bit
    def sim_fn(state):
        register.set_state(state)
        for i in range(init_rounds):
            register.clock()

        # generate keystream as normal:
        keystream = [
            output_fn._compiled(state._state) 
            for state in register.run(keystream_len)
        ]

        register.reset()
        return np.array(keystream, dtype = np.uint8)

    # access a "real" model, with potentially limited I/O access:
    # input may contain None to use the underlying secret key,
    # output may contain None to signify impossible values.
    def access_fn(state):
        register.reset()
        
        # write only to tweakable bits
        for bit in tweakable_bits:
            if state[bit] != None:
                register[bit] = state[bit]
        
        # initialization rounds:
        for i in range(init_rounds):
            register.clock()

        # generate keystream as normal:
        keystream = [
            output_fn._compiled(state._state) 
            for state in register.run(keystream_len)
        ]

        register.reset()
        return np.array(keystream, dtype = np.uint8)
    
    # test a state to see if the keystream is correct
    def test_fn(state):
        register.set_state(state)
        for i in range(init_rounds):
            register.clock()

        test_keystream = [
            output_fn._compiled(state._state) 
            for state in register.run(keystream_len)
        ]

        register.reset()
        return test_keystream == keystream

    return access_fn,sim_fn,test_fn




def insert_equation(
    lower_matrix,upper_matrix, const_vec, cube_map, # structures we modify
    maxterm, time, equation, const                  # data we update with
    ):

    linearly_independent = False
    modification_vector  = np.zeros_like(equation)
    for bit in range(len(equation)):
        if equation[bit] == 1:
            modification_vector[bit] = 1
            if bit in cube_map:
                equation ^= upper_matrix[bit]
            else:
                linearly_independent = True
                cube_map[bit] = (maxterm, time)
                const_vec[bit] = const
                lower_matrix[bit] = modification_vector
                upper_matrix[bit] = equation
                break
    return linearly_independent




# def cube_attack_offline(
#     feedback_fn, sim_fn, tweakable_vars, 
#     time_limit = None, num_tests = 20, verbose = False
#     ):

#     tweakable_vars = set(tweakable_vars)
#     start_time = time.time()
    
#     cube_map = {}
#     lower_matrix = np.eye(feedback_fn.size,dtype=np.uint8)
#     upper_matrix = np.eye(feedback_fn.size,dtype=np.uint8)
#     const_vec = np.zeros([feedback_fn.size,1],dtype=np.uint8)

#     failure_count = 0
#     already_seen = set()
#     cube_variables = set([list(tweakable_vars)[0]])
#     while True:
#         # check to make sure cubes are only checked once
#         cube = tuple(sorted(list(cube_variables)))
#         if cube in already_seen:
#             failure_count += 1
#             added_element = np.random.choice(list(tweakable_vars - cube_variables))
#             removed_element = np.random.choice(cube)
#             cube_variables.add(added_element)
#             cube_variables.remove(removed_element)

#             if failure_count > 100:
#                 if verbose:
#                     print("too many repeated cubes in random walk!")
#                 break
#             continue

#         failure_count = 0
#         already_seen.add(cube)
#         print("Cube Candidate: ", cube)

#         # get cube information:
#         equations, constants = determine_equations(sim_fn,cube,feedback_fn.size)
#         nonlinear_mask = get_nonlinear_mask(sim_fn,cube,feedback_fn.size,num_tests)
#         constant_mask = get_constant_mask(sim_fn,cube,feedback_fn.size,num_tests)

#         # counts for bookkeeping/printing:
#         useful_count = 0
#         constant_count = 0
#         nonlinear_count = 0
#         dependent_count = 0

#         for t in range(len(nonlinear_mask)):
#             # filter constant / nonlinear superpoly's
#             if constant_mask[t]:
#                 constant_count += 1
#                 continue
#             elif nonlinear_mask[t]:
#                 nonlinear_count += 1
#                 continue

#             # attempt to insert equation, and determ
#             linearly_independent = insert_equation(
#                 lower_matrix, upper_matrix, const_vec, cube_map,
#                 cube, t, equations[t], constants[t]
#             )

#             # determine whether the insert was successful
#             if linearly_independent:
#                 useful_count += 1
#             if not linearly_independent:
#                 dependent_count += 1
            
#         # add or remove elements randomly as needed:
#         #  - move up when there are any nonlinear terms
#         #  - move down when there are all constant terms
#         #  - otherwise just swap a random element
#         added_element = np.random.choice(list(tweakable_vars - cube_variables))
#         removed_element = np.random.choice(cube)
#         #print(added_element,removed_element, not np.all(constant_mask), not np.any(nonlinear_mask))
#         if not np.all(constant_mask):
#             cube_variables.add(added_element)
#         if not np.any(nonlinear_mask):
#             cube_variables.remove(removed_element)
#         #print("New: ", cube_variables)
#         # print to keep information up to date:
#         if verbose: 
#             print(
#                 f" - Useful: {useful_count} -- " +
#                 f"Constant: {constant_count} -- " +
#                 f"Nonlinear: {nonlinear_count} -- " +
#                 f"Dependent: {dependent_count}",
#             )

                
#         # this breaks out of the loop indexing the keystream by time
#         # the check at the top of this section breaks the individual cube loop
#         if all([(bit in cube_map) for bit in range(feedback_fn.size)]):
#             if verbose: print("all variables solved!")
#             break
#         if time_limit and time.time() - start_time > time_limit:
#             if verbose: print("time limit reached!")
#             break 
    
#     num_queries = 0
#     distinct_cubes = set()
#     for (cube, t) in cube_map.values():
#         if cube not in distinct_cubes:
#             num_queries += 2**len(cube)
#             distinct_cubes.add(cube)

#     if verbose:  
#         print("Number of cubes tested: ", len(already_seen))
#         print("Number of cubes found: ", len(cube_map))
#         print("Num Queries: ", num_queries)

#     output = {}
#     output['cubes'] = cube_map
#     output['lower matrix'] = lower_matrix
#     output['upper matrix'] = upper_matrix
#     output['constant vector'] = const_vec
#     return output











def cmpr_cube_summary(cmpr_fn, output_fn,tweakable_vars, analyze_sources = False):
    print('Beginning Summary:')

    tweakable_set = set(tweakable_vars)
    tweakable_counts = [len(set(block) & tweakable_set) for block in cmpr_fn.blocks]

    print('Computing Monomial Profile')
    output_anf = output_fn.translate_ANF()
    monomial_profiles = cmpr_fn.monomial_profiles()
    output_profile = output_anf.remap_constants([
        (0, MonomialProfile.logical_zero()),
        (1, MonomialProfile.logical_one())
    ]).eval_ANF(monomial_profiles)

    cube_candidates = sorted(
        output_profile.get_cube_candidates(),
        key = (lambda x: sum(x[0].counts.values()))
    )

    print('Analyzing Cube Candidates:')
    for cube_profile, target_block, num_cubes, cube_failure_prob in cube_candidates:
        # calculate actual number of tweakable cubes:
        tweakable_cube_count = 1

        for block_id in range(len(cmpr_fn.blocks)-1,-1,-1):
            if block_id in cube_profile.counts:
                # this is just product(choose(tweakable_count, cube_count))
                for i in range(cube_profile.counts[block_id]):
                    tweakable_cube_count *= (
                        (tweakable_counts[block_id] - i) /
                        (cube_profile.counts[block_id] - i)
                    )
                
        # round float to get integer approximation for number of actual cubes
        # rounding errors should not be too significant here, as only general size is needed.
        tweakable_cube_count = round(tweakable_cube_count)
        if tweakable_cube_count == 0:
            print('Cube Profile: ', cube_profile, "- Not possible with current tweakable set.")
            continue

        
        if analyze_sources:
            # reconstruct the parent term
            parent_term = TermSet(
                {k:v for k,v in cube_profile.totals.items()},
                {k:v for k,v in cube_profile.counts.items()})
            parent_term.counts[target_block] += 1

            sources = []
            for term in output_anf.args:
                term_profile = term.remap_constants([
                    (0, MonomialProfile.logical_zero()),
                    (1, MonomialProfile.logical_one())
                ]).eval_ANF(monomial_profiles)

                # check if parent_term == output_term
                for output_monomial in term_profile.terms:
                    if ((output_monomial.totals == parent_term.totals) and
                        (output_monomial.counts == parent_term.counts)
                    ):
                        sources.append(term)

        print('Cube Profile: ', cube_profile)
        print('Target Block: ', target_block, '- Target Block Size:', len(cmpr_fn.blocks[target_block]))
        print('Number of Cube Candidates (before restriction): ', num_cubes)
        print('Number of Cube Candidates (restricted to tweakable bits): ', tweakable_cube_count)
        if analyze_sources:
            print('Source terms: ')
            for source_term in sources:
                print(f" - {source_term.dense_str()}")
                for var in source_term.args:
                    var_str = str(monomial_profiles[var.index])
                    if len(var_str) >= 80:
                        var_str = var_str[:80]
                        var_str += f'... ({len(monomial_profiles[var.index].terms)} terms)'
                    print(f"   - {var.index}: {var_str}")
        print('\n')
        
    print("Summary Finished!")




# lazy product implementation for faster skipping of unusable sets :)
# attribution: https://discuss.python.org/t/a-product-function-which-supports-large-infinite-iterables/5753
def iproduct(*iterables):
    N = len(iterables)
    saved = [[] for _ in range(N)]  # All the items that we have seen of each iterable.
    exhausted = set()               # The set of indices of iterables that have been exhausted.

    idx = -1
    while True:
        idx = (idx+1) % N
        if idx in exhausted:  # dont increment exhausted iterators
            continue

        try:
            item = next(iterables[idx])
            # yield to products involving the new item:
            yield from product(*saved[:idx], [item], *saved[idx+1:]) 
            saved[idx].append(item)

        # Product is empty or all iterables exhausted.
        except StopIteration:
            exhausted.add(idx)
            if not saved[idx] or len(exhausted) == N:  
                return
    yield ()  # There are no iterables.

def cmpr_cube_attack_offline(
    cmpr_fn, output_fn, sim_fn, tweakable_vars,
    time_limit = None, num_tests = 20, verbose = False, print_depth=0
    ):
    print_skipped_cubes = False
    if verbose:
        print(f"{indent(print_depth)}Starting offline phase (Naive Algebraic Attack):")

    start_time = time.time()
    
    cube_map = {}
    lower_matrix = np.eye(cmpr_fn.size,dtype=np.uint8)
    upper_matrix = np.eye(cmpr_fn.size,dtype=np.uint8)
    const_vec = np.zeros([cmpr_fn.size,1],dtype=np.uint8)

    # break up tweakable variables by block and compute cube candidates:
    tweakable_set = set(tweakable_vars)
    tweakable_blocks = [set(block) & tweakable_set for block in cmpr_fn.blocks]

    if verbose:
        print(f"{indent(print_depth+1)}using monomial profile optimization: True")
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating larger monomial profile:")
        mp_time = time.time()

    monomial_profile = output_fn.translate_ANF().remap_constants([
        (0, MonomialProfile.logical_zero()),
        (1, MonomialProfile.logical_one())
    ]).eval_ANF(cmpr_fn.monomial_profiles())

    if verbose:
        print(f"{indent(print_depth+1)}Monomial profile computed:")
        print(f"{indent(print_depth+1)}Time: {time.time() - mp_time} s")
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating variable_map:")
        cube_cand_time = time.time()

    cube_candidates = sorted(
        monomial_profile.get_cube_candidates(),
        key = (lambda x: sum(x[0].counts.values()))
    )

    if verbose:
        print(f"{indent(print_depth+1)}Cube candidates computed:")
        print(f"{indent(print_depth+1)}Time: {time.time() - cube_cand_time} s")
        print(f"{indent(print_depth+1)}")
        print(f"{indent(print_depth+1)}Identifying Cubes:")

    # Maxterm search
    maxterm_count = 0
    for cube_profile, target_block, num_cubes in cube_candidates:
        if verbose:
           profile_prefix = f"{indent(print_depth+2)}Cube Candidate Profile: {cube_profile}"

        # check to see if the block is saturated:
        block_already_saturated = True
        for t in target_block:
            for bit in cmpr_fn.blocks[t]:
                if not (bit in cube_map):
                    block_already_saturated = False
                    break

        # create the iterators and calculate some statistics:
        tweakable_cube_count = 1
        variable_iterators = []

        # only compute tweakable bits for blocks which are not saturated
        if not block_already_saturated:
            loop_nums = []
            for block_id in range(len(cmpr_fn.blocks)-1,-1,-1):
                if block_id in cube_profile.counts:
                    
                    variable_iterators.append(combinations(tweakable_blocks[block_id],cube_profile.counts[block_id]))

                    num_loops = 1
                    for i in range(cube_profile.counts[block_id]):
                        num_loops *= (
                            (len(tweakable_blocks[block_id])-i)/
                            (cube_profile.counts[block_id]-i)
                        )
                    loop_nums.append(num_loops)
                    tweakable_cube_count *= num_loops

        # round float to get integer approximation for number of actual cubes
        tweakable_cube_count = round(tweakable_cube_count)
        variable_iterators= [x[1] for x in sorted(zip(loop_nums,variable_iterators), key = lambda x:x[0])]

        # output message for empty cube profiles:
        if tweakable_cube_count == 0:
            if print_skipped_cubes and verbose: 
                print(f'{profile_prefix}: skipped (not possible with current tweakable bits)')
            continue

        # output message for cube profiles we won't use but could:
        if block_already_saturated:
            if print_skipped_cubes and verbose: 
                print(f'{profile_prefix}: skipped (target block already saturated)')
            continue

        # test the individual cubes/maxterms:
        already_printed = False
        for cube_idx, var_selections in enumerate(iproduct(*variable_iterators)):
            # break out of the specific cube candidate loop if needed
            block_already_saturated = True
            for t in target_block:
                for bit in cmpr_fn.blocks[t]:
                    if not (bit in cube_map):
                        block_already_saturated = False
                        break

            if block_already_saturated:
                break
            if time_limit and time.time() - start_time > time_limit:
                break 
            
            # print only inside the loop to make sure there are actual cubes
            # depending on the tweakable set, this iterator may be empty
            if not already_printed:
                already_printed = True
                if verbose:
                    print(f'{profile_prefix}:')
                    print(f'{indent(print_depth+2)}Target Block: {target_block} - Target Block Size: {sum(len(cmpr_fn.blocks[t]) for t in target_block)}')
                    print(f'{indent(print_depth+2)}Number of Cube Candidates (before restriction): {num_cubes}')
                    print(f'{indent(print_depth+2)}Number of Cube Candidates (restricted to tweakable bits): {tweakable_cube_count}')

            maxterm_count += 1
            maxterm = tuple(chain(*var_selections))
            if verbose:
                print(f'\r{indent(print_depth+3)}Cube {cube_idx+1}/{tweakable_cube_count}: {maxterm}',end='')

            # counts for bookkeeping/printing:
            useful_count = 0
            constant_count = 0
            nonlinear_count = 0
            dependent_count = 0

            # cube information:
            equations, constants = determine_equations(sim_fn,maxterm,cmpr_fn.size)
            nonlinear_mask = get_nonlinear_mask(sim_fn,maxterm,cmpr_fn.size,num_tests)
            constant_mask = get_constant_mask(sim_fn,maxterm,cmpr_fn.size,num_tests)
            
            for t in range(len(constant_mask)):
                # filter constant / nonlinear superpoly's
                if constant_mask[t]:
                    constant_count += 1
                    continue
                # elif nonlinear_mask[t]:
                #     nonlinear_count += 1
                #     continue

                # attempt to insert equation, and determ
                linearly_independent = insert_equation(
                    lower_matrix, upper_matrix, const_vec, cube_map,
                    maxterm, t, equations[t], constants[t]
                )

                # determine whether the insert was successful
                if linearly_independent:
                    useful_count += 1
                if not linearly_independent:
                    dependent_count += 1
                
                # print to keep information up to date:
                # if verbose: 
                    # print(
                        
                    #     f"Useful: {useful_count} -- " +
                    #     f"Constant: {constant_count} -- " +
                    #     f"Nonlinear: {nonlinear_count} -- " +
                    #     f"Dependent: {dependent_count}",
                    # end='')

                # this breaks out of the loop indexing the keystream by time
                # the check at the top of this section breaks the individual cube loop
                block_already_saturated = True
                for t in target_block:
                    for bit in cmpr_fn.blocks[t]:
                        if not (bit in cube_map):
                            block_already_saturated = False
                            break

                if block_already_saturated:
                    if verbose: print(f'\n{indent(print_depth+2)}Target Block Saturated!')
                    break
                if time_limit and time.time() - start_time > time_limit:
                    break

            # flush the print statements with a newline
            # if verbose: print(f'{indent(print_depth+2)}')

        # This check breaks out of the monomial profile loop
        # no block saturated check because those are profile-specific
        if time_limit and time.time() - start_time > time_limit:
            break  


    num_queries = 0
    distinct_cubes = set()
    for (cube, t) in cube_map.values():
        if cube not in distinct_cubes:
            num_queries += 2**len(cube)
            distinct_cubes.add(cube)

    if verbose:
        print(f'{indent(print_depth+1)}Finished equation generation: ')
        print(f'{indent(print_depth+1)}Time: {time.time() - cube_cand_time} s')
        print(f'{indent(print_depth+1)}')
        print(f'{indent(print_depth+1)}Number of equations found: {len(cube_map)}')
        print(f'{indent(print_depth+1)}Number of cube tested: {maxterm_count}')
        print(f'{indent(print_depth+1)}Number of queries in attack: {num_queries}')
        print(f'Offline phase complete -- Total time: ', time.time() - start_time)
    
    output = {}
    output['cubes'] = cube_map
    output['lower matrix'] = lower_matrix
    output['upper matrix'] = upper_matrix
    output['constant vector'] = const_vec
    return output




u8 = numba.types.uint8
@numba.njit(u8[:,:](u8[:,:],u8[:,:],u8[:,:]))
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




def cube_attack_online(access_fn, test_fn, state_size, known_bits, cube_data, verbose = False, print_depth=0):
    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Cube Attack):")

    cube_map = cube_data['cubes']
    lower_matrix = cube_data['lower matrix']
    upper_matrix = cube_data['upper matrix']
    total_matrix = (lower_matrix @ upper_matrix) % 2
    consts = cube_data ['constant vector']

    # create copies to prevent known-variable reduction from deleting important information
    secret_bits = [i for i in range(state_size) if i not in known_bits]
    guess_bits = [i for i in secret_bits if i not in cube_map]
    if verbose: 
        print(f"{indent(print_depth+1)}Guessing Bits: {tuple(sorted(guess_bits))}")

    # if no cube bits, then cube attack is slower than brute force:
    # exit immediately
    if not cube_map:
        raise ValueError(
            'No cubes given; consider either providing ' + 
            'cubes for the attack or a brute force approach'
        )

    # using known vars and equations found, re-do the LU factorization:
    # this is slightly slower, but saves cube evaluations
    # for big cubes, this could save a lot of time
    reduced_cube_map = {bit: (None,None) for bit in set(known_bits) | set(guess_bits)}
    reduced_consts = np.zeros_like(consts,dtype=np.uint8)
    reduced_upper_matrix = np.eye(total_matrix.shape[0],dtype=np.uint8)
    reduced_lower_matrix = np.eye(total_matrix.shape[0],dtype=np.uint8)

    # Re-insert each equation and upper matrix to get new system:
    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Simplifying with known variables:")

    for eq_idx in range(len(total_matrix)):
        print(f"\r{indent(print_depth+2)}Equations simplified: {eq_idx+1}/{len(total_matrix)}", end='')
        
        # mark cube / time as None for known / guess bits
        if eq_idx in cube_map: cube,t = cube_map[eq_idx]
        else: cube,t = None,None
        
        linearly_independent = insert_equation(
            reduced_lower_matrix, reduced_upper_matrix, reduced_consts, reduced_cube_map,
            cube, t, total_matrix[eq_idx], consts[eq_idx]
        )

    if verbose:
        print(f"\n{indent(print_depth+1)}Finished Simplying:")
        print(f"{indent(print_depth+1)}Time: xxxxxx")

    # use the new reduced data going forward:
    lower_matrix = reduced_lower_matrix
    upper_matrix = reduced_upper_matrix
    total_matrix = (lower_matrix @ upper_matrix) % 2
    cube_map = reduced_cube_map
    consts = reduced_consts

    # fill in known values:
    known_values = np.zeros([state_size,1],dtype=np.uint8)
    cube_background = np.array([None] * state_size)
    for bit,val in known_bits.items():
        cube_background[bit] = val
        known_values[bit] = val
    
    # assume guess bits are 0 for the base cube calculation
    for bit in guess_bits:
        cube_background[bit] = 0

    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Summing cubes to generate equations:")

    # calculate the base cube values
    start_time = time.time()
    query_count = 0
    base_cube_values = np.zeros([state_size,1],dtype=np.uint8)


    if verbose:
        eq_idx = 0
        total_eqs = len([x for x in cube_map if x != None])

    # only calculate each cube once and re-use for different times:
    cube_cache = {}
    for bit, (cube,t) in cube_map.items():
        if verbose:
            eq_idx += 1
            print(f'\r{indent(print_depth+2)}Equations generated: {eq_idx}/{total_eqs}', end = '')

        if cube != None:
            if cube not in cube_cache:
                query_count += 2**len(cube)
                cube_cache[cube] = evaluate_super_poly(access_fn, cube, cube_background)
            base_cube_values[bit] = consts[bit] ^ cube_cache[cube][t]
    query_time = time.time() - start_time
    if verbose:
        print(f'\n{indent(print_depth+1)}Finished summing cubes:')
        print(f'{indent(print_depth+1)}Time: xxxxxx')
    
    # guess assignment of the guess bits and solve
    print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Starting to Guess:")
    guess_start_time = time.time()
    
    found = False
    guess_count = 0
    total_values = np.zeros([state_size,1],dtype=np.uint8)
    for assignment in product((0,1), repeat = len(guess_bits)):
        guess_count += 1
        if verbose:
            print(f"\r{indent(print_depth+2)}Guess count: {guess_count}",end='')

        # reset total values to default state (guess 0, known and base cube values in place):
        total_values[:] = known_values | base_cube_values

        # cubes are calculated with the ground truth state, 
        # this means changes to the guesses dont change the cube values
        # and changes need to be made to the value vector (other than setting the guess values)
        for i in range(len(assignment)):
            if assignment[i]:
                total_values[guess_bits[i]] ^= 1

        # Solve the matrix to recover the initial state:
        state_candidate = lu_solve(
            lower_matrix,
            upper_matrix,
            total_values
        )[:,0]
        
        # test if candidate is correct
        if test_fn(state_candidate):
            found = True
            #print("MATCH FOUND", state_candidate)
            break

    guess_time = time.time() - start_time
    if verbose:
        print(f"\n{indent(print_depth+1)}Guessing Finished:")
        print(f"{indent(print_depth+1)}Time: {time.time() - guess_start_time} s")

    if found:
        if verbose:
            print(f'{indent(print_depth+1)}\n{indent(print_depth+1)}Solution Found!')
            print(f'{indent(print_depth+1)}Total number of Queries: {query_count}')
            print(f'{indent(print_depth+1)}Time: {query_time}')
        return state_candidate
    else:
        return None










# returns a vector of outputs
def evaluate_super_poly(sim_fn, index_set, state, verbose=False):
    # input sanitization:
    state_copy = state.copy()

    # get the form of the cube:
    xor_total = np.zeros_like(sim_fn(state_copy))

    # sum over the cube:
    for assigment in list(product(range(2),repeat=len(index_set))):
        for n in range(len(assigment)):
            state_copy[index_set[n]] = assigment[n]
        a = sim_fn(state_copy.copy())
        if verbose: print(assigment, a)
        xor_total ^= a
    return xor_total

def get_nonlinear_mask(sim_fn, index_set, state_size, num_tests):
    offset = np.zeros(state_size,'uint8')
    nonlinear_mask = np.zeros_like(sim_fn(offset))

    for n in range(num_tests):
        state = np.random.randint(0,2,state_size,'uint8')
        delta = np.random.randint(0,2,state_size,'uint8')
        diff = state ^ delta

        # BLR test for a nonlinear relationship:
        nonlinear_mask |= (
            evaluate_super_poly(sim_fn,index_set,state.copy()) ^ 
            evaluate_super_poly(sim_fn,index_set,delta.copy()) ^
            evaluate_super_poly(sim_fn,index_set,diff.copy()) ^
            evaluate_super_poly(sim_fn,index_set,offset.copy())
        )

    return nonlinear_mask

def get_constant_mask(sim_fn, index_set, state_size, num_tests):
    state = np.zeros(state_size,'uint8')
    comparison_vector = evaluate_super_poly(sim_fn,index_set,state.copy())
    constant_mask = np.ones_like(comparison_vector)
    comparison_vector ^= 1

    for n in range(num_tests):
        state = np.random.randint(0,2,state_size,'uint8')
        constant_mask &= (
            comparison_vector ^
            evaluate_super_poly(sim_fn,index_set,state.copy())
        )

    return constant_mask

def determine_equations(fn, index_set, state_size, target_set = None):
    if target_set == None: target_set = range(state_size)

    state = np.zeros(state_size, dtype = np.uint8)
    consts = evaluate_super_poly(fn,index_set,state.copy())

    keystream_len = len(consts)
    coefs = np.zeros([state_size,keystream_len], dtype=np.uint8)
    
    for i in target_set:
        state[i] = 1
        coefs[i] = consts ^ evaluate_super_poly(fn,index_set,state.copy())
        state[i] = 0

    # transpose coefs to be [time, bit] instead.
    return coefs.T, consts

