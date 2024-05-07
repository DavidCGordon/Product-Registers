import time

import numba
import numpy as np
from itertools import combinations,chain,cycle,product,tee

from ProductRegisters import FeedbackRegister
from ProductRegisters.BooleanLogic import BooleanFunction

#TODO: Memory reuse optimizations to prevent loads of unecessary copying/writing
#TODO: Accept more general functions / implement a general cube attack.

def CMPR_access_fns(register, tweakable_bits, init_rounds=100,keystream_len=100):
    if not hasattr(register.fn,'_compiled'):
        register.fn.compile()

    for i in range(init_rounds):
        register.clock_compiled()
    keystream = [state[0] for state in register.run_compiled(keystream_len)]
    register.reset()

    # given a full state, simulate that state to get the bit
    def sim_fn(state):
        register._state = state
        for i in range(init_rounds):
            register.clock_compiled()
            output = register._state[0]

        register.reset()
        return output

    def access_fn(state):
        register.reset()
        # write only to tweakable bits
        for bit in tweakable_bits:
            if state[bit] != None:
                register._state[bit] = state[bit]

        for i in range(init_rounds):
            register.clock_compiled()
            output = register._state[0]

        register.reset()
        return output
    
    # make our own register for simulation
    def test_fn(state):
        register._state = state
        for i in range(init_rounds):
            register.clock_compiled()
        test_keystream = [s[0] for s in register.run_compiled(keystream_len)]
        register.reset()
        return test_keystream == keystream

    return access_fn,sim_fn,test_fn












def cmpr_cube_summary(cmpr_fn,tweakable_vars):
    print('Beginning Summary:')

    tweakable_set = set(tweakable_vars)
    tweakable_counts = [len(set(block) & tweakable_set) for block in cmpr_fn.blocks]

    print('Computing Monomial Profile')
    cube_candidates = sorted(
        cmpr_fn.monomial_profiles()[0].get_cube_candidates(),
        key = (lambda x: sum(c for c in x[0].counts.values()))
    )

    print('Verifying Cube Candidates:')
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
        if tweakable_cube_count > 0:
            print('Cube Profile: ', cube_profile)
            print('Target Block: ', target_block, 'Target Block Size:', len(cmpr_fn.blocks[target_block]))
            print('Number of Cube Candidates (before restriction): ', num_cubes)
            print('Number of Cube Candidates (restricted to tweakable bits): ', tweakable_cube_count)
            print('Cube Failure Probability: ', cube_failure_prob)
            print('\n')

    print("Summary Finished!")


# allows for faster skipping of unusable sets :)
# attribution: https://discuss.python.org/t/a-product-function-which-supports-large-infinite-iterables/5753
def iproduct(*iterables, repeat=1):
    iterables = [item for row in zip(*(tee(iterable, repeat) for iterable in iterables)) for item in row]
    N = len(iterables)
    saved = [[] for _ in range(N)]  # All the items that we have seen of each iterable.
    exhausted = set()               # The set of indices of iterables that have been exhausted.
    for i in cycle(range(N)):
        if i in exhausted:  # Just to avoid repeatedly hitting that exception.
            continue
        try:
            item = next(iterables[i])
            yield from product(*saved[:i], [item], *saved[i+1:])  # Finite product.
            saved[i].append(item)
        except StopIteration:
            exhausted.add(i)
            if not saved[i] or len(exhausted) == N:  # Product is empty or all iterables exhausted.
                return
    yield ()  # There are no iterables.
            

def cmpr_cube_attack_offline(cmpr_fn, sim_fn, tweakable_vars, time_limit = None, num_tests = 20, verbose = False):
    start_time = time.time()
    
    # create a state vector for coef/constant calculations:
    # doesnt matter what it is, but this avoids a ton of creating new vectors:
    coef_state = np.zeros(cmpr_fn.size,dtype=np.uint8)
    # create vars to hold outputs
    
    cube_map = {}
    lower_matrix = np.eye(cmpr_fn.size,dtype=np.uint8)
    upper_matrix = np.eye(cmpr_fn.size,dtype=np.uint8)
    const_vec = np.zeros([cmpr_fn.size,1],dtype=np.uint8)

    # break up tweakable variables by block and compute cube candidates:
    tweakable_set = set(tweakable_vars)
    tweakable_blocks = [set(block) & tweakable_set for block in cmpr_fn.blocks]
    cube_candidates = sorted(
        cmpr_fn.monomial_profiles()[0].get_cube_candidates(),
        key = (lambda x: sum(x[0].counts.values()))
    )

    # Maxterm search
    count = 0
    for cube_profile, target_block, num_cubes, cube_failure_prob in cube_candidates:
        if verbose: print('\nCube Profile: ', cube_profile)

        # check to see if cube is less efficient than brute force:
        # if sum(cube_profile.counts.values()) >= len(cmpr_fn.blocks[target_block]):
        #     if verbose: print(' - Cube skipped (less efficient than brute force)')
        #     continue

        # check to see if the block is saturated:
        block_already_saturated = all([(bit in cube_map) for bit in cmpr_fn.blocks[target_block]])

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
            if verbose: print(' - Cube skipped (not possible with current tweakable bits)')
            continue
        
        # output message for cube profiles we won't use but could:
        if block_already_saturated:
            if verbose: print(' - Cube skipped (target block already saturated)')
            continue

        already_printed = False
        # test the cubes:
        for var_selections in iproduct(*variable_iterators):
            if time_limit and time.time() - start_time > time_limit:
                break 
            
            # print only inside the loop to make sure there are actual cubes
            # depending on the tweakable set, this iterator may be empty
            if not already_printed:
                already_printed = True
                if verbose:
                    print('Target Block: ', target_block, 'Target Block Size:', len(cmpr_fn.blocks[target_block]))
                    print('Number of Cube Candidates (before restriction): ', num_cubes)
                    print('Number of Cube Candidates (restricted to tweakable bits): ', tweakable_cube_count)
                    print('Cube Failure Probability: ', cube_failure_prob)
            
            count += 1
            maxterm = list(chain(*var_selections))
            if verbose: print('Cube Candidate: ', maxterm, end='')

            if not is_useful(sim_fn,maxterm,cmpr_fn.size,num_tests):
                if verbose: print('\t -  Superpoly not linear')
                continue

            # calculate coefficients
            coef_vector, const =  determine_equation(sim_fn,maxterm,coef_state,target_set=cmpr_fn.blocks[target_block])

            # calculate lower and upper entries
            linearly_independent = False
            modification_vector  = np.zeros_like(coef_vector)
            for bit in range(len(coef_vector)):
                if coef_vector[bit] == 1:
                    modification_vector[bit] = 1
                    if bit in cube_map:
                        coef_vector ^= upper_matrix[bit]
                    else:
                        linearly_independent = True
                        cube_map[bit] = maxterm
                        const_vec[bit] = const
                        lower_matrix[bit] = modification_vector
                        upper_matrix[bit] = coef_vector
                        if verbose: 
                            print("\tTotal: ", len(cube_map))
                        break
            if verbose: 
                if not linearly_independent:
                    print('\t -  Not linearly independent')

            # if you added a pivot, check to see if the block is saturated now:
            block_already_saturated = all([(bit in cube_map) for bit in cmpr_fn.blocks[target_block]])
            if block_already_saturated:
                break
        
        # second time check to avoid printing everything
        if time_limit and time.time() - start_time > time_limit:
                break 
        
    num_queries = 0
    for cube in cube_map.values():
        num_queries += 2**len(cube)

    if verbose:  
        print("Number of cubes tested: ", count)
        print("Number of cubes found: ", len(cube_map))
        print("Num Queries: ", num_queries)

    output = {}
    output['cubes'] = cube_map
    output['equation matrix'] = (lower_matrix @ upper_matrix) % 2
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








def cube_attack_online(access_fn, test_fn, state_size, known_bits, cube_data, verbose = False):
    cube_map = cube_data['cubes']
    total_matrix = cube_data['equation matrix']
    consts = cube_data ['constant vector']

    # create copies to prevent known-variable reduction from deleting important information
    secret_bits = [i for i in range(state_size) if i not in known_bits]
    guess_bits = [i for i in secret_bits if i not in cube_map]
    if verbose: print("Guessing Bits: ", guess_bits)

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
    reduced_cube_map = {bit: None for bit in set(known_bits) | set(guess_bits)}
    reduced_consts = np.zeros_like(consts,dtype=np.uint8)
    upper_matrix = np.eye(total_matrix.shape[0],dtype=np.uint8)
    lower_matrix = np.eye(total_matrix.shape[0],dtype=np.uint8)

    coef_vector = np.zeros(state_size,dtype=np.uint8)
    modification_vector  = np.zeros_like(coef_vector)   

    # Rewrite lower and upper matrix to get new system:
    for eq_idx in range(len(total_matrix)):
        coef_vector[:] = total_matrix[eq_idx]
        modification_vector[:] = 0

        for bit in range(state_size):
            if coef_vector[bit] == 1:
                modification_vector[bit] = 1
                if bit in reduced_cube_map:
                    coef_vector ^= upper_matrix[bit]
                else:
                    reduced_cube_map[bit] = cube_map[eq_idx]
                    reduced_consts[bit] = consts[eq_idx]
                    lower_matrix[bit] = modification_vector
                    upper_matrix[bit] = coef_vector
                    # print(coef_vector,modification_vector)
                    break

    # use the new reduced data going forward:
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

    # calculate the base cube values
    start_time = time.time()
    query_count = 0
    base_cube_values = np.zeros([state_size,1],dtype=np.uint8)
    for bit,cube in cube_map.items():
        if cube != None:
            query_count += 2**len(cube)
            base_cube_values[bit] = consts[bit] ^ evaluate_super_poly(access_fn, cube, cube_background)
    query_time = time.time() - start_time
    

    # guess assignment of guess_bits and solve
    start_time = time.time()

    found = False
    guess_count = 0
    total_values = np.zeros([state_size,1],dtype=np.uint8)
    for assignment in product((0,1), repeat = len(guess_bits)):
        guess_count += 1

        # reset total values to default state (guess 0, known and base cube values in place):
        total_values[:] = known_values | base_cube_values

        # because cubes are calculated with guesses = 0, we have to account for this by adding the whole column
        # this flips both the guess bit, but also all cube values which depended on it, leading to correct values
        for i in range(len(assignment)):
            if assignment[i]:
                total_values ^= total_matrix[:,np.newaxis,guess_bits[i]]

        # Solve the matrix to recover the initial state:
        state_candidate = lu_solve(
            lower_matrix,
            upper_matrix,
            total_values
        )[:,0]
        
        # test if candidate is correct
        # print("Values: ",total_values.T[0])
        # print("Answer: ",state_candidate.T, '\n')
        if test_fn(state_candidate):
            found = True
            break

    guess_time = time.time() - start_time
    if verbose:
        print('Query count:\t', query_count, '\tQuery time: ', query_time)
        print('Guess count:\t', guess_count, '\tGuess time: ', guess_time)

    if found:
        return state_candidate
    else:
        return None,








    
# for efficiency, the user must ensure fn is callable:
def evaluate_super_poly(fn, index_set, state):
    # input sanitization:
    state_copy = state.copy()

    # sum over the cube:
    xor_total = 0
    for assigment in list(product(range(2),repeat=len(index_set))):
        for n in range(len(assigment)):
            state_copy[index_set[n]] = assigment[n]
        xor_total ^= fn(state_copy)
    return xor_total


def is_constant(fn, index_set, state_size, num_tests):
    outputs = []
    for n in range(num_tests):
        state = np.random.randint(0,2, state_size,'uint8')
        outputs.append(evaluate_super_poly(
            fn,
            index_set,
            state
        ))

    return (np.any(outputs) and not np.all(outputs))

def is_linear(fn, index_set, state_size, num_tests):
    for n in range(num_tests):
        offset = np.zeros(state_size,'uint8')
        state = np.random.randint(0,2,state_size,'uint8')
        delta = np.random.randint(0,2,state_size,'uint8')
        diff = state ^ delta

        # test for a nonlinear relationship, if found, break early
        if (
            evaluate_super_poly(fn,index_set,state) ^
            evaluate_super_poly(fn,index_set,delta) ^
            evaluate_super_poly(fn,index_set,diff) ^
            evaluate_super_poly(fn,index_set,offset)
        ):
            return False
    return True


# combines the above tests to be slightly more efficient:
def is_useful(fn, index_set, state_size, num_tests):
    val = None
    constant = True

    for n in range(num_tests):
        offset = np.zeros(state_size,'uint8')
        state = np.random.randint(0,2,state_size,'uint8')
        delta = np.random.randint(0,2,state_size,'uint8')
        diff = state ^ delta

        # test for a nonlinear relationship, if found, break early (False)
        output = evaluate_super_poly(fn,index_set,state)
        if (
            output ^ 
            evaluate_super_poly(fn,index_set,delta) ^
            evaluate_super_poly(fn,index_set,diff) ^
            evaluate_super_poly(fn,index_set,offset)
        ):
            return False
        
        # test for non-constant relationship, if found, break early (True)
        if val == None:
            val = output
        elif output != val:
            constant = False

    # This is reached if the function is linear or a constant:
    return not(constant)










def determine_coefficient(fn, index_set, state, bit):
    state = np.zeros_like(state)
    first_value = evaluate_super_poly(fn, index_set, state)
    state[bit] = 1
    second_value = evaluate_super_poly(fn, index_set, state)
    return first_value ^ second_value

def determine_constant(fn, index_set, state):
    return evaluate_super_poly(
        fn,index_set,np.zeros_like(state)
    )


# combines the above tests to be slightly more efficient:
def determine_equation(fn,index_set,state, target_set = None):
    if target_set == None: target_set = range(len(state))
    state_copy = np.zeros_like(state)
    coefs = np.zeros_like(state)
    const = evaluate_super_poly(fn,index_set,np.zeros_like(state))

    for i in target_set:
        state_copy[i] = 1
        coefs[i] = const ^ evaluate_super_poly(fn,index_set,state_copy)
        state_copy[i] = 0

    return coefs, const

