"""Cube attacks on CMPRs, with equations of any degree the monomial profile allows.

Summing the output over every assignment of a set I of tweakable bits (a cube)
leaves the superpoly of I: a polynomial in the remaining bits, whose monomials
are those of the output that contain T_I, with T_I removed (see
`docs/theory/Cube Equation Generation.md` for the identity). Each keystream
position t gives one: an equation P_{I,t}(key) = S_{I,t}, where the offline
phase recovers P_{I,t} by cube sums on a simulated register and the online
phase measures S_{I,t} by the same cube sum on the target.

The monomial profile bounds each superpoly's degree before any register is run
(`MonomialProfile.get_cube_candidates`), so the offline phase computes every
coefficient up to that bound exactly and candidates are tried lowest degree
first. The resulting equations are ordinary polynomial equations over the
initial state, so the online phase feeds them to any equation store and solver,
as NAA, RAA and FAA do -- LU linearization for the linear ones, Groebner
reduction where the degree makes linearization wasteful.

The IV is public: the bits the attacker knows (`known_bits`) hold the same
values in both phases, so a superpoly is a polynomial in the unknown bits alone.
"""
import time
from itertools import chain, combinations, product

import numpy as np

from PyPR.BooleanLogic import AND, VAR, XOR

from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Cryptanalysis.Components.Adapters.online_insertion import make_online_inserter
from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.Cryptanalysis.Components.EquationStores.FilteringEqStore import (
    FilteringEqStore,
)
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore


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
    if getattr(register.fn, '_compiled', None) is None:
        register.fn.compile()
    if getattr(output_fn, '_compiled', None) is None:
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

    return access_fn, sim_fn, test_fn


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

    # lowest superpoly degree first, then smallest cube
    cube_candidates = sorted(
        output_profile.get_cube_candidates(),
        key = (lambda x: (x[3], sum(x[0].counts.values())))
    )

    print('Analyzing Cube Candidates:')
    for cube_profile, target_blocks, num_cubes, degree in cube_candidates:
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


        # only populated and only read when analyze_sources is set
        sources = []
        if analyze_sources:
            # the output terms that feed this superpoly: those whose profile has
            # a term containing the cube with variables to spare
            for term in output_anf.args:
                term_profile = term.remap_constants([
                    (0, MonomialProfile.logical_zero()),
                    (1, MonomialProfile.logical_one())
                ]).eval_ANF(monomial_profiles)

                is_source = False
                for output_monomial in term_profile.terms:
                    diffs = [
                        (output_monomial.counts[i] if i in output_monomial.counts else 0) -
                        (cube_profile.counts[i] if i in cube_profile.counts else 0)
                        for i in set(output_monomial.counts) | set(cube_profile.counts)
                    ]
                    if all(x >= 0 for x in diffs) and sum(diffs) > 0:
                        is_source = True
                if is_source:
                    sources.append(term)

        print('Cube Profile: ', cube_profile)
        print('Superpoly Degree (at most): ', degree)
        print('Target Blocks: ', target_blocks, '- Target Block Size:', sum(len(cmpr_fn.blocks[t]) for t in target_blocks))
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
    cmpr_fn, output_fn, sim_fn, tweakable_vars, known_bits,
    max_degree = None, time_limit = None, verbose = False, print_depth=0
    ):
    """Find cubes and recover their superpolys as equations on the unknown bits.

    Candidates come from the output's monomial profile, lowest superpoly degree
    first. For a cube I whose superpoly has degree at most d over the unknown
    bits U of its target blocks, the coefficient of each monomial x^m (m a subset
    of U, |m| <= d) follows by Moebius inversion from the superpoly's values at
    the points e_s, |s| <= d:

        P(e_s) = XOR over u subset of s of coef(u),  so
        coef(s) = P(e_s) XOR (XOR over u proper subset of s of coef(u)).

    Each P(e_s) is one cube sum on the simulated register with the known bits at
    their values and the bits of s set. Because the profile bounds the degree,
    these coefficients are the whole superpoly, not an approximation of it.

    An equation is kept when it is not constant and is linearly independent of
    those kept so far, judged by an `LUEqStore` over the linearized monomials.
    A candidate is skipped once every monomial it could produce is a pivot of
    that store: nothing it yields could be independent.

    :param cmpr_fn: The target register's feedback function.
    :type cmpr_fn: CMPR
    :param output_fn: The output function.
    :type output_fn: BooleanFunction
    :param sim_fn: Simulates the register from a full initial state, returning
        the keystream after the initialization rounds (from `access_fns`).
    :type sim_fn: Callable[[np.ndarray], np.ndarray]
    :param tweakable_vars: Bits the attacker can set, from which cubes are drawn.
    :type tweakable_vars: list[int]
    :param known_bits: Values of the bits the attacker knows -- at least every
        tweakable bit. The online phase must use the same values.
    :type known_bits: dict[int, int]
    :param max_degree: Skip candidates whose superpoly degree bound exceeds this.
        The cost of one cube grows as the number of monomials up to its degree,
        so this caps the work per cube. Defaults to no cap.
    :type max_degree: int | None
    :param time_limit: Seconds after which the search stops, defaults to None.
    :type time_limit: float | None
    :param verbose: Print progress, defaults to False.
    :type verbose: bool
    :param print_depth: Indentation level for verbose output, defaults to 0.
    :type print_depth: int
    :raises ValueError: If a tweakable bit has no known value.
    :return: attack data for `cube_attack_online`.
    :rtype: dict
    """
    unvalued = sorted(set(tweakable_vars) - set(known_bits))
    if unvalued:
        raise ValueError(
            f"tweakable bits {unvalued} have no value in known_bits: the non-cube "
            "tweakable bits are held at their known values in both phases"
        )

    print_skipped_cubes = False
    if verbose:
        print(f"{indent(print_depth)}Starting offline phase (Cube Attack):")

    start_time = time.time()

    # break up tweakable variables by block and compute cube candidates:
    tweakable_set = set(tweakable_vars)
    tweakable_blocks = [set(block) & tweakable_set for block in cmpr_fn.blocks]

    # the superpoly background: known bits at their values, unknown bits zero
    background = np.zeros(cmpr_fn.size, dtype=np.uint8)
    for bit, value in known_bits.items():
        background[bit] = value

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
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Calculating cube candidates:")

    cube_cand_time = time.time()

    # lowest superpoly degree first, then smallest cube
    cube_candidates = sorted(
        (candidate for candidate in monomial_profile.get_cube_candidates()
         if max_degree is None or candidate[3] <= max_degree),
        key = (lambda x: (x[3], sum(x[0].counts.values())))
    )

    if verbose:
        print(f"{indent(print_depth+1)}Cube candidates computed:")
        print(f"{indent(print_depth+1)}Time: {time.time() - cube_cand_time} s")
        print(f"{indent(print_depth+1)}")
        print(f"{indent(print_depth+1)}Identifying Cubes:")

    # judges independence only; the equations themselves are kept below
    rank_tracker = LUEqStore()
    equations = []

    # Maxterm search
    maxterm_count = 0
    # rebuilt per candidate when verbose; only ever printed under that guard
    profile_prefix = ""
    for cube_profile, target_blocks, num_cubes, degree in cube_candidates:
        if verbose:
           profile_prefix = f"{indent(print_depth+2)}Cube Candidate Profile: {cube_profile} (degree <= {degree})"

        # the unknown bits this candidate's superpolys can involve, and the
        # points e_s (|s| <= degree) that determine them, smallest first
        region_bits = sorted(
            bit for t in target_blocks for bit in cmpr_fn.blocks[t] if bit not in known_bits
        )
        points = [
            s for size in range(min(degree, len(region_bits)) + 1)
            for s in combinations(region_bits, size)
        ]

        # saturated: every monomial it could produce is already a pivot
        region_saturated = True
        for s in points[1:]:
            if not (s in rank_tracker.comb_to_idx and
                    rank_tracker.solved_for[rank_tracker.comb_to_idx[s]]):
                region_saturated = False
                break

        # create the iterators and calculate some statistics:
        tweakable_cube_count = 1
        variable_iterators = []
        # zipped with variable_iterators below, outside the guard that used
        # to bind this, so a saturated region reached that line unbound
        loop_nums = []

        # only compute tweakable bits for regions which are not saturated
        if not region_saturated:
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
        if region_saturated:
            if print_skipped_cubes and verbose:
                print(f'{profile_prefix}: skipped (target region already saturated)')
            continue

        # test the individual cubes/maxterms:
        already_printed = False
        for cube_idx, var_selections in enumerate(iproduct(*variable_iterators)):
            if region_saturated:
                break
            if time_limit and time.time() - start_time > time_limit:
                break

            # print only inside the loop to make sure there are actual cubes
            # depending on the tweakable set, this iterator may be empty
            if not already_printed:
                already_printed = True
                if verbose:
                    print(f'{profile_prefix}:')
                    print(f'{indent(print_depth+2)}Target Blocks: {target_blocks} - Unknown Bits: {len(region_bits)}')
                    print(f'{indent(print_depth+2)}Number of Cube Candidates (before restriction): {num_cubes}')
                    print(f'{indent(print_depth+2)}Number of Cube Candidates (restricted to tweakable bits): {tweakable_cube_count}')

            maxterm_count += 1
            maxterm = tuple(chain(*var_selections))
            if verbose:
                print(f'\r{indent(print_depth+3)}Cube {cube_idx+1}/{tweakable_cube_count}: {maxterm}',end='')

            # superpoly values at each point, then Moebius inversion to the
            # coefficients, both as vectors over the keystream positions
            coefs = {}
            for s in points:
                state = background.copy()
                for bit in s:
                    state[bit] = 1
                coef = evaluate_super_poly(sim_fn, maxterm, state)
                for size in range(len(s)):
                    for u in combinations(s, size):
                        coef ^= coefs[u]
                coefs[s] = coef

            for t in range(len(coefs[()])):
                monomials = [s for s in points[1:] if coefs[s][t]]
                # a constant superpoly says nothing about the unknown bits
                if not monomials:
                    continue

                equation = XOR(*[AND(*[VAR(bit) for bit in s]) for s in monomials])
                if rank_tracker.insert_equation(equation, identifier=(maxterm, t), translate_ANF=False):
                    equations.append((maxterm, t, monomials, int(coefs[()][t])))

                    region_saturated = True
                    for s in points[1:]:
                        if not (s in rank_tracker.comb_to_idx and
                                rank_tracker.solved_for[rank_tracker.comb_to_idx[s]]):
                            region_saturated = False
                            break
                    if region_saturated:
                        if verbose: print(f'\n{indent(print_depth+2)}Target Region Saturated!')
                        break

                if time_limit and time.time() - start_time > time_limit:
                    break

        # This check breaks out of the monomial profile loop
        # no saturation check because regions are profile-specific
        if time_limit and time.time() - start_time > time_limit:
            break

    # every state bit is a variable, so solvers can read the state back, and the
    # constant is a column, so a consistent store can hold the right-hand sides
    comb_to_idx: dict[tuple[int, ...], int] = {(): 0}
    for bit in range(cmpr_fn.size):
        comb_to_idx[(bit,)] = len(comb_to_idx)
    for _maxterm, _t, monomials, _constant in equations:
        for s in monomials:
            if s not in comb_to_idx:
                comb_to_idx[s] = len(comb_to_idx)

    num_queries = 0
    distinct_cubes = set()
    for maxterm, _t, _monomials, _constant in equations:
        if maxterm not in distinct_cubes:
            num_queries += 2**len(maxterm)
            distinct_cubes.add(maxterm)

    if verbose:
        print(f'{indent(print_depth+1)}Finished equation generation: ')
        print(f'{indent(print_depth+1)}Time: {time.time() - cube_cand_time} s')
        print(f'{indent(print_depth+1)}')
        print(f'{indent(print_depth+1)}Number of equations found: {len(equations)}')
        print(f'{indent(print_depth+1)}Number of cubes tested: {maxterm_count}')
        print(f'{indent(print_depth+1)}Number of queries in attack: {num_queries}')
        print('Offline phase complete -- Total time: ', time.time() - start_time)

    output = {}
    output['equations'] = equations
    output['known bits'] = dict(known_bits)
    output['comb to idx map'] = comb_to_idx
    output['idx to comb map'] = {idx: comb for comb, idx in comb_to_idx.items()}
    output['num variables'] = len(comb_to_idx)
    return output


def cube_attack_online(
    feedback_fn, output_fn, access_fn, test_fn, attack_data,
    time_limit=None, verbose=False, print_depth=0,
    solver=None, online_store=None,
):
    """Measure each cube on the target and solve the resulting equations.

    For every equation P_{I,t}(x) = S_{I,t} from the offline phase, S_{I,t} is
    the cube sum over I of the target's keystream at position t, with the known
    bits at the values the offline phase used. The known bits enter as the
    equations x_i = v_i, so the store's system pins the whole initial state and
    any store and solver pairing that RAA and FAA accept works here.

    The target is reached only through `access_fn` and `test_fn`, as an
    attacker would reach it: the solver checks each candidate state with
    `test_fn`, so the target's keystream and initialization rounds never leave
    the interface.

    :param feedback_fn: The target register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The output function.
    :type output_fn: BooleanFunction
    :param access_fn: Queries the target with the tweakable bits set as given
        (from `access_fns`).
    :type access_fn: Callable[[np.ndarray], np.ndarray]
    :param test_fn: Whether a candidate initial state reproduces the target's
        keystream (from `access_fns`).
    :type test_fn: Callable[[np.ndarray], bool]
    :param attack_data: The output of `cmpr_cube_attack_offline`.
    :type attack_data: dict
    :param time_limit: Seconds after which to give up, defaults to None.
    :type time_limit: float | None
    :param verbose: Print progress, defaults to False.
    :type verbose: bool
    :param print_depth: Indentation level for verbose output, defaults to 0.
    :type print_depth: int
    :param solver: Defaults to `LUSolver()`.
    :type solver: LUSolver | GaussElimSolver | GrobnerSolver | SplitGrobnerSolver | None
    :param online_store: Defaults to a consistent `LUEqStore` over the attack's monomials.
    :type online_store: BaseEqStore | None
    :raises ValueError: If the offline phase found no equations.
    :return: The recovered initial state, or None.
    :rtype: list[int] | None
    """
    if verbose:
        print(f"{indent(print_depth)}Starting online phase (Cube Attack):")
    start_time = time.time()

    equations = attack_data['equations']
    known_bits = attack_data['known bits']
    comb_to_idx = attack_data['comb to idx map']
    idx_to_comb = attack_data['idx to comb map']
    num_vars = attack_data['num variables']

    # if no cubes, then cube attack is slower than brute force:
    # exit immediately
    if not equations:
        raise ValueError(
            'No cubes given; consider either providing ' +
            'cubes for the attack or a brute force approach'
        )

    if solver is None:
        solver = LUSolver()
    if online_store is None:
        online_store = LUEqStore(comb_to_idx, consistent=True)

    # the known IV; access_fn writes only the tweakable bits, so the other
    # known bits come from the target as they are
    background = np.array([None] * feedback_fn.size)
    for bit, value in known_bits.items():
        background[bit] = value

    if verbose:
        print(f"{indent(print_depth+1)}Summing cubes to generate equations:")

    insert_eq, finalize = make_online_inserter(
        online_store, idx_to_comb,
        total_eqs=len(known_bits) + len(equations), num_vars=num_vars,
        verbose=verbose, print_depth=print_depth+2,
    )

    eq_idx = 0
    determined = False
    for bit, value in sorted(known_bits.items()):
        coef_vector = np.zeros([num_vars], dtype=np.uint8)
        coef_vector[comb_to_idx[(bit,)]] = 1
        coef_vector[comb_to_idx[()]] = value
        determined = insert_eq(coef_vector, eq_idx)
        eq_idx += 1
        if determined:
            break

    # only calculate each cube once and re-use for different times:
    query_count = 0
    cube_sums = {}
    if not determined:
        for maxterm, t, monomials, constant in equations:
            if maxterm not in cube_sums:
                query_count += 2**len(maxterm)
                cube_sums[maxterm] = evaluate_super_poly(access_fn, maxterm, background)

            coef_vector = np.zeros([num_vars], dtype=np.uint8)
            for s in monomials:
                coef_vector[comb_to_idx[s]] = 1
            coef_vector[comb_to_idx[()]] = constant ^ cube_sums[maxterm][t]

            if insert_eq(coef_vector, eq_idx):
                break
            eq_idx += 1
            if time_limit and (time.time() - start_time >= time_limit):
                if verbose:
                    print(f"\n{indent(print_depth+1)}Time limit reached during cube summation.")
                break

    finalize()

    if verbose:
        if isinstance(online_store, FilteringEqStore):
            solved_count = online_store.num_determined
        else:
            solved_count = online_store.num_eqs
        print(f"\n{indent(print_depth+1)}Finished summing cubes:")
        print(f"{indent(print_depth+1)}Queries: {query_count}")
        print(f"{indent(print_depth+1)}Variables Solved: {solved_count}/{num_vars}")
        print(f"{indent(print_depth+1)}Time: {time.time() - start_time} s")

    if time_limit and (time.time() - start_time >= time_limit):
        if verbose:
            print(f"{indent(print_depth)}Online phase timed out -- Total time: {time.time() - start_time} s")
        return None

    if verbose:
        print(f"{indent(print_depth+1)}\n{indent(print_depth+1)}Starting solve:")

    initial_state, _, _ = solver.solve(
        online_store, feedback_fn, output_fn, None,
        verify=test_fn,
        verbose=verbose, _print_depth=print_depth+1,
    )

    if verbose:
        print(f"{indent(print_depth)}Online phase complete -- Total time: {time.time() - start_time} s")

    return initial_state


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
