from typing import Any

from PyPR import FeedbackRegister
from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.BooleanLogic import BooleanFunction

from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey_iterator

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import CubeEqGenerator, get_var_map
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import SubstitutionEqGenerator

from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver

import numpy as np
import numba
import time
import random

u8 = numba.types.uint8
u64 = numba.types.uint64

# small helper function to help pretty-print:
def indent(n:int) -> str:
    return ("|   " * n)



def FAA_offline(
    feedback_fn: FeedbackFunction, 
    annihilator: BooleanFunction, 
    multiple: BooleanFunction, 
    init_rounds: int, 
    margin: int,
    time_limit: int, 
    verbose: bool = False,
    _print_depth: int = 0,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles: list[MonomialProfile] | None = None,
    variable_blocks: list[list[int]] | None = None
):

    if verbose:
        print(f"{indent(_print_depth)}Starting offline phase (Fast Algebraic Attack):")
    start_time = time.time()

    # compute equations for the annihilator:
    if monomial_profiles != None and variable_blocks != None:

        if verbose:
            print(f"{indent(_print_depth+1)}using monomial profile optimization: True")
            print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Calculating monomial profile for annihilator:")
            mp_a_time = time.time()

        annihilator_mp = annihilator.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)

        if verbose:
            print(f"{indent(_print_depth+1)}Monomial profile computed:")
            print(f"{indent(_print_depth+1)}Time: {time.time() - mp_a_time} s")
            print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Calculating monomial profile for low degree multiple:")
            mp_m_time = time.time()

        # Precompute LC for low degree multiple:
        multiple_mp = multiple.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)
        max_LC = multiple_mp.upper()

        if verbose:
            print(f"{indent(_print_depth+1)}Monomial profile computed:")
            print(f"{indent(_print_depth+1)}Time: {time.time() - mp_m_time} s")
            print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Calculating variable_map:")
            var_map_time = time.time()
        
        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, annihilator_mp, variable_blocks, complete_subsets = True
        )

        if verbose:
            print(f"{indent(_print_depth+1)}Variable map computed:")
            print(f"{indent(_print_depth+1)}Time: {time.time() - var_map_time} s")
            print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Calculating linear relation:")
            lin_rel_time = time.time()
       
        # use berlekamp_massey to get the exact relation
        feedback_fn.compile()
        multiple_compiled = multiple.compile()
        test_register = FeedbackRegister(random.randint(0,2**feedback_fn.size-1), feedback_fn)
        max_count = 1000*((2*max_LC+256)//1000 + 1)
        count = 0

        for curr_LC, curr_relation in berlekamp_massey_iterator(
            seq = (multiple_compiled(state._state) for state in test_register.run(2*max_LC+256)),
            yield_rate=1000
        ):
            count += 1000
            if verbose:
                print(
                    f"\r{indent(_print_depth+2)}Bits processed: {count} / {max_count}" +
                    f"  --  Linear Complexity: {curr_LC} / {max_LC}", 
                    end=''
                )

            linear_complexity = curr_LC
            linear_relation = curr_relation



        # flip linear relation, due to dot product vs convolution
        linear_relation = linear_relation[::-1]
        margin += linear_complexity

        if verbose:
            print(f"\n{indent(_print_depth+1)}Linear relation found:")
            print(f"{indent(_print_depth+1)}Linear complexity: {linear_complexity}")
            print(f"{indent(_print_depth+1)}Time: {time.time()-lin_rel_time} s")

        # use precomputed maps for faster eq generation and storage
        annihilator_eqs = EqStore(variable_indices)
        eq_gen = CubeEqGenerator(
            feedback_fn, annihilator, (len(variable_indices) + margin), 
            variable_indices, verbose=True, _print_depth=_print_depth+2
        )

        # additional flags
        check_ranks = False
        max_LC = len(variable_indices)

    else:
        if verbose:
            print(f"{indent(_print_depth+1)}Using monomial profile optimization: False")

        # create dynamic storage and generation
        annihilator_eqs = EqStore()
        annihilator_LU = LUEqStore()
        annihilator_LU.link(annihilator_eqs)
        
        # ensure all variables are in the eq store:
        for v in range(len(feedback_fn)):
            annihilator_eqs._update_known_monomials(tuple([v]))
            annihilator_LU._update_known_monomials(tuple([v]))

        eq_gen = SubstitutionEqGenerator(
            feedback_fn, annihilator, 2**feedback_fn.size
        )

        # have to check ranks, since number of
        # variables isnt known ahead of time
        check_ranks = True
        count_into_margin = 0

        # Precompute LC for low degree multiple:
        # because max_LC isnt known, test until there are no changes:
        feedback_fn.compile()
        test_register = FeedbackRegister(random.getrandbits(feedback_fn.size), feedback_fn)
        

        if verbose:
            print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Calculating linear complexity dynamically:")
            lin_rel_time = time.time()

        count = 0
        curr_LC = 0
        curr_relation = []
        for linear_complexity, linear_relation in berlekamp_massey_iterator(
            seq = (multiple.eval(state) for state in test_register.run(2**(feedback_fn.size))),
            yield_rate=1000
        ):
            if verbose:
                print(f"\r{indent(_print_depth+2)}Bits processed (thousands): {count} -- Linear Complexity: {curr_LC}", end='')

            # check lengths first for more efficient short circuit:
            if (linear_complexity == curr_LC) and (len(linear_relation) == len(curr_relation)) and np.all(linear_relation == curr_relation):
                break

            count += 1
            curr_LC = linear_complexity
            curr_relation = linear_relation
        
        # flip linear relation, due to dot product vs convolution
        linear_relation = linear_relation[::-1]
        margin += linear_complexity

        if verbose:
            print(f"\n{indent(_print_depth+1)}Linear relation found:")
            print(f"{indent(_print_depth+1)}Linear complexity: {curr_LC}")
            print(f"{indent(_print_depth+1)}Time: {time.time()-lin_rel_time} s")

    
    if verbose:
        print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Generating Equations:")
        eq_time = time.time()
   
    # main equation loop
    for t, ann_eq in enumerate(eq_gen): #type: ignore  (to narrow types correctly)
        ann_eq: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]

        # don't generate equations for initialization rounds
        if t < init_rounds: continue
        
        annihilator_eqs.insert_equation(ann_eq, identifier = t)

        if verbose: 
            print(f'\r{indent(_print_depth+2)}Equations Found: {annihilator_eqs.num_eqs} / {annihilator_eqs.num_vars + margin}',end='')

        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print(f"\n{indent(_print_depth)}Time limit reached!")
            break

        # break step only necessary for dynamic stores:
        # reduces speed a fair bit, due to extra insert
        if check_ranks:
            ann_independent = annihilator_LU.insert_equation(ann_eq, identifier = t)
        
            # continue for margin more steps after hitting linear 
            # recurrent phase (not perfect but better than nothing)
            if not (ann_independent):
                count_into_margin += 1
                if count_into_margin == margin:
                    break

    if verbose:
        print(f'\r{indent(_print_depth+2)}Equations Found: {annihilator_eqs.num_eqs} / {annihilator_eqs.num_vars + margin}',end='\n')
        print(f"{indent(_print_depth+1)}Finished equation generation: ")
        print(f"{indent(_print_depth+1)}Time: {time.time() - eq_time} s")
        print(f"Offline phase complete -- Total time: ", time.time() - start_time)

    output = {}
    output['idx to comb map'] = annihilator_eqs.idx_to_comb
    output['comb to idx map'] = annihilator_eqs.comb_to_idx
    output['annihilator equations'] = annihilator_eqs.equations[:annihilator_eqs.num_eqs,:annihilator_eqs.num_vars]
    #output['annihilator consts'] = annihilator_eqs.constants[:annihilator_eqs.num_eqs]
    output['linear relation'] = linear_relation
    output['num variables'] = annihilator_eqs.num_vars
    output['keystream needed'] = max(annihilator_eqs.equation_ids.values()) + 1
    output['margin'] = margin - (linear_complexity)

    return output


@numba.njit(u8[:](u64,u8[:],u8[:,:],u8[:]))
def sum_over_linear_relationship(
    start_idx: int, 
    keystream: np.ndarray[tuple[int],np.dtype[np.uint8]], 
    equations: np.ndarray[tuple[int,int],np.dtype[np.uint8]], 
    linear_relation: np.ndarray[tuple[int],np.dtype[np.uint8]]
):
    coef_vector = np.zeros((equations.shape[1],), dtype="uint8")

    for i in range(len(linear_relation)):
        if (keystream[start_idx + i])==1 and (linear_relation[i]==1):
            for j in range(equations.shape[1]):
                coef_vector[j] ^= equations[start_idx + i, j]

    return coef_vector


# Dont need known bits: this is because each equation is cheap (relative to cube attacks)
# and the known bits doesnt /really/ help with the monomials (without a big loop), so it
# doesnt shrink the system that much, but does introduce a lot of overhead.
def FAA_online(
    feedback_fn: FeedbackFunction,
    output_fn: BooleanFunction,
    keystream: list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]],
    attack_data: dict[str,Any],
    test_length: int = 1000,
    time_limit: float | None = None,
    verbose: bool = False,
    _print_depth: int = 0,
    solver=None,
    online_store=None,
):
    if solver is None:
        solver = LUSolver()

    if verbose:
        print(f"{indent(_print_depth)}Starting online phase (Fast Algebraic Attack):")
    start_time = time.time()

    if type(keystream) == np.ndarray:
        pass
    elif type(keystream) == list:
        keystream = np.array(keystream, dtype = 'uint8')
    else:
        raise ValueError(f"Keystream must be a list or u8 ndarray, not {type(keystream)}")

    # unpack attack_data
    num_vars = attack_data['num variables']
    num_eqs = attack_data['keystream needed']

    annihilator_eqs = attack_data['annihilator equations']
    linear_relation = attack_data['linear relation']

    comb_to_idx = attack_data['comb to idx map']
    idx_to_comb = attack_data['idx to comb map']

    if online_store is None:
        online_store = LUEqStore(comb_to_idx, consistent=True)

    if verbose:
        print(f"{indent(_print_depth+1)}Starting Equation Substitution:")

    from PyPR.Cryptanalysis.Components.Adapters.online_insertion import make_online_inserter

    total_online_eqs = num_eqs - len(linear_relation)
    insert_eq, finalize = make_online_inserter(
        online_store, idx_to_comb,
        total_eqs=total_online_eqs, num_vars=num_vars,
        verbose=verbose, print_depth=_print_depth+2,
    )

    for eq_idx in range(total_online_eqs):
        coef_vector = sum_over_linear_relationship(
            eq_idx, keystream, annihilator_eqs, linear_relation
        )
        if insert_eq(coef_vector, eq_idx):
            break
        if time_limit and (time.time() - start_time >= time_limit):
            if verbose:
                print(f"\n{indent(_print_depth+1)}Time limit reached during substitution.")
            break

    finalize()

    if verbose:
        if hasattr(online_store, 'solved_vars') and isinstance(online_store.solved_vars, dict):
            solved_count = len(online_store.solved_vars)
        else:
            solved_count = online_store.num_eqs
        print(f"\n{indent(_print_depth+1)}Finished substituting key stream:")
        print(f"{indent(_print_depth+1)}Variables Solved: {solved_count}/{num_vars}")
        print(f"{indent(_print_depth+1)}Time: {time.time() - start_time} s")

    if time_limit and (time.time() - start_time >= time_limit):
        if verbose:
            print(f"{indent(_print_depth)}Online phase timed out -- Total time: {time.time() - start_time} s")
        return None

    if verbose:
        print(f"{indent(_print_depth+1)}\n{indent(_print_depth+1)}Starting solve:")

    initial_state, _, _ = solver.solve(
        online_store, feedback_fn, output_fn, keystream,
        test_length=test_length,
        verbose=verbose, _print_depth=_print_depth+1,
    )

    if verbose:
        print(f"{indent(_print_depth)}Online phase complete -- Total time: {time.time() - start_time} s")

    return initial_state
 