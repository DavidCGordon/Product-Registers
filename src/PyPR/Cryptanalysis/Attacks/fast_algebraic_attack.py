import random
import time
from dataclasses import dataclass

import numba
import numpy as np
from numba.core import types as nb_types

from PyPR import FeedbackRegister
from PyPR.Reporting import format_duration, get_logger

from PyPR.BooleanLogic import BooleanFunction

from PyPR.FeedbackFunctions import FeedbackFunction

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey_iterator
from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile

from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import (
    CubeEqGenerator,
    get_var_map,
)
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import (
    SubstitutionEqGenerator,
)
from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.FilteringEqStore import (
    FilteringEqStore,
)
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

log = get_logger(__name__)

u8 = nb_types.uint8
u64 = nb_types.uint64


@dataclass
class FAAOfflineData:
    """What the offline phase of the fast algebraic attack hands to the online phase.

    The offline phase records the annihilator's equation at each clock time,
    as coefficient rows over a monomial index, together with a linear relation
    on the low-degree multiple's output sequence. The online phase sums
    annihilator equations along that relation, weighted by keystream bits,
    which cancels the multiple's contribution, and solves.

    :ivar idx_to_comb: Monomial index to monomial (a sorted tuple of state bits).
    :ivar comb_to_idx: Monomial to monomial index.
    :ivar annihilator_equations: One row per clock time: the annihilator's
        coefficients at that time.
    :ivar linear_relation: The linear relation on the multiple's sequence,
        reversed for use as a dot product rather than a convolution.
    :ivar num_vars: The number of monomials indexed.
    :ivar keystream_needed: How many keystream bits the online phase needs.
    :ivar margin: The extra equations generated beyond the linear complexity.
    """
    idx_to_comb: dict[int, tuple[int, ...]]
    comb_to_idx: dict[tuple[int, ...], int]
    annihilator_equations: np.ndarray
    linear_relation: np.ndarray
    num_vars: int
    keystream_needed: int
    margin: int


@log.stage("Offline phase (Fast Algebraic Attack)")
def FAA_offline(
    feedback_fn: FeedbackFunction,
    annihilator: BooleanFunction,
    multiple: BooleanFunction,
    init_rounds: int,
    margin: int,
    time_limit: int,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles: list[MonomialProfile] | None = None,
    variable_blocks: list[list[int]] | None = None
) -> FAAOfflineData:

    start_time = time.time()

    # compute equations for the annihilator:
    if monomial_profiles != None and variable_blocks != None:
        log.info("Monomial profile optimization: on")

        log.step("Monomial profile for the annihilator")
        annihilator_mp = annihilator.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)

        log.step("Monomial profile for the low degree multiple")
        # Precompute LC for low degree multiple:
        multiple_mp = multiple.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)
        max_LC = multiple_mp.upper()

        log.step("Variable map")
        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, annihilator_mp, variable_blocks, complete_subsets = True
        )

        log.step("Linear relation")
        # use berlekamp_massey to get the exact relation
        feedback_fn.compile()
        multiple_compiled = multiple.compile()
        test_register = FeedbackRegister(random.randint(0,2**feedback_fn.size-1), feedback_fn)
        max_count = 1000*((2*max_LC+256)//1000 + 1)
        bits = log.progress("Bits processed", total=max_count, unit="bits")
        # berlekamp_massey_iterator always yields at least once, even for an
        # empty sequence, so the loop replaces these before anything reads them
        linear_complexity = 0
        linear_relation = np.array([], dtype='uint8')

        for curr_LC, curr_relation in berlekamp_massey_iterator(
            seq = (multiple_compiled(state._state) for state in test_register.run(2*max_LC+256)),
            yield_rate=1000
        ):
            bits.update(1000)
            if bits.shown:
                bits.set_status(f"Linear complexity: {curr_LC}/{max_LC}")

            linear_complexity = curr_LC
            linear_relation = curr_relation
        bits.close()

        # flip linear relation, due to dot product vs convolution
        linear_relation = linear_relation[::-1]
        margin += linear_complexity
        log.info("Linear complexity of the multiple: %d", linear_complexity)

        # use precomputed maps for faster eq generation and storage
        annihilator_eqs = EqStore(variable_indices)
        eq_gen = CubeEqGenerator(
            feedback_fn, annihilator, (len(variable_indices) + margin),
            variable_indices,
        )

        # No rank tracker on this path: the monomial profile already bounds the
        # variable count, so there is nothing to check ranks against. Its
        # absence *is* the flag -- see the guard in the equation loop below.
        annihilator_LU = None
        count_into_margin = 0
        total = len(variable_indices) + margin

    else:
        log.info("Monomial profile optimization: off")

        # create dynamic storage and generation
        annihilator_eqs = EqStore()
        annihilator_LU = LUEqStore()
        annihilator_LU.link(annihilator_eqs)

        # ensure all variables are in the eq store:
        for v in range(len(feedback_fn)):
            annihilator_eqs._update_known_monomials((v,))
            annihilator_LU._update_known_monomials((v,))

        eq_gen = SubstitutionEqGenerator(
            feedback_fn, annihilator, 2**feedback_fn.size
        )

        # have to check ranks, since number of
        # variables isnt known ahead of time
        count_into_margin = 0
        total = None

        # Precompute LC for low degree multiple:
        # because max_LC isnt known, test until there are no changes:
        feedback_fn.compile()
        test_register = FeedbackRegister(random.getrandbits(feedback_fn.size), feedback_fn)

        log.step("Linear relation")
        thousands = log.progress("Bits processed (thousands)")
        curr_LC = 0
        curr_relation = []
        # berlekamp_massey_iterator always yields at least once, even for an
        # empty sequence, so the loop replaces these before anything reads them
        linear_complexity = 0
        linear_relation = np.array([], dtype='uint8')
        for linear_complexity, linear_relation in berlekamp_massey_iterator(
            seq = (multiple.eval(state) for state in test_register.run(2**(feedback_fn.size))),
            yield_rate=1000
        ):
            if thousands.shown:
                thousands.set_status(f"Linear complexity: {curr_LC}")

            # check lengths first for more efficient short circuit:
            if (linear_complexity == curr_LC) and (len(linear_relation) == len(curr_relation)) and np.all(linear_relation == curr_relation):
                break

            thousands.update()
            curr_LC = linear_complexity
            curr_relation = linear_relation
        thousands.close()

        # flip linear relation, due to dot product vs convolution
        linear_relation = linear_relation[::-1]
        margin += linear_complexity
        log.info("Linear complexity of the multiple: %d", curr_LC)

    log.step("Generating equations")
    found = log.progress("Equations found", total=total, unit="eq")

    # main equation loop
    for t, ann_eq in enumerate(eq_gen): #type: ignore  (to narrow types correctly)
        ann_eq: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]

        # don't generate equations for initialization rounds
        if t < init_rounds: continue

        annihilator_eqs.insert_equation(ann_eq, identifier = t)
        found.update_to(annihilator_eqs.num_eqs)
        if found.shown and total is None:
            found.set_status(f"of {annihilator_eqs.num_vars} monomials so far, plus {margin} margin")

        if time_limit and (time.time() - start_time >= time_limit):
            found.close()
            log.warning("Time limit reached after %s, with %d equations",
                        format_duration(time.time() - start_time), annihilator_eqs.num_eqs)
            break

        # break step only necessary for dynamic stores:
        # reduces speed a fair bit, due to extra insert
        if annihilator_LU is not None:
            ann_independent = annihilator_LU.insert_equation(ann_eq, identifier = t)

            # continue for margin more steps after hitting linear
            # recurrent phase (not perfect but better than nothing)
            if not (ann_independent):
                count_into_margin += 1
                if count_into_margin == margin:
                    found.close()
                    log.info("Linear complexity reached; %d margin equations generated", margin)
                    break

    found.close()
    keystream_needed = max(annihilator_eqs.equation_ids.values()) + 1
    log.info("Keystream needed: %d bits; monomials: %d", keystream_needed, annihilator_eqs.num_vars)

    return FAAOfflineData(
        idx_to_comb = annihilator_eqs.idx_to_comb,
        comb_to_idx = annihilator_eqs.comb_to_idx,
        annihilator_equations = annihilator_eqs.equations[:annihilator_eqs.num_eqs,:annihilator_eqs.num_vars],
        linear_relation = linear_relation,
        num_vars = annihilator_eqs.num_vars,
        keystream_needed = keystream_needed,
        margin = margin - (linear_complexity),
    )


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
@log.stage("Online phase (Fast Algebraic Attack)")
def FAA_online(
    feedback_fn: FeedbackFunction,
    output_fn: BooleanFunction,
    keystream: list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]],
    attack_data: FAAOfflineData,
    test_length: int = 1000,
    time_limit: float | None = None,
    solver=None,
    online_store=None,
):
    if solver is None:
        solver = LUSolver()

    start_time = time.time()

    if type(keystream) == np.ndarray:
        pass
    elif type(keystream) == list:
        keystream = np.array(keystream, dtype = 'uint8')
    else:
        raise ValueError(f"Keystream must be a list or u8 ndarray, not {type(keystream)}")

    # unpack attack_data
    num_vars = attack_data.num_vars
    num_eqs = attack_data.keystream_needed

    annihilator_eqs = attack_data.annihilator_equations
    linear_relation = attack_data.linear_relation

    comb_to_idx = attack_data.comb_to_idx
    idx_to_comb = attack_data.idx_to_comb

    if online_store is None:
        online_store = LUEqStore(comb_to_idx, consistent=True)

    log.step("Equation substitution")

    from PyPR.Cryptanalysis.Components.Adapters.online_insertion import (
        make_online_inserter,
    )

    total_online_eqs = num_eqs - len(linear_relation)
    insert_eq, finalize = make_online_inserter(
        online_store, idx_to_comb,
        total_eqs=total_online_eqs, num_vars=num_vars,
    )

    for eq_idx in range(total_online_eqs):
        coef_vector = sum_over_linear_relationship(
            eq_idx, keystream, annihilator_eqs, linear_relation
        )
        if insert_eq(coef_vector, eq_idx):
            break
        if time_limit and (time.time() - start_time >= time_limit):
            log.warning("Time limit reached during substitution")
            break

    finalize()

    if isinstance(online_store, FilteringEqStore):
        solved_count = online_store.num_determined
    else:
        solved_count = online_store.num_eqs
    log.info("Variables solved by substitution: %d/%d", solved_count, num_vars)

    if time_limit and (time.time() - start_time >= time_limit):
        log.warning("Online phase timed out after %s", format_duration(time.time() - start_time))
        return None

    log.step("Solving")
    initial_state, _, _ = solver.solve(
        online_store, feedback_fn, output_fn, keystream,
        test_length=test_length,
    )
    return initial_state
