
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from PyPR.Reporting import format_duration, get_logger

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

if TYPE_CHECKING:
    from PyPR.BooleanLogic import BooleanFunction

log = get_logger(__name__)


@dataclass
class RAAOfflineData:
    """What the offline phase of the reduced algebraic attack hands to the online phase.

    For each clock time t the offline phase records the annihilator equation
    `g_t` and the multiple equation `h_t`, as coefficient rows over a shared
    monomial index; the online phase combines them with the keystream bit
    `z_t` as `h_t + z_t g_t` and solves.

    :ivar idx_to_comb: Monomial index to monomial (a sorted tuple of state bits).
    :ivar comb_to_idx: Monomial to monomial index.
    :ivar annihilator_equations: One row per clock time: the annihilator's
        coefficients at that time.
    :ivar multiple_equations: One row per clock time: the multiple's
        coefficients at that time.
    :ivar num_vars: The number of monomials indexed.
    :ivar keystream_needed: How many keystream bits the online phase needs.
    :ivar margin: The extra equations generated beyond the linear complexity.
    """
    idx_to_comb: dict[int, tuple[int, ...]]
    comb_to_idx: dict[tuple[int, ...], int]
    annihilator_equations: np.ndarray
    multiple_equations: np.ndarray
    num_vars: int
    keystream_needed: int
    margin: int


@log.stage("Offline phase (Reduced Algebraic Attack)")
def RAA_offline(
    feedback_fn, annihilator, multiple,
    init_rounds, margin,
    time_limit,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles = None,
    variable_blocks = None
    ) -> RAAOfflineData:

    start_time = time.time()

    if monomial_profiles != None and variable_blocks != None:
        log.info("Monomial profile optimization: on")

        log.step("Monomial profile")
        selected = max((annihilator, multiple), key = lambda x: x.degree())
        selected_mp = selected.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)
        log.debug("Linear complexity bound from the monomial profile: %d", selected_mp.upper())

        log.step("Variable map")
        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, selected_mp, variable_blocks, complete_subsets = True
        )

        # use precomputed maps for faster eq generation and storage
        annihilator_eqs = EqStore(variable_indices)
        multiple_eqs = EqStore(variable_indices)

        # No rank trackers on this path: the monomial profile already bounds
        # the variable count. Their absence *is* the flag -- see the guard in
        # the equation loop below.
        annihilator_LU = None
        multiple_LU = None
        count_into_margin = 0

        eq_gen = CubeEqGenerator(
            feedback_fn, [annihilator, multiple], (len(variable_indices) + margin),
            variable_indices, #output_map = variable_indices
        )
        total = len(variable_indices) + margin

    else:
        log.info("Monomial profile optimization: off")

        # create dynamic storage and generation
        annihilator_eqs = EqStore()
        annihilator_LU = LUEqStore()
        multiple_eqs = EqStore()
        multiple_LU = LUEqStore()

        # link all eq stores
        multiple_eqs.link(annihilator_eqs)
        annihilator_eqs.link(multiple_eqs)
        annihilator_LU.link(annihilator_eqs)
        annihilator_eqs.link(annihilator_LU)
        multiple_LU.link(multiple_eqs)
        multiple_eqs.link(multiple_LU)

        # ensure all variables are in the eq store:
        for v in range(len(feedback_fn)):
            annihilator_eqs._update_known_monomials((v,))

        eq_gen = SubstitutionEqGenerator(
            feedback_fn, [annihilator, multiple], 2**feedback_fn.size
        )

        # have to check ranks, since number of
        # variables isnt known ahead of time
        count_into_margin = 0
        total = None

    log.step("Generating equations")
    found = log.progress("Equations found", total=total, unit="eq")

    # main equation loop
    for t, (ann_eq, mult_eq) in enumerate(eq_gen): # type: ignore (to narrow types correctly)
        ann_eq:  BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]
        mult_eq: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]

        # don't generate equations for initializatipon rounds
        if t < init_rounds: continue

        annihilator_eqs.insert_equation(ann_eq, identifier = t)
        multiple_eqs.insert_equation(mult_eq, identifier = t)
        found.update_to(multiple_eqs.num_eqs)
        if found.shown and total is None:
            found.set_status(f"of {multiple_eqs.num_vars} monomials so far, plus {margin} margin")

        if time_limit and (time.time() - start_time >= time_limit):
            found.close()
            log.warning("Time limit reached after %s, with %d equations",
                        format_duration(time.time() - start_time), multiple_eqs.num_eqs)
            break

        # break step only necessary for dynamic stores
        if annihilator_LU is not None and multiple_LU is not None:
            ann_indep = annihilator_LU.insert_equation(ann_eq, identifier = t)
            mult_indep = multiple_LU.insert_equation(mult_eq, identifier = t)

            # continue for margin more steps after both have hit their
            # linear recurrence phase (not perfect but better than nothing)
            if not (ann_indep or mult_indep):
                count_into_margin += 1
                if count_into_margin == margin:
                    found.close()
                    log.info("Linear complexity reached; %d margin equations generated", margin)
                    break

    found.close()
    keystream_needed = max(multiple_eqs.equation_ids.values()) + 1
    log.info("Keystream needed: %d bits; monomials: %d", keystream_needed, multiple_eqs.num_vars)

    return RAAOfflineData(
        idx_to_comb = multiple_eqs.idx_to_comb,
        comb_to_idx = multiple_eqs.comb_to_idx,
        annihilator_equations = annihilator_eqs.equations[:annihilator_eqs.num_eqs,:annihilator_eqs.num_vars],
        multiple_equations = multiple_eqs.equations[:multiple_eqs.num_eqs,:multiple_eqs.num_vars],
        num_vars = multiple_eqs.num_vars,
        keystream_needed = keystream_needed,
        margin = margin,
    )


# Dont need known bits: this is because each equation is cheap (relative to cube attacks)
# and the known bits doesnt /really/ help with the monomials (without a big loop), so it
# doesnt shrink the system that much, but does introduce a lot of overhead.
@log.stage("Online phase (Reduced Algebraic Attack)")
def RAA_online(
    feedback_fn, output_fn, keystream, attack_data: RAAOfflineData,
    test_length=1000, time_limit=None,
    solver=None, online_store=None,
):
    if solver is None:
        solver = LUSolver()

    start_time = time.time()

    # unpack attack_data
    num_vars = attack_data.num_vars
    num_eqs = attack_data.keystream_needed

    annihilator_eqs = attack_data.annihilator_equations
    multiple_eqs = attack_data.multiple_equations

    comb_to_idx = attack_data.comb_to_idx
    idx_to_comb = attack_data.idx_to_comb

    if online_store is None:
        online_store = LUEqStore(comb_to_idx, consistent=(() in comb_to_idx))

    log.step("Equation substitution")

    from PyPR.Cryptanalysis.Components.Adapters.online_insertion import (
        make_online_inserter,
    )

    insert_eq, finalize = make_online_inserter(
        online_store, idx_to_comb,
        total_eqs=num_eqs, num_vars=num_vars,
    )

    for eq_idx in range(num_eqs):
        coef_vector = np.zeros([num_vars], dtype="uint8")
        coef_vector ^= multiple_eqs[eq_idx]
        coef_vector ^= keystream[eq_idx] * annihilator_eqs[eq_idx]

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
