
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
from PyPR.Cryptanalysis.Components.EquationSolving.Grob_Solver import GrobnerSolver
from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

if TYPE_CHECKING:
    from PyPR.BooleanLogic import BooleanFunction

log = get_logger(__name__)


@dataclass
class NAAOfflineData:
    """What the offline phase of the naive algebraic attack hands to the online phase.

    The offline phase finds the clock times whose output equations are linearly
    independent, and stores their LU factorization; the online phase fills in
    the observed keystream bits at those times and solves.

    :ivar guess_vars: The monomials no equation solved for, as `(index, monomial)`
        pairs; the online solver guesses them.
    :ivar equation_times: For each solved monomial index, the clock time whose
        equation solved it.
    :ivar idx_to_comb: Monomial index to monomial (a sorted tuple of state bits).
    :ivar comb_to_idx: Monomial to monomial index.
    :ivar upper_matrix: The U factor over the monomials.
    :ivar lower_matrix: The L factor over the monomials.
    :ivar keystream_needed: How many keystream bits the online phase needs.
    """
    guess_vars: list[tuple[int, tuple[int, ...]]]
    equation_times: dict[int, int]
    idx_to_comb: dict[int, tuple[int, ...]]
    comb_to_idx: dict[tuple[int, ...], int]
    upper_matrix: np.ndarray
    lower_matrix: np.ndarray
    keystream_needed: int


@log.stage("Offline phase (Naive Algebraic Attack)")
def NAA_offline(
    feedback_fn, output_fn, init_rounds,
    time_limit,

    # both are needed to specify a monomial layout:
    # this enables optimizations
    monomial_profiles = None,
    variable_blocks = None
) -> NAAOfflineData:
    start_time = time.time()
    if monomial_profiles != None and variable_blocks != None:
        log.info("Monomial profile optimization: on")

        log.step("Monomial profile")
        # A map with all subsets filled in, to sum over cubes
        output_mp = output_fn.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ]).eval_ANF(monomial_profiles)

        log.step("Variable map")
        # A map with all subsets filled in, to sum over cubes
        variable_indices = get_var_map(
            feedback_fn, output_mp, variable_blocks, complete_subsets = True, include_constant=True
        )

        # use precomputed maps for faster eq generation and storage
        eqs = LUEqStore(variable_indices)

        eq_gen = CubeEqGenerator(
            feedback_fn, output_fn, 2**feedback_fn.size,
            variable_indices
        )
        total = len(variable_indices)

    else:
        log.info("Monomial profile optimization: off")

        eqs = LUEqStore()
        # ensure all variables are in the eq store:
        for v in range(len(feedback_fn)):
            eqs._update_known_monomials((v,))

        eq_gen = SubstitutionEqGenerator(feedback_fn, output_fn, 2**feedback_fn.size)
        # monomials are discovered as equations arrive, so there is no total
        total = None

    log.step("Generating equations")
    found = log.progress("Equations found", total=total, unit="eq")

    # main loop:
    for t, equation in enumerate(eq_gen): #type: ignore  (to narrow types correctly)
        equation: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]

        if t < init_rounds:
            continue

        linearly_independent = eqs.insert_equation(
            equation,
            identifier = t,
            translate_ANF = False
        )
        found.update_to(eqs.num_eqs)
        if found.shown and total is None:
            found.set_status(f"of {eqs.num_vars} monomials so far")

        if not linearly_independent:
            # all equations from this point are linearly dependent.
            found.close()
            log.info("Linear complexity reached at t=%d", t)
            break

        if time_limit and (time.time() - start_time >= time_limit):
            found.close()
            log.warning("Time limit reached after %s, with %d/%d equations",
                        format_duration(time.time() - start_time), eqs.num_eqs, eqs.num_vars)
            break

    found.close()
    not_solved = [(x,eqs.idx_to_comb[x]) for x in range(eqs.num_vars) if x not in eqs.equation_ids]
    keystream_needed = max(eqs.equation_ids.values()) + 1
    log.info("Keystream needed: %d bits; unsolved monomials: %d", keystream_needed, len(not_solved))

    return NAAOfflineData(
        guess_vars = not_solved,
        equation_times = eqs.equation_ids,
        idx_to_comb = eqs.idx_to_comb,
        comb_to_idx = eqs.comb_to_idx,
        upper_matrix = eqs.upper_matrix[:eqs.num_vars,:eqs.num_vars],
        lower_matrix = eqs.lower_matrix[:eqs.num_vars,:eqs.num_vars],
        keystream_needed = keystream_needed,
    )


@log.stage("Online phase (Naive Algebraic Attack)")
def NAA_online(
    feedback_fn, output_fn, keystream, attack_data: NAAOfflineData,
    test_length=1000, time_limit=None,
    solver=None,
):
    if isinstance(solver, GrobnerSolver):
        raise ValueError(
            "NAA with GrobnerSolver is not supported: NAA's offline phase produces "
            "LU matrices with separate constants, which Gröbner solving cannot use. "
            "Use LUSolver or GaussElimSolver instead."
        )

    start_time = time.time()

    # unpack attack_data
    var_map = attack_data.equation_times
    upper_matrix = attack_data.upper_matrix
    lower_matrix = attack_data.lower_matrix
    num_vars = len(upper_matrix)
    comb_to_idx = attack_data.comb_to_idx

    # reconstruct an LUEqStore from the stored matrices
    solved_store = LUEqStore(comb_to_idx)
    solved_store.upper_matrix[:num_vars,:num_vars] = upper_matrix
    solved_store.lower_matrix[:num_vars,:num_vars] = lower_matrix
    for v in range(num_vars):
        if v in var_map:
            solved_store.solved_for[v] = 1
    solved_store.num_eqs = sum(1 for v in range(num_vars) if v in var_map)
    solved_store.equation_ids = {v: var_map[v] for v in var_map}

    # build constant vector from keystream
    additional_constants = np.zeros([num_vars], dtype=np.uint8)
    for v in range(num_vars):
        if v in var_map:
            additional_constants[v] = keystream[var_map[v]]

    if solver is None:
        solver = LUSolver(additional_constants=additional_constants)
    else:
        solver.additional_constants = additional_constants

    if time_limit and (time.time() - start_time >= time_limit):
        log.warning("Time limit reached before solving")
        return None

    log.step("Solving")
    initial_state, _, _ = solver.solve(
        solved_store, feedback_fn, output_fn, keystream,
        test_length=test_length,
    )
    return initial_state
