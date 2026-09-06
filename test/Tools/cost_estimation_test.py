"""Tests for CostEstimation: estimate_cost_comp and estimate_cost_cube.

These tests verify basic contractual properties of the cost estimators:
positivity, ordering relative to register complexity, determinism under fixed
seeds, and the expected regime differences between CMPR (composition favours
cube) and filter-generator (composition favours symbolic) architectures.

For full benchmark plots (run times can reach tens of minutes), see
experiments/equation_generation_profiling/cost_estimation/experiment.py.
"""
import random

import numpy as np

from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.BooleanLogic import VAR, AND, CONST
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template
from PyPR.Tools.CostEstimation import estimate_cost_comp, estimate_cost_cube
from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile, TermSet


M2  = MPR(2,  "3")
M3  = MPR(3,  "5")
M5  = MPR(5,  "12")
M7  = MPR(7,  "65")
M31 = MPR(31, "400002B8")


def test_comp_cost_positive():
    """estimate_cost_comp returns a positive value for a chained CMPR."""
    random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__()])
    C.generateChaining(template=arman_template(3))
    cost = estimate_cost_comp(C, VAR(0), C.monomial_profiles())
    assert cost > 0, f"Expected positive comp cost, got {cost}"


def test_cube_cost_positive():
    """estimate_cost_cube returns a positive value for a chained CMPR."""
    random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__()])
    C.generateChaining(template=arman_template(3))
    cost = estimate_cost_cube(C, VAR(0), C.monomial_profiles())
    assert cost > 0, f"Expected positive cube cost, got {cost}"


def test_comp_greater_than_cube_for_cmpr():
    """For a chained CMPR the composition cost should exceed the cube cost."""
    random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__()])
    C.generateChaining(template=arman_template(3))
    mps = C.monomial_profiles()
    comp = estimate_cost_comp(C, VAR(0), mps)
    cube = estimate_cost_cube(C, VAR(0), mps)
    assert comp > cube, f"Expected comp ({comp:,}) > cube ({cube:,}) for CMPR"


def test_unlinear_register_comp_less_than_cube():
    """For a plain MPR with a dense AND output the cube cost exceeds composition.

    The filter-generator regime: simple linear update, expensive output.
    """
    mps = [MonomialProfile([TermSet({0: M31.size}, {0: 1})]) for _ in range(M31.size)]
    output_fn = AND(*(VAR(i) for i in range(10)))
    comp = estimate_cost_comp(M31, output_fn, mps)
    cube = estimate_cost_cube(M31, output_fn, mps)
    assert cube > comp, f"Expected cube ({cube:,}) > comp ({comp:,}) for filter generator"


def test_zero_output_gives_zero_cube_cost():
    """CONST(0) output has no monomials, so cube cost is 0."""
    random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__()])
    C.generateChaining(template=arman_template())
    cost = estimate_cost_cube(C, CONST(0), C.monomial_profiles())
    assert cost == 0, f"Expected 0 cube cost for CONST(0), got {cost}"


def test_deterministic_with_seed():
    """Same seed gives identical estimates on repeated calls."""
    def run(seed):
        random.seed(seed)
        np.random.seed(seed)
        C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__()])
        C.generateChaining(template=arman_template(3))
        mps = C.monomial_profiles()
        return (
            estimate_cost_comp(C, VAR(0), mps),
            estimate_cost_cube(C, VAR(0), mps),
        )

    assert run(42) == run(42), "Estimates should be deterministic for fixed seed"
    assert run(1) != run(2), "Different seeds should (generally) give different estimates"


def test_larger_cmpr_higher_costs():
    """Adding more chained components increases both cost estimates."""
    random.seed(0)
    C_small = CMPR([M7.__copy__(), M3.__copy__()])
    C_small.generateChaining(template=arman_template())
    mps_s = C_small.monomial_profiles()
    comp_s = estimate_cost_comp(C_small, VAR(0), mps_s)
    cube_s = estimate_cost_cube(C_small, VAR(0), mps_s)

    random.seed(0)
    C_large = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__()])
    C_large.generateChaining(template=arman_template())
    mps_l = C_large.monomial_profiles()
    comp_l = estimate_cost_comp(C_large, VAR(0), mps_l)
    cube_l = estimate_cost_cube(C_large, VAR(0), mps_l)

    assert comp_l > comp_s, "Larger CMPR should have higher comp cost"
    assert cube_l > cube_s, "Larger CMPR should have higher cube cost"
