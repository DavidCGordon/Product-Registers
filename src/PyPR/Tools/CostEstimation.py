"""Asymptotic cost estimators for ANF equation generation.

Two methods for generating the ANF coefficients of a CMPR's output are compared:
  - ANF composition: symbolically compose the feedback function tree
  - Cube-based generation: compute coefficients via cube sums over register evaluations

Both functions return a count of *coefficient flips* — a hardware-agnostic unit that
lets you compare asymptotic scaling without the large constant-factor advantage of
the cube method (JIT, array-level XOR) obscuring the picture.

See Theory/Cube Equation Generation.md for the mathematical background.
"""
from __future__ import annotations

import math
from typing import TYPE_CHECKING

from PyPR.BooleanLogic.Gates import AND, OR, NAND, NOR, XOR, XNOR
from PyPR.Tools.RootCounting.MonomialProfile import MonomialProfile, TermSet

if TYPE_CHECKING:
    from PyPR.BooleanLogic import BooleanFunction
    from PyPR.FeedbackFunctions import FeedbackFunction


def estimate_cost_comp(
    feedback_fn: "FeedbackFunction",
    output_fn: "BooleanFunction",
    monomial_profiles: list[MonomialProfile],
) -> int:
    """Estimate the ANF composition cost in coefficient flips.

    Walks the binarized function tree of every feedback function and the output
    function, computing the monomial profile and cumulative cost bottom-up via
    dynamic programming. Multiplication cost is ``(|M1| * |M2|) / 2`` and addition
    cost is ``min(|M1|, |M2|)``. These assume O(1) coefficient lookup, which is
    favorable to composition (see Theory doc).

    :param feedback_fn: The register's feedback function.
    :type feedback_fn: FeedbackFunction
    :param output_fn: The observed output function (e.g. ``VAR(0)``).
    :type output_fn: BooleanFunction
    :param monomial_profiles: Per-bit monomial profiles for the register state,
        as returned by ``CMPR.monomial_profiles()``.
    :type monomial_profiles: list[MonomialProfile]
    :return: Estimated total coefficient flips across all functions.
    :rtype: int
    """
    mp_cache: dict = {}
    cost_cache: dict = {}
    total = 0

    for fn in list(feedback_fn.fn_list) + [output_fn]:
        new_fn = fn.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one()),
        ]).binarize()

        stack = [new_fn]
        last = None

        while stack:
            node = stack[-1]

            if node is False:
                last = stack.pop()
                continue

            if node in cost_cache:
                last = stack.pop()
                continue

            if last is False:
                mp_cache[node] = node._eval_ANF(mp_cache, monomial_profiles)
                left_cost = cost_cache[node.args[0]]
                right_cost = cost_cache[node.args[1]]
                cost_cache[node] = left_cost + right_cost

                left_upper = mp_cache[node.args[0]].upper()
                right_upper = mp_cache[node.args[1]].upper()

                if type(node) in (AND, OR, NAND, NOR):
                    cost_cache[node] += (left_upper * right_upper) // 2
                elif type(node) in (XOR, XNOR):
                    cost_cache[node] += min(left_upper, right_upper)

                last = stack.pop()
                continue

            if node.is_leaf():
                mp_cache[node] = node._eval_ANF(mp_cache, monomial_profiles)
                cost_cache[node] = 0
                last = stack.pop()
                continue

            stack.append(False)
            for child in reversed(node.args):
                stack.append(child)

        total += cost_cache[new_fn]

    return total


def estimate_cost_cube(
    feedback_fn: "FeedbackFunction",
    output_fn: "BooleanFunction",
    monomial_profiles: list[MonomialProfile],
) -> int:
    """Estimate the cube-based equation generation cost in coefficient flips.

    For each output monomial of degree *d*, the graded-split optimisation reduces
    the per-monomial cost from O(2^d) to::

        (2 ** (d // 2)) * (2 + d % 2) - 1

    The total is summed over all monomials in the output's monomial profile,
    using the degree-rollover iteration from ``MonomialProfile.get_monomials()``.

    :param feedback_fn: The register's feedback function (used only to extract
        the block structure via ``monomial_profiles``).
    :type feedback_fn: FeedbackFunction
    :param output_fn: The observed output function (e.g. ``VAR(0)``).
    :type output_fn: BooleanFunction
    :param monomial_profiles: Per-bit monomial profiles for the register state.
    :type monomial_profiles: list[MonomialProfile]
    :return: Estimated total coefficient flips.
    :rtype: int
    """
    mp = output_fn.remap_constants([
        (0, MonomialProfile.logical_zero()),
        (1, MonomialProfile.logical_one()),
    ]).eval_ANF(monomial_profiles)

    def _rect(term: TermSet) -> tuple[int, ...]:
        max_key = max(list(term.totals.keys()) + [0])
        row = [0] * (max_key + 1)
        for k, c in term.counts.items():
            row[k] = c
        return tuple(row)

    def _safe_get(row: tuple[int, ...], d: int) -> int:
        return row[d] if d < len(row) else 0

    rects = [_rect(term) for term in mp.terms]
    if not rects:
        return 0

    dim = max(len(r) for r in rects)

    totals = [0] * dim
    for term in mp.terms:
        for k, v in term.totals.items():
            totals[k] = v

    # stable sort: primary ascending on each dimension > 0, then descending on 0
    for d in range(1, dim):
        rects = sorted(rects, key=lambda r: _safe_get(r, d))
    rects = sorted(rects, key=lambda r: _safe_get(r, 0), reverse=True)

    estimate = 0
    rect_idx = 0
    curr_vec = [0] * dim

    while rect_idx < len(rects):
        curr_vec[0] += 1

        rollover_idx = 0
        rollover = list(curr_vec)

        while any(_safe_get(rects[rect_idx], i) < rollover[i] for i in range(dim)):
            if rollover_idx < dim - 1:
                rollover[rollover_idx] = 0
                rollover[rollover_idx + 1] += 1
                rollover_idx += 1
            else:
                rollover = list(curr_vec)
                rollover_idx = 0
                rect_idx += 1
                if rect_idx >= len(rects):
                    break

        curr_vec = rollover
        if rect_idx >= len(rects):
            break

        degree = sum(curr_vec)
        count = 1
        for i in range(dim):
            count *= math.comb(totals[i], curr_vec[i])

        cost = (2 ** (degree // 2)) * (2 + degree % 2) - 1
        estimate += count * cost

    return estimate
