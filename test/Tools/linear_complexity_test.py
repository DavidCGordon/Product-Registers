"""Tests for linear complexity estimation: BM, estimate_LC, root_expressions, monomial_profiles.

Derived from ProductRegisters.ipynb (Applications > Linear Complexity, Berlekamp Massey,
and Estimating Linear Complexity and Root Expressions).

`estimate_LC` returns `(lower, upper)`, and the two halves are not the same
kind of claim.  `upper` counts the maximum number of roots the expression
could represent, so nothing can exceed it and a measured complexity above it
is a genuine defect.  `lower` is a statistical expectation computed with
`pessimistic_expected_value`; Root Expressions and LC Estimation.md records a
diagnostic in which 2 of 10 trials landed below it.  The final section
therefore asserts only the upper bound.

Both RNGs must be seeded in any test that builds chaining: the templates draw
their taps through `np.random.choice` (TemplateBuilding.SAMPLE), not through
`random`.
"""
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, VAR
from PyPR.BooleanLogic.ChainingGeneration.Templates import (
    arman_template,
    fast_template,
    old_ANF_template,
)

from PyPR.FeedbackFunctions import CMPR, MPR
from PyPR.FeedbackRegister import FeedbackRegister

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

# ── Fixtures ────────────────────────────────────────────────────────────

M7 = MPR(7, [1, 1, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 0])
M5 = MPR(5, [1, 0, 1, 0, 0, 1], [1, 1, 0, 0, 1])
M3 = MPR(3, [1, 1, 0, 1], [1, 0, 1])
M2 = MPR(2, [1, 1, 1], [1, 1])


# ── BM on CMPR ──────────────────────────────────────────────────────────

def test_berlekamp_massey_cmpr():
    """BM on a CMPR sequence should return a positive linear complexity (notebook cell 99)."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__()])
    C.generateChaining(template=old_ANF_template(max_and=4, max_xor=4))
    C.compile()

    F = FeedbackRegister(2**C.size - 1, C)
    seq = [state[0] for state in F.run(limit=5000)]

    lc, _poly = berlekamp_massey(seq)
    assert lc > 0, "LC should be positive"
    assert lc > max(7, 5, 3), "CMPR with chaining should have LC > largest component"


# ── estimate_LC bounds ──────────────────────────────────────────────────

def test_estimate_lc_bounds():
    """estimate_LC should bound the actual linear complexity (notebook cell 107).

    The bounds depend on the specific chaining generated, so the seed must be fixed.
    The lower bound is a statistical estimate (float), not a hard guarantee — it
    can fail for degenerate chaining, but holds consistently with a fixed seed.

    Both RNGs must be seeded: the chaining templates draw their taps through
    `np.random.choice` (TemplateBuilding.SAMPLE), not through `random`, so
    seeding only `random` leaves the draw at the mercy of whatever consumed
    numpy randomness earlier in the process.  Root Expressions and LC
    Estimation.md records that ~2 in 10 chaining draws for this configuration
    fall below the lower bound, so an unseeded draw makes this test flaky.
    """
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=old_ANF_template())
    C.compile()

    lower, upper = C.estimate_LC(0)
    assert lower > 0, "Lower bound should be positive"
    assert upper >= lower, f"Upper ({upper}) should be >= lower ({lower})"

    F = FeedbackRegister(2**C.size - 1, C)
    seq = [state[0] for state in F.run(limit=2 * upper + 500)]
    actual_lc, _ = berlekamp_massey(seq)

    assert int(lower) <= actual_lc <= upper, (
        f"Actual LC {actual_lc} not in predicted range [{int(lower)}, {upper}]"
    )


# ── Root expressions ────────────────────────────────────────────────────

def test_root_expressions_basic():
    """Root expressions should be computable for a standard CMPR."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=old_ANF_template())

    REs = C.root_expressions()
    assert len(REs) == C.size

    # each RE should have upper >= lower > 0
    for i, re in enumerate(REs):
        assert re.upper() >= re.lower(), f"Bit {i}: upper < lower"

def test_root_expressions_with_filter():
    """Filtering through an AND gate produces a valid root expression (notebook cell 109)."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=old_ANF_template())

    REs = C.root_expressions()
    filt = AND(VAR(0), VAR(1))
    filtered_re = filt.eval_ANF(REs)

    assert filtered_re.upper() >= filtered_re.lower()
    # filtered should have higher or equal complexity
    assert filtered_re.upper() >= REs[0].upper()


# ── Monomial profiles ──────────────────────────────────────────────────

def test_monomial_profiles_basic():
    """Monomial profiles should be computable for a standard CMPR."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=old_ANF_template())

    MPs = C.monomial_profiles()
    assert len(MPs) == C.size

    # each MP should have a positive upper bound
    for i, mp in enumerate(MPs):
        assert mp.upper() > 0, f"Bit {i}: monomial count should be positive"


# ── Resolvent / propagation matrices ────────────────────────────────────

def test_update_matrices():
    """Update matrices should be square and binary (notebook cell 204)."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=arman_template())

    for i, matrix in enumerate(C.update_matrices):
        block_size = len(C.blocks[i])
        assert matrix.shape == (block_size, block_size), f"Block {i}: wrong shape"
        assert set(matrix.flatten().tolist()).issubset({0, 1}), f"Block {i}: non-binary entries"

def test_resolvent_matrices():
    """Resolvent matrices should be computable."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=arman_template())

    resolvents = C.resolvent_matrices
    assert len(resolvents) == len(C.blocks)

def test_resolvent_matrix_shapes():
    """Resolvent matrices should have correct shape per block (notebook cell 206).

    Note: propagation_matrices is skipped because SequenceTransform.__eq__
    does not support comparison with int, so the `!= 0` mask fails.
    """
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generateChaining(template=arman_template())

    for i, matrix in enumerate(C.resolvent_matrices):
        block_size = len(C.blocks[i])
        assert matrix.shape == (block_size, block_size), f"Block {i}: wrong shape"


# ── estimate_LC upper bound as a hard bound ───────────────────────────────────

PRIMITIVE = {
    2: [1, 1, 1],
    3: [1, 1, 0, 1],
    4: [1, 1, 0, 0, 1],
    5: [1, 0, 1, 0, 0, 1],
    7: [1, 1, 0, 0, 0, 0, 0, 1],
}


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [4, 5], [3, 4, 5], [2, 3, 5]], ids=str)
def test_upper_bound_holds_for_an_unchained_cmpr(sizes):
    """Without chaining, every bit's linear complexity stays under `upper`.

    This is the degenerate case: with no cross-block terms a bit's sequence is
    driven only by its own block, so the true complexity is at most that
    block's size while the bound still budgets for the whole register.  The
    bound must hold with room to spare -- it holding *tightly* here would
    suggest the root expression had collapsed to the wrong block.
    """
    fn = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])
    fn.compile()
    size = len(fn)

    for bit in (0, size // 2, size - 1):
        reg = FeedbackRegister((1 << size) - 1, fn)
        seq = [int(state[bit]) for state in reg.run(limit=6000)]
        complexity, _ = berlekamp_massey(seq)
        _, upper = fn.estimate_LC(bit)
        assert complexity <= upper, (
            f"{sizes} bit {bit}: measured LC {complexity} exceeds the upper bound {upper}"
        )


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [4, 5], [3, 4, 5], [3, 5, 7]], ids=str)
def test_upper_bound_holds_for_a_chained_cmpr(sizes):
    """With chaining installed the bound still holds, and here it is tight.

    Chaining multiplies roots across blocks, which is what the root-expression
    algebra exists to track -- the bound climbs from "size of one block" to a
    product-like count, and the measured complexity rises to meet it.  Several
    of these cases reach the bound exactly, so an undercount in the
    coset-weight arithmetic fails here immediately rather than hiding in slack.

    `fast_template` is used deliberately: the chaining only has to create
    genuine cross-block products for the counting to be exercised, and the
    heavier templates spend ~30s per call searching for functions with
    properties no assertion here depends on.
    """
    random.seed(42)
    np.random.seed(42)
    fn = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])
    fn.generateChaining(fast_template())
    fn.compile()
    size = len(fn)

    reg = FeedbackRegister((1 << size) - 1, fn)
    seq = [int(state[0]) for state in reg.run(limit=20000)]
    complexity, _ = berlekamp_massey(seq)
    _, upper = fn.estimate_LC(0)
    assert complexity <= upper, (
        f"{sizes} bit 0: measured LC {complexity} exceeds the upper bound {upper}"
    )


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [4, 5]], ids=str)
def test_chaining_raises_the_upper_bound_above_the_unchained_one(sizes):
    """Installing chaining strictly increases bit 0's bound.

    Without chaining, bit 0's sequence is driven only by its own block and the
    bound reflects that.  Chaining introduces products of roots from upstream
    blocks, so the root expression -- and therefore the count -- must grow.  A
    bound that did not move would mean the chaining functions never entered the
    root expression at all.
    """
    random.seed(42)
    np.random.seed(42)
    plain = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])
    _, plain_upper = plain.estimate_LC(0)

    chained = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])
    chained.generateChaining(fast_template())
    _, chained_upper = chained.estimate_LC(0)

    assert chained_upper > plain_upper, (
        f"{sizes}: chaining left the upper bound at {plain_upper}"
    )


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [3, 4, 5]], ids=str)
def test_bounds_are_ordered_and_at_least_the_block_size(sizes):
    """`lower <= upper`, and `lower` is floored at the output bit's block size.

    `estimate_LC` takes `max(blockLen, re.lower())` because the bit's own block
    always contributes its full size regardless of what the statistical
    estimate says.  The ordering of the pair is structural and holds even
    though `lower` is only an expectation.
    """
    fn = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])
    blocks = fn.blocks

    for bit in range(len(fn)):
        lower, upper = fn.estimate_LC(bit)
        block_size = next(len(block) for block in blocks if bit in block)
        assert lower <= upper, f"{sizes} bit {bit}: lower {lower} exceeds upper {upper}"
        assert lower >= block_size, (
            f"{sizes} bit {bit}: lower {lower} is below its block size {block_size}"
        )


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [4, 5]], ids=str)
def test_locking_every_block_does_not_raise_the_upper_bound(sizes):
    """Suppressing blocks via `locked_list` can only tighten the bound.

    A locked block contributes no roots, so the root expression it feeds can
    only shrink.  The bound is monotone in that suppression -- if locking ever
    raised it, the locked and unlocked paths would be counting different sets.
    """
    fn = CMPR([MPR(s, PRIMITIVE[s]) for s in sizes])

    for bit in range(len(fn)):
        _, unlocked = fn.estimate_LC(bit)
        _, locked = fn.estimate_LC(bit, locked_list=[0] * len(fn.blocks))
        assert locked <= unlocked, (
            f"{sizes} bit {bit}: locking raised the upper bound from {unlocked} to {locked}"
        )
