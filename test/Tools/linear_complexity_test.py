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
import itertools
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, VAR, XOR
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
    C.generate_chaining(template=old_ANF_template(max_and=4, max_xor=4))
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
    C.generate_chaining(template=old_ANF_template())
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
    C.generate_chaining(template=old_ANF_template())

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
    C.generate_chaining(template=old_ANF_template())

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
    C.generate_chaining(template=old_ANF_template())

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
    C.generate_chaining(template=arman_template())

    for i, matrix in enumerate(C.update_matrices):
        block_size = len(C.blocks[i])
        assert matrix.shape == (block_size, block_size), f"Block {i}: wrong shape"
        assert set(matrix.flatten().tolist()).issubset({0, 1}), f"Block {i}: non-binary entries"

def test_resolvent_matrices():
    """Resolvent matrices should be computable."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generate_chaining(template=arman_template())

    resolvents = C.resolvent_matrices
    assert len(resolvents) == len(C.blocks)

def test_propagation_matrices_default_update_all_ones():
    """With the default update every bit of a block reaches every other bit.

    The resolvent is (I + UD)^{-1} = sum_t U^t D^t, and U is multiplication
    by x on F_2[x]/P. x is invertible of finite order N, so U^t e_j = x^{t+j}
    reaches every x^i (docs/theory/Resolvent Analysis.md, Proposition).
    """
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generate_chaining(template=arman_template())

    for i, mask in enumerate(C.propagation_matrices):
        block_size = len(C.blocks[i])
        assert mask.shape == (block_size, block_size)
        assert mask.all(), f"Block {i}: a zero in the propagation matrix"


def test_propagation_matrices_other_update_polynomials():
    """A non-default update need not reach every bit.

    Entry (i, j) is 1 exactly when x^i appears in u^t x^j for some t >= 0.
    With u = 1 the update is the identity, so only the diagonal is reached.
    With P = x^4 + x + 1 and u = x + x^3 three entries stay zero.
    """
    identity = CMPR([MPR(4, "13", [1, 0, 0, 0])]).propagation_matrices[0]
    assert (identity == np.eye(4, dtype=np.uint8)).all()

    partial = CMPR([MPR(4, "13", [0, 1, 0, 1])]).propagation_matrices[0]
    assert int(partial.sum()) == 13


def test_resolvent_matrix_shapes():
    """Resolvent matrices should have correct shape per block (notebook cell 206)."""
    random.seed(0)
    np.random.seed(0)
    C = CMPR([M7.__copy__(), M5.__copy__(), M3.__copy__(), M2.__copy__()])
    C.generate_chaining(template=arman_template())

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
    fn.generate_chaining(fast_template())
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
    chained.generate_chaining(fast_template())
    _, chained_upper = chained.estimate_LC(0)

    assert chained_upper > plain_upper, (
        f"{sizes}: chaining left the upper bound at {plain_upper}"
    )


@pytest.mark.parametrize("sizes", [[3, 4], [3, 5], [3, 4, 5]], ids=str)
def test_bounds_are_ordered_and_at_least_the_block_size(sizes):
    """`lower <= upper`, and `lower` is floored at the output bit's block size.

    `estimate_LC` takes `max(block_len, re.lower())` because the bit's own block
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


@pytest.mark.parametrize("sizes", [(5, 3), (7, 4, 3), (7, 5, 3, 2)])
def test_mesh_path_agrees_with_the_default_path(sizes):
    """The mesh optimization is a faster route to the same root expressions.

    Both bounds are compared on every bit. Two things the mesh path has to
    get right, and once did not, are covered by the seeds below: a constant
    term in a block's own chaining adds the root 1, and a constant in an
    upstream block is inherited, because any product of upstream terms
    contains 1 * 1 = 1.
    """
    polys = {7: "83", 5: "25", 4: "13", 3: "B", 2: "7"}
    for max_and in (1, 2, 3):
        for seed in range(5):
            random.seed(seed)
            np.random.seed(seed)
            fn = CMPR([MPR(n, polys[n]) for n in sizes])
            fn.generate_chaining(template=arman_template(max_and=max_and))
            mesh = fn.root_expressions()
            default = fn.root_expressions(force_default=True)
            for bit in range(fn.size):
                assert (mesh[bit].upper(), mesh[bit].lower()) == (default[bit].upper(), default[bit].lower()), (
                    f"{sizes} max_and={max_and} seed={seed} bit {bit}"
                )


def test_mesh_path_agrees_with_the_default_path_under_locking():
    """The two paths also agree for every pattern of locked blocks.

    A locked block (falsy flag) adds no roots of its own. The mesh path once
    wrote the source block's roots regardless of its flag, and turned a
    locked block into a field of size 0.
    """
    polys = {7: "83", 4: "13", 3: "B"}
    for max_and in (1, 3):
        random.seed(0)
        np.random.seed(0)
        fn = CMPR([MPR(n, polys[n]) for n in (7, 4, 3)])
        fn.generate_chaining(template=arman_template(max_and=max_and))
        for flags in itertools.product([0, 1], repeat=3):
            mesh = fn.root_expressions(list(flags))
            default = fn.root_expressions(list(flags), force_default=True)
            for bit in range(fn.size):
                assert (mesh[bit].upper(), mesh[bit].lower()) == (default[bit].upper(), default[bit].lower()), (
                    f"max_and={max_and} flags={flags} bit {bit}"
                )


def test_locking_a_block_removes_its_size_as_the_lower_floor():
    """estimate_LC floors the lower bound at the bit's block size only if that block is unlocked."""
    random.seed(0)
    np.random.seed(0)
    fn = CMPR([MPR(7, "83"), MPR(4, "13"), MPR(3, "B")])
    fn.generate_chaining(template=arman_template(max_and=1))
    source_bit = fn.blocks[0][0]
    unlocked_lower, _ = fn.estimate_LC(source_bit)
    locked_lower, locked_upper = fn.estimate_LC(source_bit, locked_list=[0, 1, 1])
    # the source block takes no chaining input, so locked it contributes nothing
    assert unlocked_lower >= 7
    assert locked_upper == 0
    assert locked_lower == 0


def test_mesh_path_requires_pairwise_coprime_sizes(monkeypatch):
    """Sizes 6 and 4 differ but share a factor, so root_expressions takes the default path.

    The mesh multiplies root counts across fields, which is exact only when
    the fields meet in F_2 (Root Expressions, Proposition 3); F_64 and F_16
    share F_4.
    """
    import PyPR.Tools.RootCounting.MeshOptimization as mesh_optimization

    calls = []
    original = mesh_optimization.re_compute_single_mesh
    monkeypatch.setattr(
        mesh_optimization, "re_compute_single_mesh",
        lambda *args: calls.append(args) or original(*args),
    )
    for sizes, polys in [((6, 4), {6: "43", 4: "13"}), ((5, 4), {5: "25", 4: "13"})]:
        random.seed(0)
        np.random.seed(0)
        fn = CMPR([MPR(n, polys[n]) for n in sizes])
        fn.generate_chaining(template=arman_template(max_and=2))
        calls.clear()
        fn.root_expressions()
        used_mesh = len(calls) > 0
        assert used_mesh == (sizes == (5, 4)), f"{sizes}: mesh used = {used_mesh}"


def test_mesh_path_on_a_block_without_chaining_input():
    """A block that takes no input from upstream keeps only its own roots and monomials.

    The mesh reads a block's degree as the rate for trading weight upstream,
    and that degree counts the block's own linear term, so a block with no
    upstream input would otherwise be treated as chaining linearly.
    """
    fn = CMPR([MPR(7, "83"), MPR(5, "25"), MPR(3, "B")])  # no chaining generated
    for block in fn.blocks:
        bit = block[0]
        assert fn.root_expressions()[bit].upper() == len(block)
        assert fn.monomial_profiles()[bit].upper() == len(block)


def test_monomial_profile_mesh_path_agrees_with_the_default_path():
    """Mesh and default monomial profiles give the same counts."""
    polys = {7: "83", 5: "25", 4: "13", 3: "B", 2: "7"}
    for sizes in [(5, 3), (7, 4, 3), (7, 5, 3, 2)]:
        for max_and in (1, 3):
            random.seed(0)
            np.random.seed(0)
            fn = CMPR([MPR(n, polys[n]) for n in sizes])
            fn.generate_chaining(template=arman_template(max_and=max_and))
            mesh = fn.monomial_profiles()
            default = fn.monomial_profiles(force_default=True)
            for bit in range(fn.size):
                assert mesh[bit].upper() == default[bit].upper(), f"{sizes} max_and={max_and} bit {bit}"


def test_monomial_profile_mesh_path_with_repeated_and_one_bit_blocks():
    """Monomials are indexed by register, so the mesh needs no condition on block sizes.

    Repeated sizes and one-bit blocks run through the mesh and give the same
    counts as the default path.
    """
    for components in [
        [MPR(4, "13"), MPR(4, "19")],
        [MPR(5, "25"), MPR(5, "29"), MPR(3, "B")],
        [MPR(5, "25"), MPR(1, "3"), MPR(3, "B")],
    ]:
        for seed in range(2):
            random.seed(seed)
            np.random.seed(seed)
            fn = CMPR([component.__copy__() for component in components])
            fn.generate_chaining(template=arman_template(max_and=2))
            mesh = fn.monomial_profiles()
            default = fn.monomial_profiles(force_default=True)
            for bit in range(fn.size):
                assert mesh[bit].upper() == default[bit].upper(), f"{[len(b) for b in fn.blocks]} seed={seed} bit {bit}"


def test_mesh_path_locking_without_chaining_constants():
    """With no constant in any chaining function, a locked block must not supply the root 1.

    The templates often add a constant, which hides the case; this register's
    chaining is written by hand without one.
    """
    fn = CMPR([MPR(7, "83"), MPR(4, "13"), MPR(3, "B")])
    source, middle, last = fn.blocks
    for i, bit in enumerate(middle):
        fn.fn_list[bit] = XOR(fn.fn_list[bit], AND(VAR(source[i]), VAR(source[i + 1])))
    for i, bit in enumerate(last):
        fn.fn_list[bit] = XOR(fn.fn_list[bit], VAR(middle[i]))
    for flags in itertools.product([0, 1], repeat=3):
        mesh = fn.root_expressions(list(flags))
        default = fn.root_expressions(list(flags), force_default=True)
        for bit in range(fn.size):
            assert (mesh[bit].upper(), mesh[bit].lower()) == (default[bit].upper(), default[bit].lower()), (
                f"flags={flags} bit {bit}"
            )


def test_chaining_nonlinear_in_own_bits_is_rejected():
    """A block whose update multiplies its own bits is outside the model, so the estimators refuse it.

    On such a (7, 4, 3) register Berlekamp-Massey measured 90 and 81 where
    the mesh gave 32 and 35 and the default path 4 and 7: neither is a bound.
    """
    fn = CMPR([MPR(7, "83"), MPR(4, "13"), MPR(3, "B")])
    source, middle, last = fn.blocks
    for i, bit in enumerate(middle):
        fn.fn_list[bit] = XOR(fn.fn_list[bit], AND(VAR(middle[(i + 1) % len(middle)]), VAR(source[i])))
    for i, bit in enumerate(last):
        fn.fn_list[bit] = XOR(fn.fn_list[bit], VAR(middle[i]))
    with pytest.raises(ValueError, match="multiplies its own bits"):
        fn.root_expressions()
    with pytest.raises(ValueError, match="multiplies its own bits"):
        fn.monomial_profiles()
    with pytest.raises(ValueError, match="multiplies its own bits"):
        fn.estimate_LC(0)

