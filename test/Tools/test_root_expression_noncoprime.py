"""Root expressions over fields whose degrees share a factor.

`RootExpression.upper()` counts a product of coset classes from different fields
as the product of their sizes (the box), which is exact only when the degrees
are pairwise coprime (Proposition 3 of docs/theory/Root Expressions and LC
Estimation.md). When they share a factor, `upper()` and `lower()` cap such a
product at the subgroup the fields generate, and `upper()` logs a warning.
These tests check that the cap keeps `upper()` an upper bound, against brute
force and against Berlekamp-Massey on CMPRs with non-coprime block sizes.
"""
import itertools
import math
import random

import numpy as np
import pytest

from PyPR.BooleanLogic.ChainingGeneration.Templates import (
    arman_template,
    fast_template,
    old_ANF_template,
)

from PyPR.FeedbackFunctions import CMPR, MPR
from PyPR.FeedbackRegister import FeedbackRegister

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey
from PyPR.Tools.RootCounting import RootExpression as root_expression_module
from PyPR.Tools.RootCounting.Combinatorics import binsum
from PyPR.Tools.RootCounting.JordanSet import JordanSet
from PyPR.Tools.RootCounting.RootExpression import RootExpression


def test_coprime_product_is_counted_as_a_box_without_warning(monkeypatch):
    # degrees 3 and 5 are coprime: distinct pairs give distinct products
    # (Proposition 3), so the count is the product of the class sizes.
    # Warnings are recorded from the module's logger: PyPR prints through its
    # own console, which capsys does not reliably see.
    warnings = []
    monkeypatch.setattr(root_expression_module.log, "warning", lambda msg, *args: warnings.append(msg % args))
    re = RootExpression({(3, 5): {JordanSet({3: 1, 5: 2}, {1})}})
    assert re.upper() == binsum(3, 1) * binsum(5, 2)
    assert warnings == []


def test_documented_non_coprime_example_is_capped_to_its_true_size(monkeypatch):
    # The example after Proposition 3: <2.1><4.2> has 2 * 10 = 20 tuples but
    # only 15 distinct products. GF(4)* and GF(16)* generate a subgroup of
    # order lcm(3, 15) = 15, which caps the box at the true size.
    warnings = []
    monkeypatch.setattr(root_expression_module.log, "warning", lambda msg, *args: warnings.append(msg % args))
    re = RootExpression({(2, 4): {JordanSet({2: 1, 4: 2}, {1})}})
    assert re.upper() == 15
    assert len(warnings) == 1
    assert "(2, 4)" in warnings[0]
    assert re.lower() <= re.upper()


@pytest.mark.parametrize("degrees", [(2, 4), (2, 6), (3, 6), (4, 6), (2, 4, 3), (3, 6, 5), (2, 3, 6)])
def test_non_coprime_upper_bound_is_never_below_the_true_count(degrees):
    # Brute force inside GF(2^L), L = lcm of the degrees, with exponents taken
    # relative to a primitive element. GF(2^e) sits there as the exponents
    # divisible by (2^L - 1)/(2^e - 1). A field's class <e.w> is defined
    # relative to its own generator beta, here a random generator of GF(2^e)*
    # (scale * k with k a unit), not necessarily the compatible one: a CMPR's
    # blocks choose their polynomials independently. Its elements are beta^x
    # with x an e-bit string of weight 1..w. upper() must count at least the
    # distinct products, and at most the plain box.
    random.seed(11)
    L = math.lcm(*degrees)
    mod = 2**L - 1
    for _ in range(6):
        weights = {e: random.randint(1, e) for e in degrees}
        classes = []
        for e in degrees:
            order = 2**e - 1
            unit = random.choice([k for k in range(1, order + 1) if math.gcd(k, order) == 1])
            scale = (mod // order) * unit
            classes.append([
                (x * scale) % mod
                for x in range(1, 2**e)
                if 1 <= x.bit_count() <= weights[e]
            ])
        products = {sum(t) % mod for t in itertools.product(*classes)}

        re = RootExpression({tuple(sorted(degrees)): {JordanSet(dict(weights), {1})}})
        box = math.prod(binsum(e, weights[e]) for e in degrees)
        assert len(products) <= re.upper() <= box, (degrees, weights)
        assert re.lower() <= re.upper()


# ── Ground truth: Berlekamp-Massey on CMPRs with non-coprime block sizes ──────

# Primitive polynomials, as coefficient lists with index i holding x^i.
PRIMITIVE = {
    2: [[1, 1, 1]],
    3: [[1, 1, 0, 1], [1, 0, 1, 1]],
    4: [[1, 1, 0, 0, 1], [1, 0, 0, 1, 1]],
    5: [[1, 0, 1, 0, 0, 1], [1, 0, 0, 1, 0, 1]],
    6: [[1, 1, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 1, 1], [1, 1, 0, 1, 1, 0, 1]],
    8: [[1, 0, 1, 1, 1, 0, 0, 0, 1], [1, 0, 1, 1, 0, 1, 0, 0, 1]],
    9: [[1, 0, 0, 0, 1, 0, 0, 0, 0, 1], [1, 0, 0, 1, 0, 0, 0, 0, 0, 1]],
}
TEMPLATES = {
    "old_ANF": lambda: old_ANF_template(max_and=3, max_xor=3),
    "fast": fast_template,
    "arman": arman_template,
}

# 20 keystreams of 9 to 16 bits. Products of classes from different fields
# only arise from the third block on: a downstream block's chaining multiplies
# bits of the block before it, and those bits mix the roots of the two blocks
# above. So every case is a chain (a, b, c, ...) whose first two block sizes
# share a factor: subfield pairs 2|4, 2|6, 3|6, 3|9, 4|8, 2|8, and pairs that
# share a factor without either dividing the other, (4,6), (6,4), (8,6). A
# source block much smaller than the next is avoided: old_ANF_template's
# rejection sampling can fail to finish there.
CASES = [
    ((4, 2, 3), "fast", 0), ((6, 2, 3), "fast", 0), ((6, 3, 2), "arman", 1),
    ((6, 3, 2), "old_ANF", 0), ((6, 3, 4), "fast", 0), ((6, 3, 5), "arman", 2),
    ((6, 4, 3), "arman", 0), ((6, 4, 3), "old_ANF", 2), ((6, 4, 2), "fast", 0),
    ((6, 4, 5), "arman", 1), ((4, 6, 2), "arman", 2), ((4, 6, 3), "fast", 1),
    ((8, 4, 3), "arman", 0), ((8, 2, 3), "fast", 0), ((8, 6, 2), "arman", 1),
    ((3, 6, 2), "arman", 0), ((9, 3, 2), "arman", 0), ((6, 4, 3, 2), "fast", 1),
    ((4, 2, 3, 5), "arman", 0), ((8, 4, 3), "fast", 1),
]


@pytest.mark.slow
@pytest.mark.parametrize(("sizes", "template", "seed"), CASES)
def test_non_coprime_cmpr_upper_bound_holds_against_berlekamp_massey(sizes, template, seed):
    # Blocks get independently drawn primitive polynomials and update
    # polynomials, so their generators are not compatible across the subfield
    # embeddings. Both RNGs are seeded: the templates draw through np.random.
    random.seed(seed)
    np.random.seed(seed)
    blocks = []
    for size in sizes:
        update = [0] * size
        while not any(update):
            update = [random.randint(0, 1) for _ in range(size)]
        blocks.append(MPR(size, random.choice(PRIMITIVE[size]), update))
    C = CMPR(blocks)
    C.generate_chaining(template=TEMPLATES[template]())

    # a bit of the most downstream block, which receives chaining
    bit = C.blocks[-1][0]
    lower, upper = C.estimate_LC(bit)
    # the case exercises the cap: its root expression multiplies classes from
    # fields whose degrees share a factor
    fields = C.root_expressions()[bit].root_table
    assert any(math.gcd(a, b) > 1 for key in fields for a, b in itertools.combinations(key, 2))

    # BM needs 2 * LC samples, and LC <= upper if the bound holds
    F = FeedbackRegister(2**C.size - 1, C)
    seq = [state[bit] for state in F.run(limit=2 * upper + 20, compiled=False)]
    lc, _ = berlekamp_massey(seq)
    assert lc <= upper
    # lower is a statistical estimate, not a bound (it often exceeds lc for
    # non-coprime blocks); it must not exceed upper
    assert lower <= upper
