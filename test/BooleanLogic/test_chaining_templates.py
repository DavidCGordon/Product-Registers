"""Tests for the chaining-template combinators and the templates built from them.

A template is a small program that samples a Boolean function. Combinators nest:
GATE applies a gate to what its sources sample, the filters (NONCONSTANT,
DISTINCT, UNIQUE) resample until a condition holds, and OPTIONAL sometimes
contributes nothing at all.

The one non-obvious rule is `allow_empty_return`. For OPTIONAL it would mean
"you may drop"; for the filters it means "give up instead of retrying". Because
both read the same flag, no combinator can grant it to OPTIONAL without changing
how the filters beneath it behave -- which is why OPTIONAL drops regardless of
the flag, and why these tests pin both sides.
"""
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, VAR
from PyPR.BooleanLogic.ChainingGeneration.TemplateBuilding import (
    GATE,
    NONCONSTANT,
    OPTIONAL,
    VALUE,
)
from PyPR.BooleanLogic.ChainingGeneration.Templates import (
    fast_template,
    three_majority_template,
)

from PyPR.FeedbackFunctions import CMPR, MPR

# ── OPTIONAL ─────────────────────────────────────────────────────────────────

def test_optional_drops_its_input_from_a_gate():
    """OPTIONAL's home: GATE(AND, [a, OPTIONAL(b)]) is sometimes AND(a)."""
    random.seed(0)
    template = GATE({"gate_class": AND},
                    sources=[VALUE(VAR(0)), OPTIONAL({"drop_chance": 1.0}, source=VALUE(VAR(1)))])
    arities = {len(template.sample()[0].args) for _ in range(10)}
    assert arities == {1}

def test_optional_that_never_drops_passes_its_source_through():
    random.seed(0)
    template = GATE({"gate_class": AND},
                    sources=[VALUE(VAR(0)), OPTIONAL({"drop_chance": 0.0}, source=VALUE(VAR(1)))])
    arities = {len(template.sample()[0].args) for _ in range(10)}
    assert arities == {2}

def test_optional_drops_whatever_the_flag_says():
    """The flag can't gate the drop: nothing that hosts OPTIONAL grants it."""
    random.seed(0)
    optional = OPTIONAL({"drop_chance": 1.0}, source=VALUE(VAR(0)))
    assert optional._sample(allow_empty_return=False) == []
    assert optional._sample(allow_empty_return=True) == []

def test_a_filter_rejects_an_optional_that_drops():
    """NONCONSTANT demands an output from its source; an empty one is a
    composition it doesn't support, reported rather than silently accepted."""
    random.seed(0)
    template = NONCONSTANT({"attempt_limit": 3},
                           source=OPTIONAL({"drop_chance": 1.0}, source=VALUE(VAR(0))))
    with pytest.raises(ValueError, match="expects an output from its source"):
        template.sample()


# ── three_majority_template ──────────────────────────────────────────────────
# Both requirements are consumed per chained block. An int means the same value
# in every block; anything else is the per-block sequence. Only the int form used
# to work -- a list left the per-block names unbound and raised NameError.

def _three_majority_chaining(correlation_immunity, algebraic_degree):
    """Chaining sampled from a seeded template, as comparable strings."""
    np.random.seed(5)
    random.seed(5)
    cmpr = CMPR([MPR(7, "12"), MPR(7, "12"), MPR(7, "12")])
    functions = three_majority_template(correlation_immunity, algebraic_degree)(cmpr)
    return {index: fn.pretty_str() for index, fn in functions.items()}

def test_three_majority_accepts_per_block_sequences():
    as_ints = _three_majority_chaining(1, 1)
    assert _three_majority_chaining([1, 1, 1], [1, 1, 1]) == as_ints
    assert _three_majority_chaining((1, 1, 1), (1, 1, 1)) == as_ints

def test_three_majority_accepts_a_mix_of_int_and_sequence():
    assert _three_majority_chaining(1, [1, 1, 1]) == _three_majority_chaining(1, 1)

def test_three_majority_rejects_requirements_its_blocks_cannot_meet():
    """(ci + 1)(2 ad - 1) variables are needed; a 7-bit block can't host 3 and 4."""
    with pytest.raises(ValueError, match="Unable to complete 3-Maj"):
        _three_majority_chaining(3, 4)

def test_three_majority_accepts_one_entry_per_chained_block():
    """The last block feeds nothing, so its entry may be left off."""
    assert _three_majority_chaining([1, 1], [1, 1]) == _three_majority_chaining(1, 1)

@pytest.mark.parametrize("short_or_long", [[1], [1, 1, 1, 1]], ids=["short", "long"])
def test_three_majority_rejects_a_sequence_of_the_wrong_length(short_or_long):
    """A short list used to fail partway through generation with an IndexError."""
    with pytest.raises(ValueError, match="entries, but this CMPR has 3 blocks"):
        _three_majority_chaining(short_or_long, 1)


# ── reproducibility ──────────────────────────────────────────────────────────

def test_chaining_is_reproducible_when_both_generators_are_seeded():
    """Templates sample from Python's `random` (drop chances, probabilistic
    templates) and from numpy's generator (`SAMPLE` picks sources with
    np.random.choice). Seeding only one leaves the chaining different in every
    process, so a test that generates chaining must seed both."""
    chainings = []
    for _ in range(2):
        random.seed(11)
        np.random.seed(11)
        C = CMPR([MPR(7, "65"), MPR(5, "12"), MPR(3, "5")])
        C.generateChaining(template=fast_template())
        chainings.append([fn.dense_str() for fn in C.fn_list])
    assert chainings[0] == chainings[1]
