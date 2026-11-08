"""Tests for BooleanFunction.merge_redundant.

merge_redundant flattens nested associative gates (XOR, AND, OR) into one gate,
drops one-argument gates, and leaves subfunctions -- nodes with more than one
parent -- whole. The order of a merged gate's arguments is unspecified, so these
tests compare arguments as sets and functions by truth table.
"""
import itertools
import random

import pytest

from PyPR.BooleanLogic import AND, NAND, NOR, NOT, OR, VAR, XNOR, XOR, BooleanFunction


@pytest.mark.parametrize("in_place", [False, True])
def test_merging_preserves_the_function_on_random_dags(in_place):
    """Whatever is merged, the truth table must not change. The DAGs mix every
    case the merge distinguishes: nested gates of the same kind, one-argument
    gates, repeated arguments, gates that do not merge, and shared nodes."""
    random.seed(7)
    num_vars = 5
    inputs = list(itertools.product((0, 1), repeat=num_vars))

    for _ in range(150):
        nodes: list[BooleanFunction] = [VAR(i) for i in range(num_vars)]
        for _ in range(random.randint(3, 25)):
            gate = random.choice([XOR, XOR, AND, OR, NAND, NOT])
            if gate is NOT:
                nodes.append(NOT(random.choice(nodes)))
            else:
                # drawing with replacement from every earlier node gives repeated
                # arguments and nodes with several parents
                arity = random.choice([1, 2, 2, 3, 5])
                nodes.append(gate(*random.choices(nodes, k=arity)))
        fn = nodes[-1]
        expected = [fn.eval(x) for x in inputs]

        if in_place:
            merged = fn.copy().merge_redundant(in_place=True)
        else:
            merged = fn.merge_redundant()
            # the original is left as it was
            assert [fn.eval(x) for x in inputs] == expected
        assert [merged.eval(x) for x in inputs] == expected


def test_a_nested_chain_becomes_one_gate():
    a, b, c, d = VAR(0), VAR(1), VAR(2), VAR(3)
    for gate in (XOR, AND, OR):
        merged = gate(gate(gate(a, b), c), d).merge_redundant()
        assert type(merged) is gate
        assert set(merged.args) == {a, b, c, d}
        assert len(merged.args) == 4


def test_repeated_terms_cancel_in_xor_and_collapse_in_and_or():
    a, b, c = VAR(0), VAR(1), VAR(2)
    # x ^ x = 0: a appears twice across the nesting and cancels
    assert set(XOR(XOR(a, b), XOR(a, c)).merge_redundant().args) == {b, c}
    # x & x = x and x | x = x: a is kept once
    for gate in (AND, OR):
        merged = gate(gate(a, b), gate(a, c)).merge_redundant()
        assert set(merged.args) == {a, b, c}
        assert len(merged.args) == 3


def test_an_xor_that_cancels_completely_is_the_constant_zero():
    """The empty sum is 0. It used to come back as an XOR with no arguments,
    which cannot be evaluated."""
    a, b = VAR(0), VAR(1)
    merged = AND(XOR(XOR(a, b), XOR(b, a)), b).merge_redundant()
    for x in itertools.product((0, 1), repeat=2):
        assert merged.eval(x) == 0


def test_a_single_term_left_over_needs_no_gate():
    a, b = VAR(0), VAR(1)
    assert XOR(XOR(a, b), a).merge_redundant() is b
    assert AND(a, AND(a, a)).merge_redundant() is a


@pytest.mark.parametrize("gate", [NAND, NOR, XNOR])
def test_a_one_argument_negated_gate_keeps_its_negation(gate):
    """A one-argument XOR, AND or OR is its argument and can be dropped. A
    one-argument NAND, NOR or XNOR is the argument's negation, and dropping it
    the same way inverted the function."""
    a, b = VAR(0), VAR(1)
    fn = gate(AND(a, b))
    merged = fn.merge_redundant()
    assert type(merged) is NOT
    for x in itertools.product((0, 1), repeat=2):
        assert merged.eval(x) == fn.eval(x)


def test_a_subfunction_is_kept_whole_and_stays_shared():
    """A node with two parents is not flattened into either: both still need
    it, and after merging both must refer to one object."""
    a, b, c, d = VAR(0), VAR(1), VAR(2), VAR(3)
    shared = XOR(a, b)
    merged = AND(XOR(shared, c), XOR(shared, d)).merge_redundant()

    left, right = merged.args
    kept = [arg for arg in left.args if type(arg) is XOR]
    assert len(kept) == 1
    assert set(kept[0].args) == {a, b}
    assert any(arg is kept[0] for arg in right.args)


def test_a_result_passed_through_a_one_argument_gate_is_not_consumed():
    """XOR(shared) passes shared's merged result straight through, so the gate
    above it sees a result that also belongs to shared's other parent. It must
    copy those terms, not take them over, or the other parent's view of shared
    would change with it."""
    a, b, c = VAR(0), VAR(1), VAR(2)
    shared = XOR(a, b)
    fn = AND(XOR(XOR(shared), c), shared)
    merged = fn.merge_redundant()

    absorbed, kept = merged.args
    assert set(absorbed.args) == {a, b, c}
    assert set(kept.args) == {a, b}
    for x in itertools.product((0, 1), repeat=3):
        assert merged.eval(x) == fn.eval(x)


def test_a_deep_chain_merges_without_recursion():
    """Loop unrolling builds chains thousands deep. Each level folds one term
    into the terms taken over from the level below, so the chain is merged in a
    single pass with no recursion."""
    depth = 20_000
    variables = [VAR(i) for i in range(depth)]
    fn = variables[0]
    for var in variables[1:]:
        fn = XOR(fn, var)

    merged = fn.merge_redundant()
    assert len(merged.args) == depth
    assert set(merged.args) == set(variables)
