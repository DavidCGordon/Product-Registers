"""Tests for PyPR.BooleanLogic.Optimization and its use by FeedbackFunction.compile.

`optimize` rewrites several functions as one DAG with fewer two-input
operations. Two things must hold whatever it does: every output computes the
same function of the bits, and the operation count never goes up. The rest of
the tests pin down each rewrite on an example small enough to count by hand.
"""
import itertools
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import (
    AND,
    CONST,
    NAND,
    NOR,
    NOT,
    OR,
    VAR,
    XNOR,
    XOR,
    BooleanANF,
    BooleanFunction,
    Optimization,
)
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template
from PyPR.BooleanLogic.Optimization import operation_count, optimize
from PyPR.BooleanLogic.SAT import TseytinFuse

from PyPR.FeedbackFunctions import (
    CMPR,
    FCSR,
    MPR,
    CrossJoin,
    FeedbackFunction,
    Fibonacci,
    TFunction,
)
from PyPR.FeedbackRegister import FeedbackRegister


@pytest.mark.parametrize("pair_budget", [Optimization.PAIR_BUDGET, 40, 0])
def test_random_functions_keep_their_truth_tables_and_never_cost_more(pair_budget, monkeypatch):
    """The two promises, on DAGs mixing every gate, constants, repeated
    arguments, one-argument gates and nodes shared between outputs. A small
    pair budget makes extraction stop partway or skip large gates, which must
    be just as correct."""
    monkeypatch.setattr(Optimization, "PAIR_BUDGET", pair_budget)
    random.seed(11)
    num_vars = 6
    inputs = list(itertools.product((0, 1), repeat=num_vars))

    for _ in range(300):
        nodes: list[BooleanFunction] = [VAR(i) for i in range(num_vars)]
        nodes += [CONST(0), CONST(1)]
        for _ in range(random.randint(3, 40)):
            gate = random.choice([XOR, XOR, AND, AND, OR, NAND, NOR, XNOR, NOT])
            if gate is NOT:
                nodes.append(NOT(random.choice(nodes)))
            else:
                nodes.append(gate(*random.choices(nodes, k=random.choice([1, 2, 2, 3, 5]))))
        functions = random.sample(nodes[num_vars + 2:], k=min(4, len(nodes) - num_vars - 2))
        before = [function.dense_str() for function in functions]

        optimized = optimize(functions)

        assert len(optimized) == len(functions)
        for function, result in zip(functions, optimized, strict=True):
            for x in inputs:
                assert int(result.eval(x)) == int(function.eval(x))
        assert operation_count(optimized) <= operation_count(functions)
        # the inputs are not edited
        assert [function.dense_str() for function in functions] == before


def test_operation_count_counts_each_shared_node_once():
    a, b, c, d = (VAR(i) for i in range(4))
    shared = AND(a, b)
    # the AND once, an XOR in each output; a NOT is one, a NAND of three is 2 + 1
    assert operation_count([XOR(shared, c), XOR(shared, d)]) == 3
    assert operation_count([NOT(a), NAND(a, b, c)]) == 4
    # equal logic in separate objects is counted separately: this measures the
    # DAG as built, which is what `optimize` then improves on
    assert operation_count([XOR(AND(a, b), c), XOR(AND(a, b), d)]) == 4


def test_equal_logic_in_different_outputs_becomes_one_shared_node():
    """Structural hashing: the two ANDs are separate objects over separate VAR
    objects, with their arguments in different orders, and still merge."""
    f0 = XOR(AND(VAR(0), VAR(1)), VAR(2))
    f1 = XOR(AND(VAR(1), VAR(0)), VAR(3))
    g0, g1 = optimize([f0, f1])

    and_0 = [arg for arg in g0.args if type(arg) is AND]
    and_1 = [arg for arg in g1.args if type(arg) is AND]
    assert len(and_0) == 1
    assert len(and_1) == 1
    assert and_0[0] is and_1[0]
    assert operation_count([g0, g1]) == 3


def test_a_pair_of_terms_shared_by_several_gates_is_computed_once():
    """Pair extraction (Paar's greedy algorithm). Three sums each contain
    a ^ b: as written that is 2 XORs each, 6 in all. With t = a ^ b computed
    once, each sum is t ^ x: 1 + 3 = 4."""
    a, b, c, d, e = (VAR(i) for i in range(5))
    functions = [XOR(a, b, c), XOR(a, b, d), XOR(a, b, e)]
    assert operation_count(functions) == 6
    assert operation_count(optimize(functions)) == 4

    # the same for a common two-literal cube among products
    functions = [AND(a, b, c), AND(a, b, d), AND(a, b, e)]
    assert operation_count(optimize(functions)) == 4


def test_a_common_literal_is_factored_out_of_a_sum_of_products():
    """x&a ^ x&b ^ x&c is 3 ANDs and 2 XORs. Factored, x & (a ^ b ^ c) is 1 AND
    and 2 XORs. AND distributes over OR the same way, and OR over AND."""
    x, a, b, c = (VAR(i) for i in range(4))
    assert operation_count(optimize([XOR(AND(x, a), AND(x, b), AND(x, c))])) == 3
    assert operation_count(optimize([OR(AND(x, a), AND(x, b), AND(x, c))])) == 3
    assert operation_count(optimize([AND(OR(x, a), OR(x, b), OR(x, c))])) == 3


def test_a_product_with_other_users_is_not_factored():
    """Factoring x out of x&a ^ x&b would stop computing x&a, which the second
    output needs anyway; the product is shared instead."""
    x, a, b = VAR(0), VAR(1), VAR(2)
    product = AND(x, a)
    functions = [XOR(product, AND(x, b)), product]
    optimized = optimize(functions)
    assert operation_count(optimized) == 3
    assert any(arg is optimized[1] for arg in optimized[0].args)


@pytest.mark.parametrize(("function", "expected"), [
    (XOR(VAR(0), VAR(0)), 0),                       # x ^ x = 0
    (AND(VAR(0), NOT(VAR(0))), 0),                  # x & ~x = 0
    (OR(VAR(0), NOT(VAR(0))), 1),                   # x | ~x = 1
    (XNOR(VAR(0), VAR(0)), 1),
    (AND(VAR(0), CONST(0)), 0),
    (OR(VAR(1), CONST(1)), 1),
    (XOR(AND(VAR(0), VAR(1)), AND(VAR(1), VAR(0))), 0),
], ids=["xor-cancel", "and-contradiction", "or-tautology", "xnor-self", "and-zero", "or-one", "equal-cubes"])
def test_constant_functions_are_reduced_to_a_constant(function, expected):
    result, = optimize([function])
    assert type(result) is CONST
    assert result.value == expected


def test_negations_are_simplified():
    a, b = VAR(0), VAR(1)
    # a double negation is the variable itself, at no cost
    result, = optimize([NOT(NOT(a))])
    assert type(result) is VAR
    assert result.index == 0
    # x ^ 1 is a complement, and a complement inside an XOR moves to its output
    assert operation_count(optimize([XOR(a, CONST(1))])) == 1
    assert operation_count(optimize([XOR(NOT(a), NOT(b))])) == 1
    # De Morgan: ~a & ~b is ~(a | b), one complement where there were two
    assert operation_count([AND(NOT(a), NOT(b))]) == 3
    assert operation_count(optimize([AND(NOT(a), NOT(b))])) == 2


def test_nodes_it_does_not_understand_are_kept_over_optimized_arguments():
    """A constant that is not a bit and a gate of another class are opaque:
    they keep their class and meaning, and only their arguments are rewritten."""
    a, b, c = VAR(0), VAR(1), VAR(2)
    anf = BooleanANF([[0, 1]])
    kept_constant, = optimize([CONST(anf)])
    assert type(kept_constant) is CONST
    assert kept_constant.value == anf

    fused = TseytinFuse(AND(VAR(0), VAR(1)))(XOR(a, a, b), c)
    result, = optimize([fused])
    assert type(result) is type(fused)
    # its first argument, a ^ a ^ b, was simplified to b
    assert type(result.args[0]) is VAR
    assert result.args[0].index == 1
    for x in itertools.product((0, 1), repeat=3):
        assert int(result.eval(x)) == int(fused.eval(x))


def test_the_result_does_not_depend_on_the_run():
    """Every choice is broken by node number, never by hash order, so the same
    input gives the same program -- and so the same compiled code."""
    random.seed(4)
    pool = [tuple(random.sample(range(12), 3)) for _ in range(20)]
    functions = [XOR(*[AND(*[VAR(i) for i in cube]) for cube in random.sample(pool, 6)]) for _ in range(8)]
    first = [function.dense_str() for function in optimize(functions)]
    second = [function.dense_str() for function in optimize(functions)]
    assert first == second


def test_feedback_function_optimize_computes_the_same_update():
    random.seed(2)
    np.random.seed(2)  # chaining also draws from numpy's generator
    cmpr = CMPR([MPR(7, "65"), MPR(5, "12"), MPR(3, "5")])
    cmpr.generate_chaining(template=arman_template(max_and=3))

    optimized = cmpr.optimize()
    assert type(optimized) is FeedbackFunction
    assert optimized.size == cmpr.size
    assert operation_count(optimized.fn_list) <= operation_count(cmpr.fn_list)
    for _ in range(200):
        state = [random.randint(0, 1) for _ in range(cmpr.size)]
        assert [int(f.eval(state)) for f in optimized.fn_list] == [int(f.eval(state)) for f in cmpr.fn_list]


@pytest.mark.parametrize("name", ["cmpr", "fibonacci", "crossjoin", "fcsr", "tfunction", "unrolled", "dense"])
def test_compiled_registers_clock_the_same_with_and_without_optimization(name):
    """The compiled update must match plain evaluation of fn_list, optimized or
    not. The cases cover both paths through compile: updates the optimizer
    improves, which are emitted as a shared straight-line program, and ones it
    cannot (a plain LFSR), which are compiled as written."""
    random.seed(6)
    np.random.seed(6)  # chaining also draws from numpy's generator
    if name == "cmpr":
        fn = CMPR([MPR(7, "65"), MPR(5, "12"), MPR(3, "5")])
        fn.generate_chaining(template=arman_template(max_and=3))
    elif name == "fibonacci":
        fn = Fibonacci(7, "65")
    elif name == "crossjoin":
        fn = CrossJoin(7, "69")
        fn.generate_nonlinearity(max_ands=3, tap_density=0.75)
    elif name == "fcsr":
        fn = FCSR(8, 2**8 - 5)
    elif name == "tfunction":
        fn = TFunction(8)
    elif name == "unrolled":
        # four clocks of a chained CMPR as one update: heavy sharing between bits
        base = CMPR([MPR(5, "12"), MPR(3, "5")])
        base.generate_chaining(template=arman_template(max_and=2))
        steps = list(base.iterator(5))
        fn = FeedbackFunction(list(steps[-1]))
    else:
        # products drawn from a small pool, so outputs share cubes and literals
        pool = [tuple(sorted(random.sample(range(10), 2))) for _ in range(12)]
        fn = FeedbackFunction([
            XOR(VAR((bit + 1) % 10), *[AND(*[VAR(i) for i in cube]) for cube in random.sample(pool, 5)])
            for bit in range(10)
        ])

    seed = [random.randint(0, 1) for _ in range(fn.size)]
    expected = [str(state) for state in FeedbackRegister(seed, fn).run(64, compiled=False)]

    # compile passes optimize=True to generate_python
    fn.compile()
    assert [str(state) for state in FeedbackRegister(seed, fn).run(64)] == expected
    # the other compiled entry point, which returns a new array
    state = np.array(seed, dtype=np.uint8)
    assert list(fn._compiled(state)) == [int(f.eval(seed)) for f in fn.fn_list]

    fn.compile(optimize=False)
    assert [str(state) for state in FeedbackRegister(seed, fn).run(64)] == expected


def test_generate_python_optimizes_only_when_asked():
    """The flag is off by default, so the statements follow the function as
    written; with it on they follow the optimized form. Both must compute the
    function, which is checked by running them."""
    a, b, c = VAR(0), VAR(1), VAR(2)
    # a&b appears twice and cancels, leaving a&c
    function = XOR(AND(a, b), AND(a, c), AND(b, a))
    assert function.generate_python() == [
        "output = ((array[0] & array[1]) ^ (array[0] & array[2]) ^ (array[1] & array[0]))"
    ]
    assert function.generate_python(optimize=True) == ["output = (array[0] & array[2])"]

    for statements in (function.generate_python(), function.generate_python(optimize=True)):
        for array in itertools.product((0, 1), repeat=3):
            namespace = {"array": array}
            exec("\n".join(statements), namespace)
            assert namespace["output"] == int(function.eval(array))

    # an override names a node of the function as written, which the optimized
    # form does not contain
    with pytest.raises(ValueError, match="overrides cannot be combined"):
        function.generate_python(optimize=True, overrides={function.args[0]: "g"})


def test_generate_python_leaves_a_function_it_cannot_improve_as_written():
    function = XOR(AND(VAR(0), VAR(1)), VAR(2))
    assert function.generate_python(optimize=True) == function.generate_python()

    lfsr = Fibonacci(7, "65")
    assert lfsr.generate_python(optimize=True) == lfsr.generate_python()


def test_feedback_function_generate_python_computes_the_update():
    """One statement list for the whole register: with the flag off each bit is
    written from its own function, and with it on the bits share temporaries.
    Either way, running the statements must give the next state."""
    random.seed(6)
    pool = [tuple(sorted(random.sample(range(10), 2))) for _ in range(12)]
    fn = FeedbackFunction([
        XOR(VAR((bit + 1) % 10), *[AND(*[VAR(i) for i in cube]) for cube in random.sample(pool, 5)])
        for bit in range(10)
    ])
    plain = fn.generate_python()
    optimized = fn.generate_python(output_name="out", array_name="state", optimize=True)
    assert optimized != plain
    # products shared between bits are computed once, into temporaries
    assert any(line.startswith("t0 = ") for line in optimized)
    assert not any(line.startswith("t0 = ") for line in plain)

    for _ in range(50):
        state = [random.randint(0, 1) for _ in range(fn.size)]
        expected = [int(f.eval(state)) for f in fn.fn_list]
        namespace = {"curr_state": state, "next_state": [0] * fn.size}
        exec("\n".join(plain), namespace)
        assert namespace["next_state"] == expected
        namespace = {"state": state, "out": [0] * fn.size}
        exec("\n".join(optimized), namespace)
        assert namespace["out"] == expected
