"""The DAG walk every whole-function method is built on.

`postorder` yields each node once, after all of its arguments, and skips the
nodes a `stop` predicate claims. The methods built on it (eval, copy, compose,
code generation, the SAT walks, ids, statistics) differ only in what they
compute per node, so the contract below is what they all rely on.

The later tests are regressions from before the walk was shared: methods that
recursed without memoisation, and one that visited a shared leaf once per
parent.
"""
import itertools

from PyPR.BooleanLogic import AND, CONST, NOT, OR, VAR, XOR
from PyPR.BooleanLogic.Gates import NAND, NOR, XNOR


def test_arguments_come_before_their_gate_left_to_right():
    a, b, c = VAR(0), VAR(1), VAR(2)
    product = AND(a, b)
    fn = XOR(product, c)
    assert list(fn.postorder()) == [a, b, product, c, fn]

def test_a_shared_node_is_yielded_once():
    a, b = VAR(0), VAR(1)
    shared = AND(a, b)
    fn = XOR(shared, OR(shared, a), shared)
    nodes = list(fn.postorder())
    assert len(nodes) == len(set(nodes)) == 5
    # yielded where the depth-first search first finishes it
    assert nodes.index(shared) < nodes.index(fn.args[1])

def test_equal_but_distinct_leaves_are_distinct_nodes():
    """Nodes are identified by object, not by the variable they name."""
    assert len(list(AND(VAR(0), VAR(0)).postorder())) == 3

def test_a_leaf_yields_itself():
    leaf = CONST(1)
    assert list(leaf.postorder()) == [leaf]

def test_stop_cuts_off_a_node_and_everything_beneath_it():
    a, b, c = VAR(0), VAR(1), VAR(2)
    below = AND(a, b)
    fn = XOR(below, c)
    assert list(fn.postorder(stop=lambda n: n is below)) == [c, fn]

def test_stop_is_how_a_later_call_resumes_an_earlier_one():
    """generate_ids and tseytin_labels continue from their own previous output:
    nodes that already have an id keep it, and only new nodes are numbered."""
    shared = AND(VAR(0), VAR(1))
    ids = XOR(shared, VAR(2)).generate_ids()
    before = dict(ids)
    ids = OR(shared, VAR(3)).generate_ids(ids)
    assert all(ids[node] == index for node, index in before.items())
    assert len(ids) == len(before) + 2   # VAR(3) and the OR; the AND is reused

def test_deep_functions_do_not_hit_the_recursion_limit():
    """max_idx and idxs_used used to recurse on every argument."""
    fn = VAR(0)
    for i in range(5000):
        fn = XOR(fn, VAR(i % 7))
    assert fn.num_nodes() == 5000 + 1 + 5000
    assert fn.max_idx() == 6
    assert fn.idxs_used() == set(range(7))
    # 5001 leaves XORed together, each 1: odd parity
    assert fn.eval([1] * 7) == 1
    assert fn.copy().num_nodes() == fn.num_nodes()
    assert str(fn).startswith("XOR(XOR(")

def test_heavy_sharing_is_linear_not_exponential():
    """Each level uses the previous one twice: 2^60 paths through 61 nodes. A
    walk without memoisation (as idxs_used and max_idx were) never finishes."""
    fn = VAR(3)
    for _ in range(60):
        fn = XOR(fn, fn)
    assert fn.num_nodes() == 61
    assert fn.idxs_used() == {3}
    assert fn.max_idx() == 3
    # x xor x = 0, so every level above the first is 0
    assert fn.eval([0, 0, 0, 1]) == 0

def test_inputs_lists_each_leaf_object_once():
    v = VAR(0)
    fn = XOR(v, AND(v, VAR(1)))
    assert [id(leaf) for leaf in fn.inputs()] == [id(v), id(fn.args[1].args[1])]

def test_shift_and_remap_move_a_shared_leaf_once():
    """A VAR used by two gates was listed by inputs() once per gate, so shifting
    moved it twice: x0 + x0 x1 shifted by one came out as x2 + x2 x2."""
    v = VAR(0)
    fn = XOR(v, AND(v, VAR(1)))
    assert str(fn.shift_indices(1)) == "XOR(VAR(1),AND(VAR(1),VAR(2)))"
    assert str(fn.remap_indices({0: 1, 1: 2})) == "XOR(VAR(1),AND(VAR(1),VAR(2)))"
    # the original is untouched
    assert str(fn) == "XOR(VAR(0),AND(VAR(0),VAR(1)))"

def test_binarize_a_single_argument_negated_gate():
    """XNOR, NAND and NOR of one argument are its negation. Binarizing them
    reduced over all but the last argument -- none -- and raised."""
    for gate in (XNOR, NAND, NOR):
        fn = XOR(gate(VAR(0)), VAR(1))
        binary = fn.binarize()
        assert all(len(node.args) <= 2 for node in binary.postorder())
        for bits in itertools.product([0, 1], repeat=2):
            assert int(binary.eval(list(bits))) == int(fn.eval(list(bits))), gate.__name__
        assert isinstance(binary.args[0], NOT)
