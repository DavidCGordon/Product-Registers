"""Every gate must mean the same thing in every representation it has.

A gate is defined by `eval`, and three other paths each claim to compute the
same function: the Tseytin clauses `sat()` and `enum_models()` reason over,
the Python that `compile()` generates and numba compiles, and the C and VHDL
backends. Nothing forces them to agree -- each gate implements each path
separately -- so these tests compare every path against `eval`, gate by gate
and arity by arity.

They exist because they would have failed. Until this was written, the SAT
encoding of a one-input XNOR, NAND or NOR asserted the output *equalled* its
input, and a NAND or NOR of three or more inputs swapped its output wire with an
input wire in the final step -- so `sat()` and `functionally_equivalent()`
answered for the wrong function while `eval` and compiled code were right.
Separately, C generation raised for XNOR, NAND, NOR and NOT, and `sat()` raised
on any gate with no arguments.
"""
import itertools

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, CONST, NOT, OR, VAR, XOR
from PyPR.BooleanLogic.Gates import NAND, NOR, XNOR

_GATES = [XOR, AND, OR, XNOR, NAND, NOR]
_GATE_IDS = [gate.__name__ for gate in _GATES]


# ── the SAT encoding agrees with eval ────────────────────────────────────────

@pytest.mark.parametrize("arity", [1, 2, 3, 4, 5])
@pytest.mark.parametrize("gate", _GATES, ids=_GATE_IDS)
def test_sat_models_are_exactly_the_inputs_eval_accepts(gate, arity):
    """The satisfying assignments of the clauses are the gate's on-set.

    Checked for the gate and for its negation: the negation's models are the
    off-set, so the two together pin every input, not just the accepted ones.
    `enum_models` may omit variables that don't matter to a model, so each model
    is expanded over its missing variables before comparing.
    """
    inputs = [VAR(i) for i in range(arity)]
    for fn in (gate(*inputs), NOT(gate(*inputs))):
        expected = {
            bits for bits in itertools.product([0, 1], repeat=arity)
            if fn.eval(list(bits))
        }
        models = set()
        for model in fn.enum_models():
            missing = [i for i in range(arity) if i not in model]
            for fill in itertools.product([0, 1], repeat=len(missing)):
                full = {**{k: int(v) for k, v in model.items()}, **dict(zip(missing, fill))}
                models.add(tuple(full[i] for i in range(arity)))

        assert models == expected, (
            f"{gate.__name__} of {arity} input(s): clauses accept {sorted(models)}, "
            f"eval accepts {sorted(expected)}"
        )

def test_nested_gates_are_equivalent_to_their_anf():
    """Mixed nesting, where a wrong wire in one gate corrupts its parent's input."""
    fn = NAND(XOR(VAR(0), VAR(1), VAR(2)), NOR(VAR(1), VAR(3), VAR(4)), XNOR(VAR(0), VAR(4)))
    for bits in itertools.product([0, 1], repeat=5):
        pinned = AND(fn if fn.eval(list(bits)) else NOT(fn),
                     *[VAR(i) if b else NOT(VAR(i)) for i, b in enumerate(bits)])
        assert pinned.sat() is not None, f"clauses reject an assignment eval accepts: {bits}"


# ── generated code agrees with eval ──────────────────────────────────────────

@pytest.mark.parametrize("arity", [1, 2, 3, 4])
@pytest.mark.parametrize("gate", _GATES, ids=_GATE_IDS)
def test_compiled_gate_agrees_with_eval(gate, arity):
    """`compile()` runs the generated Python through numba, so this covers
    `_generate_python` as well as the compiled path."""
    fn = gate(*[VAR(i) for i in range(arity)])
    compiled = fn.compile()
    for bits in itertools.product([0, 1], repeat=arity):
        # uint8 arrays, as registers pass their state -- a Python list goes
        # through numba's deprecated reflected-list path
        assert int(compiled(np.array(bits, dtype=np.uint8))) == int(fn.eval(list(bits))), (
            f"{gate.__name__} of {arity}: compiled differs from eval at {bits}"
        )

def test_compiled_not_agrees_with_eval():
    compiled = NOT(VAR(0)).compile()
    assert [int(compiled(np.array([b], dtype=np.uint8))) for b in (0, 1)] == [1, 0]

@pytest.mark.parametrize(("gate", "expected"), [
    (XOR,  "output = (array[0] ^ array[1]);"),
    (AND,  "output = (array[0] & array[1]);"),
    (OR,   "output = (array[0] | array[1]);"),
    (XNOR, "output = (!(array[0] ^ array[1]));"),
    (NAND, "output = (!(array[0] & array[1]));"),
    (NOR,  "output = (!(array[0] | array[1]));"),
], ids=_GATE_IDS)
def test_c_generation(gate, expected):
    """The negated gates generated C through the wrong hook and raised."""
    assert gate(VAR(0), VAR(1)).generate_c() == [expected]

def test_c_generation_for_not():
    assert NOT(VAR(0)).generate_c() == ["output = (!(array[0]));"]

@pytest.mark.parametrize("gate", _GATES + [NOT], ids=_GATE_IDS + ["NOT"])
def test_every_backend_generates_for_every_gate(gate):
    fn = NOT(VAR(0)) if gate is NOT else gate(VAR(0), VAR(1), VAR(2))
    for backend in (fn.generate_c, fn.generate_VHDL, fn.generate_python):
        assert backend(), f"{gate.__name__}.{backend.__name__} produced nothing"


# ── degenerate arities ───────────────────────────────────────────────────────

@pytest.mark.parametrize("gate", _GATES, ids=_GATE_IDS)
def test_a_gate_with_no_arguments_is_unsatisfiable(gate):
    """An empty gate is encoded as a contradiction -- a unit clause and its
    negation. They were built as `(x)` rather than `(x,)`, which is an int, not a
    1-tuple, so the solver was handed bare integers and raised."""
    assert gate().sat() is None
    assert AND(VAR(0), gate()).sat() is None

def test_not_rejects_more_than_one_argument_wire():
    """NOT has arg_limit = 1, so this is unreachable through normal construction;
    the hook must still fail loudly rather than fall through and return None."""
    node = NOT(VAR(0))
    node.args = (VAR(0), VAR(1))          # bypass arg_limit to reach the branch
    label_map = {node: [5], node.args[0]: [2], node.args[1]: [3]}
    with pytest.raises(ValueError, match="single argument"):
        node._tseytin_clauses(label_map)


# ── labels are threaded correctly across circuits ────────────────────────────

def test_threading_after_a_lone_false_constant_starts_at_a_real_variable():
    """CONST(0) labels itself -1, which used to make the next wire 0 -- and 0 is
    not a literal at all; it terminates a clause in DIMACS."""
    variable = VAR(0)
    data = CONST(0).tseytin()
    _clauses, node_labels, _variable_labels = variable.tseytin(*data)
    assert node_labels[variable] == [2]

def test_threading_reuses_the_wires_of_a_shared_subexpression():
    """The reason tseytin() accepts its own previous output: a node that appears
    in several circuits is encoded once, on one set of wires."""
    shared = AND(VAR(0), VAR(1))
    data = XOR(shared, VAR(2)).tseytin()
    wire = data[1][shared]
    clauses, node_labels, variable_labels = OR(shared, VAR(3)).tseytin(*data)

    assert node_labels[shared] == wire
    assert variable_labels[0] == data[2][0]
    separately = len(XOR(shared, VAR(2)).tseytin()[0]) + len(OR(shared, VAR(3)).tseytin()[0])
    assert len(clauses) < separately
