"""Tests for fusion nodes: subtrees that own their own SAT encoding.

A fusion node must denote exactly the subtree it wraps, so every test here is
ultimately a differential one -- the fused form against the expanded form, by
truth table, by SAT, or by both. The fused CNF is written by hand, and nothing
but a comparison like this establishes that it agrees.
"""
import itertools
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, CONST, NOT, OR, VAR, XOR
from PyPR.BooleanLogic.SAT import TseytinFuse, _TseytinFuse

# Templates are written over VAR(0)..VAR(n-1); the node supplies those inputs
# positionally, so the node denotes the template composed with its arguments.
MAJ3 = OR(AND(VAR(0), VAR(1)), AND(VAR(1), VAR(2)), AND(VAR(0), VAR(2)))


# ── the fused encoding denotes the wrapped subtree ───────────────────────────

def test_fused_node_matches_expansion_on_every_input():
    """Exhaustive truth table: a fusion node equals the subtree it wraps."""
    fused = TseytinFuse(MAJ3)(VAR(0), VAR(1), VAR(2))
    expanded = fused.expand()
    for bits in itertools.product([0, 1], repeat=3):
        assert fused.eval(list(bits)) == expanded.eval(list(bits)), (
            f"majority disagreed on {bits}"
        )

def test_verify_confirms_the_encoding_by_sat():
    """verify() is the solver-level check that the CNF matches the subtree."""
    assert TseytinFuse(MAJ3)(VAR(0), VAR(1), VAR(2)).verify()

def test_fused_node_is_equivalent_when_embedded():
    """A fusion node inside a larger function encodes the same function."""
    maj = TseytinFuse(MAJ3)
    fused = XOR(maj(VAR(0), VAR(1), VAR(2)), NOT(VAR(3)))
    unfused = XOR(MAJ3, NOT(VAR(3)))
    assert fused.functionally_equivalent(unfused)

def test_fusion_nodes_nest():
    """An outer fusion node treats an inner one as an ordinary argument."""
    maj = TseytinFuse(MAJ3)
    conj = TseytinFuse(AND(VAR(0), VAR(1), VAR(2), VAR(3)))
    nested = maj(conj(VAR(0), VAR(1), VAR(2), VAR(3)), VAR(4), VAR(5))
    assert nested.verify()


# ── degenerate templates ─────────────────────────────────────────────────────
# A template that allocates no interior wires still has to report a result wire,
# which is the case the label bookkeeping is easiest to get wrong on.

def test_leaf_and_constant_templates_encode():
    for template, arity in [(VAR(0), 1), (CONST(1), 1), (NOT(VAR(0)), 1)]:
        node = TseytinFuse(template)(*[VAR(i) for i in range(arity)])
        assert node.verify(), f"{template} template failed to verify"


# ── the placeholder direct encodings ─────────────────────────────────────────
# An n-input AND or OR is definable in one wire and n+1 clauses; chaining
# two-input gates spends n-1 wires and 3(n-1). Both must encode the same thing.

def test_direct_and_or_use_one_wire():
    """A direct form never labels its template: the whole domain is one wire."""
    for gate in (AND, OR):
        for n in (2, 4, 8):
            template = gate(*[VAR(i) for i in range(n)])
            node = TseytinFuse(template)(*[VAR(i) for i in range(n)])
            labels, _ = node.tseytin_labels()
            clauses = node.tseytin_clauses(labels)
            # every solver variable, less the reserved constant wire (+-1)
            wires = {w for ws in labels.values() for w in ws if abs(w) != 1}
            assert len(labels[node]) == 1, (
                f"{gate.__name__}{n} claimed {len(labels[node])} wires, so it did not fuse"
            )
            # n argument wires plus the one carrying the gate's result
            assert len(wires) == n + 1, (
                f"{gate.__name__}{n} used {len(wires)} wires, expected {n + 1}"
            )
            assert len(clauses) == n + 1, (
                f"{gate.__name__}{n} produced {len(clauses)} clauses, expected {n + 1}"
            )
            assert node.verify()

def test_direct_encoding_beats_the_chain():
    """The fused form is strictly cheaper than expanding, past two inputs."""
    template = AND(*[VAR(i) for i in range(8)])
    fused = TseytinFuse(template)(*[VAR(i) for i in range(8)])
    plain = AND(*[VAR(i) for i in range(8)])

    f_labels, _ = fused.tseytin_labels()
    p_labels, _ = plain.tseytin_labels()
    # solver variables in each, less the reserved constant wire (+-1)
    fused_wires = {w for ws in f_labels.values() for w in ws if abs(w) != 1}
    plain_wires = {w for ws in p_labels.values() for w in ws if abs(w) != 1}
    assert len(fused_wires) < len(plain_wires)
    assert len(fused.tseytin_clauses(f_labels)) < len(plain.tseytin_clauses(p_labels))

def test_template_without_a_direct_form_falls_back():
    """XOR has no compact direct CNF, so it expands and claims interior wires."""
    node = TseytinFuse(XOR(VAR(0), VAR(1), VAR(2)))(VAR(0), VAR(1), VAR(2))
    labels, _ = node.tseytin_labels()
    assert len(labels[node]) > 1, "expected the fallback expansion, not a direct form"
    assert node.template not in labels, "the template's interior should not be keyed by node"
    assert node.verify()

def test_template_is_copied_on_construction():
    """Mutating the template afterwards must not reach a node built from it."""
    template = AND(VAR(0), VAR(1))
    node = TseytinFuse(template)(VAR(5), VAR(6))
    assert node.template is not template
    template.add_arguments(VAR(2))
    assert len(node.template.args) == 2, "the node's template followed the caller's mutation"
    assert node.verify()

def test_two_nodes_from_one_factory_do_not_share_wires():
    """Copies keep the label map, which is keyed by identity, unambiguous."""
    factory = TseytinFuse(MAJ3)
    a = factory(VAR(0), VAR(1), VAR(2))
    b = factory(VAR(3), VAR(4), VAR(5))
    assert a.template is not b.template
    labels, _ = XOR(a, b).tseytin_labels()
    assert set(labels[a]).isdisjoint(labels[b])
    assert XOR(a, b).functionally_equivalent(XOR(a.expand(), b.expand()))


# ── arguments are fixed by the template ──────────────────────────────────────

def test_argument_mutation_is_rejected():
    """Mutating args would desync the boundary from the template it encodes."""
    node = TseytinFuse(MAJ3)(VAR(0), VAR(1), VAR(2))
    for mutate in (node.add_arguments, node.remove_arguments):
        with pytest.raises(ValueError, match="fixed by the template"):
            mutate(VAR(3))


# ── differential over random surroundings ────────────────────────────────────

def test_fused_and_unfused_agree_on_random_functions():
    """Fusing a subtree changes the encoding, never the function."""
    random.seed(42)
    maj = TseytinFuse(MAJ3)
    for trial in range(25):
        idxs = [random.randrange(4) for _ in range(4)]
        outer = random.choice([AND, OR, XOR])
        fused = outer(maj(*[VAR(i) for i in idxs[:3]]), VAR(idxs[3]))
        unfused = outer(
            MAJ3.compose({j: VAR(idxs[j]) for j in range(3)}), VAR(idxs[3])
        )
        assert fused.functionally_equivalent(unfused), f"trial {trial}: {idxs}, {outer.__name__}"
        for bits in itertools.product([0, 1], repeat=4):
            assert fused.eval(list(bits)) == unfused.eval(list(bits)), (
                f"trial {trial}: disagreed on {bits}"
            )


# ── verify() has teeth ───────────────────────────────────────────────────────
# Every test above asserts a correct encoding verifies. That alone would pass
# against a verify() that always returned True, so these deliberately encode
# something other than the template and require it to be caught. They are the
# tests that exercise the case a hand-written fusion actually fails in.

class _EncodesTheWrongGate(_TseytinFuse):
    """Encodes a disjunction whatever the template says."""

    def _direct_skeleton(self, result):
        args = [self._arg_wires[i] for i in range(len(self.args))]
        output = max(args + [result]) + 1
        return [output], (
            [(-arg, output) for arg in args] + [tuple([-output] + args)]
        )

class _DropsOneImplication(_TseytinFuse):
    """A conjunction encoding missing `output -> arg_labels[0]`.

    Under-constrains the output wire rather than mis-defining it, which is the
    shape a hand-written CNF is most likely to get wrong: every clause present
    is correct, and the encoding is merely too weak.
    """

    def _direct_skeleton(self, result):
        args = [self._arg_wires[i] for i in range(len(self.args))]
        output = max(args + [result]) + 1
        return [output], (
            [(-output, arg) for arg in args[1:]]
            + [tuple([output] + [-arg for arg in args])]
        )

def test_verify_rejects_an_encoding_of_the_wrong_gate():
    template = AND(VAR(0), VAR(1), VAR(2))
    node = _EncodesTheWrongGate(template, VAR(0), VAR(1), VAR(2))
    assert not node.verify(), "verify() accepted a disjunction for a conjunction"

def test_verify_rejects_a_single_missing_clause():
    template = AND(VAR(0), VAR(1), VAR(2))
    node = _DropsOneImplication(template, VAR(0), VAR(1), VAR(2))
    assert not node.verify(), "verify() accepted an under-constrained output wire"

def test_a_wrong_fusion_is_caught_when_embedded_too():
    """The mismatch has to survive being wired into a larger circuit."""
    template = AND(VAR(0), VAR(1), VAR(2))
    broken = _DropsOneImplication(template, VAR(0), VAR(1), VAR(2))
    sound = TseytinFuse(template)(VAR(0), VAR(1), VAR(2))
    assert not XOR(broken, VAR(3)).functionally_equivalent(XOR(sound, VAR(3)))


# ── a fusion node inside ordinary DAG operations ─────────────────────────────
# Fusion changes only the SAT encoding. Everything else the library does to a
# DAG -- copy it, substitute into it, translate it to ANF, generate code for it,
# compile it -- must treat a fusion node as the function it denotes. The base
# walks rebuild nodes through `_copy`, whose default would have passed the first
# argument where the template goes; codegen had no hooks at all, so a fused node
# could not be compiled or used in a register.

def _fused_circuit():
    """XOR(maj(x3, x4, x5), NOT x6), fused and unfused, over seven inputs."""
    node = TseytinFuse(MAJ3)(VAR(3), VAR(4), VAR(5))
    fused = XOR(node, NOT(VAR(6)))
    reference = XOR(MAJ3.compose({0: VAR(3), 1: VAR(4), 2: VAR(5)}), NOT(VAR(6)))
    states = [list(bits) for bits in itertools.product([0, 1], repeat=7)]
    return node, fused, reference, states

def test_copying_a_dag_keeps_the_fusion_node_and_its_meaning():
    node, fused, reference, states = _fused_circuit()
    copied = fused.__copy__()

    assert isinstance(copied.args[0], _TseytinFuse)
    assert copied.args[0] is not node, "a copy must not share the original node"
    # the template is the node's own data, not a child: a deep copy copies it too
    assert copied.args[0].template is not node.template, "the copy shares the template"
    assert all(copied.eval(s) == reference.eval(s) for s in states)

def test_composing_into_a_dag_with_a_fusion_node():
    _node, fused, _reference, states = _fused_circuit()
    swapped = fused.compose({6: VAR(0)})
    expected = XOR(MAJ3.compose({0: VAR(3), 1: VAR(4), 2: VAR(5)}), NOT(VAR(0)))
    assert all(swapped.eval(s) == expected.eval(s) for s in states)

def test_translating_a_dag_with_a_fusion_node_to_anf():
    _node, fused, reference, states = _fused_circuit()
    anf = fused.translate_ANF()
    assert all(anf.eval(s) == reference.eval(s) for s in states)

def test_generated_code_is_the_templates_gates_over_the_arguments():
    node, _fused, _reference, _states = _fused_circuit()
    assert node.generate_c() == [
        "output = ((array[3] & array[4]) | (array[4] & array[5]) | (array[3] & array[5]));"
    ]
    assert node.generate_VHDL()
    assert node.generate_python()

def test_a_dag_with_a_fusion_node_compiles():
    """compile() builds on generate_python, so a fused node used to break it."""
    _node, fused, reference, states = _fused_circuit()
    compiled = fused.compile()
    for s in states:
        assert int(compiled(np.array(s, dtype=np.uint8))) == int(reference.eval(s))

def test_a_template_reading_an_input_it_is_not_given_is_rejected():
    """VAR(5) with three arguments used to become a free interior wire -- a
    different function from the one `eval` computes, which raised instead."""
    with pytest.raises(ValueError, match=r"reads inputs \[5\]"):
        TseytinFuse(XOR(VAR(0), VAR(5)))(VAR(0), VAR(1), VAR(2))
