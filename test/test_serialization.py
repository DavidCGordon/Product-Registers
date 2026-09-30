"""Tests for JSON serialization.

Covers to_JSON/from_JSON and to_file/from_file for every feedback function type
and for registers, each verified by comparing the output bit sequences of the
original and reconstructed registers; the shared-structure guarantee of
generate_JSON; compatibility with files written before the `Serializable` mixin;
and the value types BooleanANF and BooleanGF.
"""
import json
import os
import random
import tempfile
from pathlib import Path

import numpy as np
import pytest

from PyPR.JSON_Serialization import Serializable, generate_JSON, parse_JSON

from PyPR.BooleanLogic import AND, CONST, OR, VAR, XOR, BooleanANF, BooleanFunction
from PyPR.BooleanLogic.BooleanGF import BooleanGF
from PyPR.BooleanLogic.ChainingGeneration.Templates import fast_template

from PyPR.FeedbackFunctions import (
    CMPR,
    FCSR,
    MPR,
    CrossJoin,
    FeedbackFunction,
    Fibonacci,
    Galois,
    TFunction,
)
from PyPR.FeedbackRegister import FeedbackRegister

# ── Helpers ──────────────────────────────────────────────────────────────────

def _seq(reg, limit=20):
    return [int(state[0]) for state in reg.run(compiled=False, limit=limit)]

def _ff_json_round_trip(fn, seed, limit=20):
    """Serialize a FeedbackFunction to JSON and back; compare output sequences."""
    reg1 = FeedbackRegister(seed, fn)
    seq1 = _seq(reg1, limit)

    json_data = fn.to_JSON()
    fn2 = type(fn).from_JSON(json_data)

    reg2 = FeedbackRegister(seed, fn2)
    seq2 = _seq(reg2, limit)
    assert seq1 == seq2, f"{type(fn).__name__} JSON round-trip changed the output sequence"


# ── FeedbackFunction JSON round-trips ─────────────────────────────────────────

def test_mpr_json_round_trip():
    _ff_json_round_trip(MPR(5, "12"), seed=1)

def test_fibonacci_json_round_trip():
    _ff_json_round_trip(Fibonacci(5, "12"), seed=31)

def test_galois_json_round_trip():
    _ff_json_round_trip(Galois(5, "12"), seed=31)

def test_fcsr_json_round_trip():
    F = FCSR(3, 5)
    _ff_json_round_trip(F, seed=2**F.size - 1)

def test_cmpr_json_round_trip():
    random.seed(0)
    np.random.seed(0)  # chaining also draws from numpy's generator
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    M3 = MPR(3, "5")
    C = CMPR([M7, M5, M3])
    C.generateChaining(template=fast_template())
    _ff_json_round_trip(C, seed=2**C.size - 1)

def test_crossjoin_json_round_trip():
    random.seed(0)
    CJ = CrossJoin(7, "65")
    CJ.generateNonlinearity(maxAnds=3)
    _ff_json_round_trip(CJ, seed=2**7 - 1)

def test_tfunction_json_round_trip():
    T = TFunction(4)
    _ff_json_round_trip(T, seed=0, limit=16)


# ── FeedbackRegister JSON round-trip ──────────────────────────────────────────

def test_register_json_round_trip():
    """FeedbackRegister.to_JSON/from_JSON preserves the output sequence."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    seq1 = _seq(reg)

    reg.reset()
    json_data = reg.to_JSON()
    reg2 = FeedbackRegister.from_JSON(json_data)
    seq2 = _seq(reg2)

    assert seq1 == seq2

def test_register_json_round_trip_cmpr():
    """FeedbackRegister wrapping a CMPR serializes and restores correctly."""
    random.seed(0)
    np.random.seed(0)  # chaining also draws from numpy's generator
    M7 = MPR(7, "65")
    M5 = MPR(5, "12")
    C = CMPR([M7, M5])
    C.generateChaining(template=fast_template())
    seed = 2**C.size - 1

    reg = FeedbackRegister(seed, C)
    seq1 = _seq(reg, limit=30)

    reg.reset()
    json_data = reg.to_JSON()
    reg2 = FeedbackRegister.from_JSON(json_data)
    seq2 = _seq(reg2, limit=30)

    assert seq1 == seq2


# ── FeedbackRegister file round-trip ──────────────────────────────────────────

def test_register_file_round_trip():
    """to_file/from_file preserves the output sequence."""
    M = MPR(5, "12")
    reg = FeedbackRegister(1, M)
    seq1 = _seq(reg)
    reg.reset()

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = f.name
    try:
        reg.to_file(path)
        reg2 = FeedbackRegister.from_file(path)
        seq2 = _seq(reg2)
        assert seq1 == seq2
    finally:
        os.unlink(path)

def test_feedbackfunction_file_round_trip():
    """FeedbackFunction.to_file/from_file preserves the output sequence."""
    M = MPR(7, "65")
    reg1 = FeedbackRegister(1, M)
    seq1 = _seq(reg1, limit=30)

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = f.name
    try:
        M.to_file(path)
        M2 = MPR.from_file(path)
        reg2 = FeedbackRegister(1, M2)
        seq2 = _seq(reg2, limit=30)
        assert seq1 == seq2
    finally:
        os.unlink(path)


# ── Round-trip equality on the whole state space ─────────────────────────────

def test_json_round_trip_agrees_on_every_state():
    """A deserialized feedback function agrees with the original on all 2^n states.

    The sequence comparisons above follow one orbit from one seed, so they only
    exercise the bits that orbit happens to visit.  Comparing every bit's update
    function at every state is the full contract: it catches a bit whose DAG was
    dropped or rewired during serialization even when the bit-0 sequence is
    unaffected, which a single-orbit check cannot distinguish from success.

    Sizes are kept small because the comparison is exhaustive in the state space.
    """
    cases = [
        MPR(5, "12"),
        Fibonacci(5, "12"),
        Galois(5, "12"),
        FCSR(4, 13),
        TFunction(5),
        CMPR([MPR(3, [1, 1, 0, 1]), MPR(4, [1, 1, 0, 0, 1])]),
    ]

    for fn in cases:
        restored = type(fn).from_JSON(fn.to_JSON())
        size = len(fn)
        assert len(restored) == size, f"{type(fn).__name__}: round-trip changed the size"

        for value in range(2 ** size):
            state = [(value >> k) & 1 for k in range(size)]
            original_bits = [bit_fn.eval(state) for bit_fn in fn.fn_list]
            restored_bits = [bit_fn.eval(state) for bit_fn in restored.fn_list]
            assert original_bits == restored_bits, (
                f"{type(fn).__name__}: round-trip differs at state {state} "
                f"({original_bits} vs {restored_bits})"
            )


# ── Shared structure across objects ──────────────────────────────────────────

def test_structure_shared_between_objects_comes_back_shared():
    """One generate_JSON call gives each object one id, however many of the
    stored objects reach it, so parsing rebuilds a single shared object -- the
    same object, not equal copies."""
    shared = AND(VAR(0), VAR(1))
    a = XOR(shared, VAR(2))
    b = OR(shared, VAR(3), shared)

    a2, b2, shared2 = parse_JSON(json.loads(json.dumps(generate_JSON(a, b, shared))))

    assert a2.args[0] is shared2
    assert b2.args[0] is shared2
    assert b2.args[2] is shared2
    assert shared2 is not shared
    for bits in [[0, 0, 0, 0], [1, 1, 0, 1], [1, 1, 1, 0]]:
        assert a2.eval(bits) == a.eval(bits)
        assert b2.eval(bits) == b.eval(bits)

def test_a_register_and_its_function_stored_together_stay_linked():
    fn = MPR(5, "12")
    reg = FeedbackRegister(1, fn)
    reg2, fn2 = parse_JSON(json.loads(json.dumps(generate_JSON(reg, fn))))
    assert reg2.fn is fn2
    assert isinstance(fn2, MPR)

def test_the_same_object_passed_twice_is_stored_once():
    fn = XOR(VAR(0), VAR(1))
    data = generate_JSON(fn, fn)
    assert data["return order"][0] == data["return order"][1]
    first, second = parse_JSON(data)
    assert first is second


# ── Compatibility and the class registry ─────────────────────────────────────

def test_files_written_before_the_mixin_still_load():
    """The fixture was written by the per-class serialization this replaced.
    Class names are still stored as module.QualifiedName, so it loads, and every
    entry re-saves identically -- except that the CMPR entry no longer carries
    `blocks`, which has since become a cached property (derived from fn_list,
    so it is not written). The recomputed value matches the one stored."""
    stored = json.loads(Path("test/BooleanLogic/JSON Generation/test.json").read_text())
    (C,) = parse_JSON(stored)
    assert isinstance(C, CMPR)
    assert C.size == 12

    resaved = generate_JSON(C)
    assert resaved["objects"][:-1] == stored["objects"][:-1]
    old_entry, new_entry = stored["objects"][-1]["data"], resaved["objects"][-1]["data"]
    assert {k: v for k, v in old_entry.items() if k != "blocks"} == new_entry

    # drop the attribute the old file set, so the cached property recomputes it
    del C.__dict__["blocks"]
    assert C.blocks == old_entry["blocks"]

def test_every_serializable_class_is_registered_under_its_stored_name():
    for cls in [BooleanFunction, XOR, AND, VAR, CONST, BooleanANF, BooleanGF,
                FeedbackFunction, MPR, CMPR, FeedbackRegister]:
        name = f"{cls.__module__}.{cls.__qualname__}"
        assert Serializable._registry[name] is cls

def test_parsing_an_unregistered_class_names_it():
    data = XOR(VAR(0), VAR(1)).to_JSON()
    data["objects"][-1]["class"] = "somewhere.Unknown"
    with pytest.raises(TypeError, match="'somewhere.Unknown' not recognized"):
        parse_JSON(data)

def test_from_json_rejects_a_different_class():
    with pytest.raises(ValueError, match="not a subclass of PyPR.FeedbackFunctions.MPR.MPR"):
        MPR.from_JSON(Fibonacci(5, "12").to_JSON())
    # a superclass accepts its subclasses
    assert isinstance(FeedbackFunction.from_JSON(MPR(5, "12").to_JSON()), MPR)

def test_boolean_function_file_round_trip(tmp_path):
    fn = XOR(AND(VAR(0), VAR(1)), CONST(1))
    path = tmp_path / "fn.json"
    fn.to_file(str(path))
    fn2 = BooleanFunction.from_file(str(path))
    assert fn2.dense_str() == fn.dense_str()


# ── BooleanANF ───────────────────────────────────────────────────────────────

@pytest.mark.parametrize("anf", [
    BooleanANF(),                                    # the zero function
    BooleanANF([1]),                                 # the constant 1: one empty term
    BooleanANF([[0, 1], [2], []]),
    BooleanANF([["a", "b"], ["c"]]),                 # any hashable variable type
    BooleanANF([[(0, 1), (2,)], [(0, 1)]]),          # tuple variables, e.g. monomials
    BooleanANF([[1, "x", (2, 3)], [None]]),          # mixed types in one term
], ids=["zero", "one", "ints", "strings", "tuples", "mixed"])
def test_anf_round_trip(anf):
    restored = BooleanANF.from_JSON(json.loads(json.dumps(anf.to_JSON())))
    assert restored == anf

def test_anf_is_written_the_same_however_it_was_built():
    """The terms live in frozensets, whose order is arbitrary; the written form
    is ordered, so equal ANFs produce identical files."""
    first = BooleanANF([[2, 0], [1], [3, 1, 2], []])
    second = BooleanANF([[], [1, 2, 3], [1], [0, 2]])
    assert first == second
    assert json.dumps(first.to_JSON()) == json.dumps(second.to_JSON())
    assert first.to_JSON()["objects"][0]["data"] == {"terms": [[], [1], [0, 2], [1, 2, 3]]}

def test_anf_variable_that_json_cannot_hold_is_rejected():
    with pytest.raises(TypeError, match="variable of type frozenset"):
        BooleanANF([[frozenset({1})]]).to_JSON()

def test_anf_converted_from_a_function_round_trips():
    fn = XOR(AND(VAR(0), VAR(1)), OR(VAR(1), VAR(2)), CONST(1))
    anf = BooleanANF.from_BooleanFunction(fn)
    assert BooleanANF.from_JSON(anf.to_JSON()) == anf


# ── BooleanGF ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("fraction", [
    BooleanGF.zero(),
    BooleanGF.one(),
    BooleanGF.delay(),
    BooleanGF([1, 1, 0, 1], [1, 0, 1]),
    BooleanGF.from_seq([1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1]),
], ids=["zero", "one", "delay", "fraction", "from_seq"])
def test_gf_round_trip(fraction):
    restored = BooleanGF.from_JSON(json.loads(json.dumps(fraction.to_JSON())))
    assert restored == fraction

def test_gf_is_written_in_constructor_order():
    """Index i holds the coefficient of D^i, as the constructor reads it, so the
    written lists are exactly what built the fraction."""
    fraction = BooleanGF([1, 1, 0, 1], [1, 0, 1])   # (1 + D + D^3) / (1 + D^2)
    assert fraction.to_JSON()["objects"][0]["data"] == {
        "numerator": [1, 1, 0, 1], "denominator": [1, 0, 1],
    }

def test_gf_hash_agrees_with_equality():
    """BooleanGF defined __eq__ without __hash__, so it could not be a dict key --
    which is how serialization tracks objects -- and its __eq__ raised on any
    other type instead of returning NotImplemented."""
    a, b = BooleanGF([1, 1], [1]), BooleanGF([1, 1], [1])
    assert a == b
    assert hash(a) == hash(b)
    assert len({a, b, BooleanGF.one()}) == 2
    assert (BooleanGF.one() == 1) is False

def test_value_types_mix_with_functions_in_one_file():
    """Equal ANFs and fractions are one entry each, since they are immutable and
    compare by value; everything comes back in the order it was passed."""
    anf, fraction, fn = BooleanANF([[0], [1]]), BooleanGF.delay(), XOR(VAR(0), VAR(1))
    data = generate_JSON(anf, fraction, fn, BooleanANF([[1], [0]]))
    assert data["return order"][0] == data["return order"][3]

    anf2, fraction2, fn2, anf3 = parse_JSON(json.loads(json.dumps(data)))
    assert anf2 == anf
    assert anf3 is anf2
    assert fraction2 == fraction
    assert fn2.dense_str() == fn.dense_str()


# ── Objects held in node fields ──────────────────────────────────────────────

def test_const_holding_an_object_is_stored_by_reference():
    """remap_constants puts objects such as BooleanANFs into CONST nodes. They
    get entries of their own and the CONST refers to them by id, so one object
    used by several constants -- or also stored on its own -- comes back as one
    object. A CONST holding a plain bit is written exactly as before."""
    anf = BooleanANF([[0], [1]])
    f = XOR(CONST(anf), VAR(2), CONST(1))
    g = AND(CONST(anf), VAR(0))

    data = json.loads(json.dumps(generate_JSON(f, g, anf)))
    f2, g2, anf2 = parse_JSON(data)

    assert f2.args[0].value == anf
    assert f2.args[0].value is g2.args[0].value is anf2
    assert f2.args[2].value == 1

    const_entries = [e["data"] for e in data["objects"] if e["class"].endswith(".CONST")]
    assert {"args": [], "arg_limit": 0, "value": 1} in const_entries

def test_remapped_constants_round_trip():
    """The motivating case: a function whose constants were remapped to ANFs."""
    fn = XOR(AND(VAR(0), CONST(1)), CONST(0), VAR(1))
    remapped = fn.remap_constants([(0, BooleanANF([0])), (1, BooleanANF([1]))])
    restored = BooleanFunction.from_JSON(json.loads(json.dumps(remapped.to_JSON())))

    constants = [n for n in restored.postorder() if isinstance(n, CONST)]
    assert sorted(len(c.value) for c in constants) == [0, 1]   # ANFs of 0 and 1
    assert all(isinstance(c.value, BooleanANF) for c in constants)

def test_fused_node_round_trips_through_its_constructor():
    """A fusion node stores its template by reference and rebuilds the rest (its
    SAT skeleton) on parsing, so it both evaluates and encodes the same."""
    from PyPR.BooleanLogic.SAT import TseytinFuse
    fuse = TseytinFuse(OR(AND(VAR(0), VAR(1)), VAR(2)))
    fn = XOR(fuse(VAR(3), VAR(4), VAR(5)), VAR(0))

    restored = BooleanFunction.from_JSON(json.loads(json.dumps(fn.to_JSON())))

    assert type(restored.args[0]) is type(fn.args[0])
    for value in range(2 ** 6):
        bits = [(value >> i) & 1 for i in range(6)]
        assert restored.eval(bits) == fn.eval(bits)
    assert restored.tseytin()[0] == fn.tseytin()[0]
