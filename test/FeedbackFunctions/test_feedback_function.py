"""Tests for behaviour shared by every FeedbackFunction through the base class:
copying, the linearity check, and VHDL export."""
import copy
import random

import numpy as np
import pytest

from PyPR.BooleanLogic import AND, CONST, VAR, XOR
from PyPR.BooleanLogic.ChainingGeneration.Templates import fast_template

from PyPR.FeedbackFunctions import CMPR, MPR, FeedbackFunction, Fibonacci
from PyPR.FeedbackRegister import FeedbackRegister

# copy.copy, copy.deepcopy and .copy() are one operation: always a deep copy
COPIERS = [copy.copy, copy.deepcopy, lambda obj: obj.copy()]
COPIER_IDS = ["copy.copy", "copy.deepcopy", ".copy()"]


@pytest.mark.parametrize("copier", COPIERS, ids=COPIER_IDS)
def test_copy_leaves_original_intact(copier):
    """Copying used to share the original's __dict__, so the copy's reset of
    _compiled and fn_list landed on the original too."""
    C = CMPR([MPR(5, "12"), MPR(3, "5")])
    # chaining draws from both random and numpy's generator, so seed both
    random.seed(0)
    np.random.seed(0)
    C.generateChaining(template=fast_template())
    C.compile()
    original_fn_list = C.fn_list
    original_nodes = list(C.fn_list)

    C2 = copier(C)

    assert type(C2) is CMPR
    assert C2.__dict__ is not C.__dict__
    assert C.fn_list is original_fn_list
    assert C.fn_list == original_nodes
    assert C._compiled is not None

    # both run compiled: the copy shares the (immutable) numba dispatchers
    FeedbackRegister(1, C).clock()
    FeedbackRegister(1, C2).clock()

@pytest.mark.parametrize("copier", COPIERS, ids=COPIER_IDS)
def test_copy_is_independent(copier):
    """Mutable attributes (polynomial lists, fn_list, cached arrays) are not shared."""
    C = CMPR([MPR(5, "12"), MPR(3, "5")])
    C2 = copier(C)

    assert C2.fn_list is not C.fn_list
    assert all(a is not b for a, b in zip(C2.fn_list, C.fn_list, strict=True))
    assert C2.primitive_polynomials == C.primitive_polynomials

    C2.primitive_polynomials.append([1])
    C2.divisions.append(99)
    assert C2.primitive_polynomials != C.primitive_polynomials
    assert 99 not in C.divisions

@pytest.mark.parametrize("copier", COPIERS, ids=COPIER_IDS)
def test_copy_preserves_sharing_across_bits_and_attributes(copier):
    """One memo threads through the whole copy, so a node reached from two bits,
    or from an extra attribute, is copied once and stays shared in the copy."""
    shared = AND(VAR(0), VAR(1))
    F = FeedbackFunction([XOR(shared, VAR(2)), XOR(shared, VAR(0)), VAR(1)])
    F.extra = {"node": shared, "array": np.array([1, 2, 3])}

    F2 = copier(F)
    shared2 = F2.fn_list[0].args[0]

    assert shared2 is not shared
    assert F2.fn_list[1].args[0] is shared2
    assert F2.extra["node"] is shared2
    assert F2.extra["array"] is not F.extra["array"]
    assert F2.extra["array"].tolist() == [1, 2, 3]

@pytest.mark.parametrize("copier", COPIERS, ids=COPIER_IDS)
def test_copy_simulates_identically(copier):
    C = CMPR([MPR(5, "12"), MPR(3, "5")])
    # chaining draws from both random and numpy's generator, so seed both
    random.seed(0)
    np.random.seed(0)
    C.generateChaining(template=fast_template())
    C2 = copier(C)

    reg, reg2 = FeedbackRegister(1, C), FeedbackRegister(1, C2)
    states = [s[:].tolist() for s in reg.run(40, compiled=False)]
    states2 = [s[:].tolist() for s in reg2.run(40, compiled=False)]
    assert states == states2

@pytest.mark.parametrize("copier", COPIERS, ids=COPIER_IDS)
def test_register_copy(copier):
    """A register copy owns its fn and state arrays, and a compiled register
    copies to one that still runs compiled."""
    reg = FeedbackRegister(5, Fibonacci(5, "12"), compile=True)
    reg.clock()

    reg2 = copier(reg)

    assert reg2.fn is not reg.fn
    assert reg2._state is not reg._state
    assert reg2._seed is not reg._seed
    assert reg2[:].tolist() == reg[:].tolist()

    reg2.clock()
    assert reg2[:].tolist() != reg[:].tolist()
    reg.clock()
    assert reg2[:].tolist() == reg[:].tolist()

    reg2.reset()
    assert reg2[:].tolist() == reg._seed.tolist()

def test_register_deepcopy_threads_memo_into_fn():
    """Copying a register together with its function in one call keeps them linked."""
    F = Fibonacci(5, "12")
    reg = FeedbackRegister(1, F)
    reg2, F2 = copy.deepcopy((reg, F))
    assert reg2.fn is F2
    assert F2 is not F

def test_is_linear():
    """isLinear used to return from inside its loop after the first gate type.

    It is decided on algebraic degree, not gate names: an LFSR built from its
    ANF still contains AND nodes (one per monomial, single-variable included),
    so a gate-name check would call it nonlinear."""
    F = Fibonacci(5, "12")
    assert "AND" in F.gateSummary()
    assert F.isLinear()
    assert F.isLinear(allowAffine=True)
    assert MPR(5, "12").isLinear()

    # chaining multiplies bits of neighbouring MPR blocks, so the result has degree > 1
    C = CMPR([MPR(5, "12"), MPR(3, "5")])
    # chaining draws from both random and numpy's generator, so seed both
    random.seed(0)
    np.random.seed(0)
    C.generateChaining(template=fast_template())
    assert not C.isLinear()
    assert not C.isLinear(allowAffine=True)

def test_is_linear_affine_split():
    """A constant monomial makes an update affine but not linear: it no longer
    fixes the zero state."""
    affine = FeedbackFunction([VAR(1), XOR(VAR(0), CONST(1))])
    assert not affine.isLinear()
    assert affine.isLinear(allowAffine=True)

    # x0 + 1 + 1 = x0: the constants cancel in the ANF, so this one is linear
    cancelled = FeedbackFunction([VAR(1), XOR(VAR(0), CONST(1), CONST(1))])
    assert cancelled.isLinear()

def test_write_vhdl_output_uses_declared_signal(tmp_path):
    path = tmp_path / "fpr.vhd"
    Fibonacci(5, "12").write_VHDL(str(path))
    vhdl = path.read_text()

    assert "signal curr_state, next_state" in vhdl
    assert "output <= curr_state;" in vhdl
    assert "currstate" not in vhdl
