"""Tests for ResolventSolving: field_eye and field_invert over BooleanGF matrices.

field_eye creates an identity matrix; field_invert solves M^{-1} via Gaussian
elimination over the BooleanGF field.  The core invariant is M * M^{-1} = I.
"""
import pytest
import numpy as np

from PyPR.BooleanLogic.BooleanGF import BooleanGF
from PyPR.Tools.ResolventSolving import field_eye, field_invert


# ── field_eye ─────────────────────────────────────────────────────────────────

def test_field_eye_shape():
    """field_eye(field, n) returns an n×n matrix."""
    I = field_eye(BooleanGF, 4)
    assert I.shape == (4, 4)

def test_field_eye_diagonal_is_one():
    """Diagonal entries of field_eye are BooleanGF.one()."""
    n = 3
    I = field_eye(BooleanGF, n)
    for i in range(n):
        assert I[i, i].simplify() == BooleanGF.one().simplify(), (
            f"Diagonal entry [{i},{i}] is not one: {I[i, i]}"
        )

def test_field_eye_off_diagonal_is_zero():
    """Off-diagonal entries of field_eye are BooleanGF.zero()."""
    n = 3
    I = field_eye(BooleanGF, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                assert I[i, j].simplify() == BooleanGF.zero().simplify(), (
                    f"Off-diagonal entry [{i},{j}] is not zero: {I[i, j]}"
                )

def test_field_eye_size_one():
    """field_eye of size 1 is the 1×1 identity."""
    I = field_eye(BooleanGF, 1)
    assert I.shape == (1, 1)
    assert I[0, 0].simplify() == BooleanGF.one().simplify()


# ── field_invert ──────────────────────────────────────────────────────────────

def test_field_invert_2x2_satisfies_mv_identity():
    """M * field_invert(M) == I for a known 2×2 invertible BooleanGF matrix.

    Over GF(2), [[1,1],[0,1]]^2 = [[1,0],[0,1]], so this matrix is its own
    inverse.  We compare via simplify() because BooleanGF.__eq__ compares
    num/den directly without reducing to lowest terms.
    """
    # [[1, 1], [0, 1]] as constant BooleanGF elements (numerator/denominator lists)
    entries = [
        BooleanGF([1], [1]), BooleanGF([1], [1]),
        BooleanGF([0], [1]), BooleanGF([1], [1]),
    ]
    M = np.asarray(entries, dtype=BooleanGF).reshape(2, 2)
    Minv = field_invert(BooleanGF, M)
    product = M @ Minv
    # Verify M * M^-1 = I entry-by-entry after reduction to lowest terms
    n = 2
    for i in range(n):
        for j in range(n):
            expected = BooleanGF.one() if i == j else BooleanGF.zero()
            assert product[i, j].simplify() == expected.simplify(), (
                f"M * M^-1 not identity at [{i},{j}]: got {product[i, j]}"
            )

def test_field_invert_identity_matrix():
    """The inverse of the identity matrix is itself."""
    n = 3
    I = field_eye(BooleanGF, n)
    Iinv = field_invert(BooleanGF, I)
    # simplify() needed: BooleanGF.__eq__ is syntactic (num/den), not semantic
    for i in range(n):
        for j in range(n):
            expected = BooleanGF.one() if i == j else BooleanGF.zero()
            assert Iinv[i, j].simplify() == expected.simplify(), (
                f"Inverse of identity not identity at [{i},{j}]: got {Iinv[i, j]}"
            )

def test_field_invert_raises_for_nonsquare():
    """field_invert raises ValueError for non-square input."""
    non_square = np.asarray(
        [BooleanGF.one(), BooleanGF.zero(), BooleanGF.one()],
        dtype=BooleanGF
    ).reshape(1, 3)
    with pytest.raises(ValueError):
        field_invert(BooleanGF, non_square)

def test_field_invert_3x3_round_trip():
    """M * M^-1 == I for a 3×3 upper-triangular GF(2) matrix.

    Upper-triangular matrices with ones on the diagonal are always invertible
    over any field, so [[1,1,0],[0,1,1],[0,0,1]] is a safe test case.
    simplify() is required before equality: BooleanGF stores num/den symbolically
    and does not auto-reduce.
    """
    # [[1,1,0],[0,1,1],[0,0,1]] as BooleanGF constants
    rows = [
        [([1], [1]), ([1], [1]), ([0], [1])],
        [([0], [1]), ([1], [1]), ([1], [1])],
        [([0], [1]), ([0], [1]), ([1], [1])],
    ]
    entries = [BooleanGF(num, den) for row in rows for num, den in row]
    M = np.asarray(entries, dtype=BooleanGF).reshape(3, 3)
    Minv = field_invert(BooleanGF, M)
    product = M @ Minv
    n = 3
    for i in range(n):
        for j in range(n):
            expected = BooleanGF.one() if i == j else BooleanGF.zero()
            assert product[i, j].simplify() == expected.simplify(), (
                f"3×3 M * M^-1 not identity at [{i},{j}]: got {product[i, j]}"
            )
