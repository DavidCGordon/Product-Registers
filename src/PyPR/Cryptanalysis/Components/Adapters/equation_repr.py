"""Single-equation representation converters.

Bridges between the three equation representations:
  - Coefficient vector (np.ndarray uint8) — used by EqStore, LUEqStore
  - BooleanANF (frozenset of frozensets) — used by GroebnerEqStore, SymbolicEqStore
  - BooleanFunction (DAG of gates) — used by equation generators
"""
import numpy as np
from numpy.typing import NDArray

from PyPR.BooleanLogic import BooleanFunction, BooleanANF, XOR, AND, VAR, CONST


def extract_monomials(
    equation: BooleanFunction,
    translate_ANF: bool = True,
) -> list[tuple[int, ...]]:
    """Extract monomial tuples from a BooleanFunction.

    Each monomial is a sorted tuple of variable indices.
    The constant term (CONST(1)) is represented as the empty tuple ``()``.
    CONST(0) terms are dropped.

    :param equation: The boolean function to parse.
    :type equation: BooleanFunction
    :param translate_ANF: Whether to convert to ANF first.
    :type translate_ANF: bool
    :return: List of monomial tuples present in the equation.
    :rtype: list[tuple[int, ...]]
    """
    if translate_ANF:
        equation_anf = equation.translate_ANF()
    else:
        equation_anf = equation

    monomials = []
    for term in equation_anf.args:
        if type(term) != CONST:
            comb = tuple(sorted([var.index for var in term.args]))
        elif term.value == 1:
            comb = tuple()
        else:
            continue
        monomials.append(comb)
    return monomials

def boolean_function_to_coef_vector(
    equation: BooleanFunction,
    comb_to_idx: dict[tuple[int, ...], int],
    num_vars: int,
    translate_ANF: bool = True,
) -> NDArray[np.uint8]:
    """Parse a BooleanFunction into a coefficient vector over an index map.

    All monomials in the equation must already exist in ``comb_to_idx``.
    For dynamic stores, call ``_update_known_monomials`` for any new
    monomials before calling this function.

    :param equation: The boolean function to convert.
    :type equation: BooleanFunction
    :param comb_to_idx: Mapping from monomial tuple to column index.
    :type comb_to_idx: dict[tuple[int, ...], int]
    :param num_vars: Number of variables (length of the output vector).
    :type num_vars: int
    :param translate_ANF: Whether to convert to ANF first.
    :type translate_ANF: bool
    :return: Coefficient vector of length ``num_vars``.
    :rtype: NDArray[np.uint8]
    """
    monomials = extract_monomials(equation, translate_ANF)
    coef_vector = np.zeros(num_vars, dtype=np.uint8)
    for comb in monomials:
        coef_vector[comb_to_idx[comb]] = 1
    return coef_vector

def coef_vector_to_anf(
    coef_vector: NDArray[np.uint8],
    idx_to_comb: dict[int, tuple[int, ...]],
) -> BooleanANF:
    """Convert a coefficient vector to a BooleanANF.

    :param coef_vector: Binary coefficient vector (uint8, entries 0 or 1).
    :type coef_vector: NDArray[np.uint8]
    :param idx_to_comb: Mapping from column index to monomial tuple.
    :type idx_to_comb: dict[int, tuple[int, ...]]
    :return: The corresponding BooleanANF.
    :rtype: BooleanANF
    """
    terms = set()
    for idx in range(len(coef_vector)):
        if coef_vector[idx]:
            terms.add(frozenset(idx_to_comb[idx]))
    return BooleanANF(frozenset(terms), fast_init=True)

def coef_vector_to_boolean_function(
    coef_vector: NDArray[np.uint8],
    idx_to_comb: dict[int, tuple[int, ...]],
) -> BooleanFunction:
    """Convert a coefficient vector to a BooleanFunction in ANF form.

    This is the inverse of the BoolFunc→coef_vector parsing in
    ``insert_equation``. Enables CubeEqGenerator output to feed
    into symbolic stores.

    :param coef_vector: Binary coefficient vector (uint8, entries 0 or 1).
    :type coef_vector: NDArray[np.uint8]
    :param idx_to_comb: Mapping from column index to monomial tuple.
    :type idx_to_comb: dict[int, tuple[int, ...]]
    :return: A BooleanFunction in ANF (top-level XOR of ANDs).
    :rtype: BooleanFunction
    """
    return coef_vector_to_anf(coef_vector, idx_to_comb).to_BooleanFunction()
