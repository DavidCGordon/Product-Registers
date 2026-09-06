from PyPR.BooleanLogic import BooleanFunction

from PyPR.Cryptanalysis.Components.EquationStores.IndexedEqStore import IndexedEqStore
from PyPR.Cryptanalysis.Components.Adapters.equation_repr import (
    extract_monomials,
    boolean_function_to_coef_vector,
)
import numpy as np

import numba
u8 = numba.types.uint8
u64 = numba.types.uint64
b1 = numba.types.b1

from typing import Any

@numba.njit(numba.types.Tuple((b1,u64))(u8[:,:],u8[:,:],u8[:],u8[:],u8[:],u64,u64))
def _LU_reduction_consistent(
    upper: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    lower: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    extra_constants: np.ndarray[tuple[int],np.dtype[np.uint8]],
    solved_for: np.ndarray[tuple[int],np.dtype[np.uint8]],
    coef_vector: np.ndarray[tuple[int],np.dtype[np.uint8]],
    num_vars: int,
    const_idx: int
) -> tuple[bool, int]:
    """This function performs LU reduction, and inserts the results into the corresponding matrix.

    This variant is for the consistent version, which has to treat constants differently, in order
    to perform consistency checks. The index maps are required to have an index for constants, and
    the _extra_constants parameter of the equation store is used for additional bookkeeping (and thus,
    will be nonzero). The information is structured and can be used, but is complicated, and it's
    advised to avoid using it unless you know what you are doing. The structure is described below:

    When decomposing, there are two main ways to look at constants:
     - **Associated/Augmented constants:** the constants are stored as if they were independent of the \
        reduction process; when the upper/lower equation pair is placed, the constant is slotted in \
        unmodified, and will be correct when the equation is reconstructed.
     - **Accumulated constants:** the constants are modified when reducing, and the constant in a row \
        is equal to the constant of the reduced equation, rather than the constant of the reconstructed \
        equation. when reducing, in order to maintain accuracy, you must use other accumulated constants \
        so the effects of reduction tend to "fold" into each other.

    There is a natural correspondence between view of the constants to the L/U decomposition:
     - **Augmented constants correspond to 1's in the constant column of L:** If the equation E can be \
        reconstructed as the sum of some reduced vectors from U (`u_1, ... u_n`), then changing the \
        coefficient corresponding to the constant equation (i.e. `1=1`) simply changes the constant \
        added to reconstructed equation (the constant is simply added in at the end).
     - **Accumulated constants correspond to 1's in the constant column of U:** These constants are the \
        values of the rows the are present in (zeroing out each row individually, and thus zeroing out \
        the sum when many such `u_1, ... u_n` are added together). When reducing equations, \
        the naive strategy of treating the constant as a variable will naturally accumulate the \
        constants in the upper matrix (and in the inconsistent case anything after the constant \
        column will have its constants cancelled by the reduction prcess, so nothing about the lower \
        matrix needs to be considered / handled specially). This is the constant you need to know \
        to check if an equation is consistent (if the rest of the equation reduces, the accumulated \
        constant must be `0` in order to be consistent).

    When the constants are stored at index `0` (or indpendently), it makes most sense to store them as
    augmented constants: The constants are immediately stored in L, and cancelled, letting the rest of
    the reduction proceed without constants. Because index 0 is a large column for L, this also does not
    break it's triangularity, but storing the augmented constants at any other index might.

    When the constants are stored at the index `num_vars-1` (the last index), it makes the most sense
    to store the constants as accumulated constants: they will naturally accumulate in that column over
    the course of the reduction process, and are checked at the end. Additionally, because  the last
    index is a large column for U, this also does not break it's triangularity, but storing the accumulated
    constants at any earlier index might.

    For consistency checking, the second option would be much more convenient because we get the augmented
    constants, but there are some problems: We strongly want to avoid modifying the input index maps (to
    avoid breaking the usage contract / surprising the user) and the index maps can also be dynamic (so
    always using the "last" index may not even make sense since it changes). We also strongly want to
    maintain the triangularity of L and U, to avoid nasty edge cases and surprising results, which prevents
    us from just storing the accumulated constants in U at an arbitrary column. This means that we must
    handle the slightly more complicated case which interpolates between the two cases above, and we
    understand the reduction behavior for for an arbitrary `const_idx`:

    The reduction proceeds as normal from index 0 to const_idx, the accumulated constants are computed
    naturally, and stored in the upper matrix. Any equation which is reduced past the const_idx will
    have its augmented constant stored in the lower matrix when it is placed. However, this constant
    will have already been partially reduced by the upper portion of the reduction process, thus what is
    stored in the lower matrix is not the true augmented constant but the cached 'derived constant'.
    However, we need the accumulated constant for consistency checking (and we need the other accumulated
    constants in order to calculate these); this is where _extra_consts comes in. We use the portion
    after `const_idx` to store the accumulated constants we need to keep computing accumulated constants.
    For symmetry and to store additional information, we use the portion before `const_idx`to store the
    augmenting constants, which are not captured in U.

    Although understanding the verbal explanation is the best way to interact with the `_extra_constants`
    vector, the following true statements can be useful examples/shortcuts:
     - `L[:const_idx,   const_idx]` = zeroes
     - `L[const_idx+1:, const_idx]` = derived constant (augmented constant for the partially reduced equation)
     - `U[:const_idx,   const_idx]` = accumulated constant
     - `U[const_idx+1:, const_idx]` = zeroes
     - `_extra_constants[:const_idx,   const_idx]` = augmented constant
     - `_extra_constants[const_idx+1:, const_idx]` = accumulated constant

    :param upper: The Upper matrix of the Equation store. The reduced coefficients will be
        inserted here
    :type upper: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param lower: The Lower matrix of the Equation store. The reduction coefficients will be
        inserted here
    :type lower: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param extra_constants: The vector storing additional information to help with consistency
        checking
    :type extra_constants: np.ndarray[tuple[int],np.dtype[np.uint8]]
    :param solved_for: The vector describing which variables have an associated equation.
    :type solved_for: np.ndarray[tuple[int],np.dtype[np.uint8]]
    :param coef_vector: The vector of coefficients for the equation being inserted
    :type coef_vector: np.ndarray[tuple[int],np.dtype[np.uint8]]
    :param num_vars: The number of variables the equation store is currently tracking (which
        might differ from the length of the coefficient vector). This allows the loop to be
        shorter for dynamic equation stores.
    :type num_vars: int
    :param const_idx: The index which corresponds to constant terms in the equation store
    :type const_idx: int
    :raises ValueError: If the inserted equation is inconsistent
    :return: A tuple with a boolean indicating whether or not the equation was inserted,
        and an int indicating at which index the equation was inserted (if not inserted,
        the index will returned but has no meaning and shouldn't be used)
    :rtype: tuple[bool, int]
    """
    linearly_independent = False
    coef_vector = coef_vector.copy()
    modification_vector = np.zeros_like(coef_vector)

    augmented_constant = coef_vector[const_idx]
    for idx in range(const_idx):
        # skip coefficients which don't need to be cancelled:
        if coef_vector[idx] == 0:
            continue

        modification_vector[idx] = 1

        if solved_for[idx]:
            # update coefficient vectors
            for i in range(idx,num_vars):
                coef_vector[i] ^= upper[idx,i]
        else:
            linearly_independent = True
            solved_for[idx] = 1
            upper[idx] = coef_vector
            lower[idx] = modification_vector

            # top of extra constants stores the augmented constants
            # accumulated constants are naturally computed in the upper matrix
            extra_constants[idx] = augmented_constant
            break

    # cancel constants as if they were at the end, to enable checking consistency
    if linearly_independent:
        return linearly_independent, idx

    derived_constant = coef_vector[const_idx]
    for idx in range(const_idx+1,len(coef_vector)):
        # skip coefficients which don't need to be cancelled:
        if coef_vector[idx] == 0:
            continue

        modification_vector[idx] = 1

        if solved_for[idx]:
            coef_vector[const_idx] ^= extra_constants[idx]
            for i in range(idx,num_vars):
                coef_vector[i] ^= upper[idx,i]

        else:
            linearly_independent = True
            solved_for[idx] = 1
            upper[idx] = coef_vector
            lower[idx] = modification_vector

            # bottom of extra constants stores the augmented constants
            # the derived constants are stored in the lower matrix to allow proper reconstruction.
            extra_constants[idx] = coef_vector[const_idx]
            lower[idx, const_idx] = derived_constant
            upper[idx, const_idx] = 0
            break

    if (
        not linearly_independent and
        coef_vector[const_idx] == 1   # accumulated constant == 1
    ):
        raise ValueError("Inconsistent!")

    return linearly_independent, idx

@numba.njit(numba.types.Tuple((b1,u64))(u8[:,:],u8[:,:],u8[:],u8[:],u64))
def _LU_reduction(
    upper: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    lower: np.ndarray[tuple[int,int],np.dtype[np.uint8]],
    solved_for: np.ndarray[tuple[int],np.dtype[np.uint8]],
    coef_vector: np.ndarray[tuple[int],np.dtype[np.uint8]],
    num_vars: int,
) -> tuple[bool, int]:
    """This function performs LU reduction, and inserts the results into the corresponding matrix.

    This variant is for the inconsistent version. CONST(1) is treated the same as any other monomial,
    and is simply decomposed into L and U (extra constants are not used, and remain all zero).

    :param upper: The Upper matrix of the Equation store. The reduced coefficients will be
        inserted here
    :type upper: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param lower: The Lower matrix of the Equation store. The reduction coefficients will be
        inserted here
    :type lower: np.ndarray[tuple[int,int],np.dtype[np.uint8]]
    :param solved_for: The vector describing which variables have an associated equation.
    :type solved_for: np.ndarray[tuple[int],np.dtype[np.uint8]]
    :param coef_vector: The vector of coefficients for the equation being inserted
    :type coef_vector: np.ndarray[tuple[int],np.dtype[np.uint8]]
    :param num_vars: The number of variables the equation store is currently tracking (which
        might differ from the length of the coefficient vector). This allows the loop to be
        shorter for dynamic equation stores.
    :type num_vars: int
    :return: A tuple with a boolean indicating whether or not the equation was inserted,
        and an int indicating at which index the equation was inserted (if not inserted,
        the index will returned but has no meaning and shouldn't be used)
    :rtype: tuple[bool, int]
    """
    linearly_independent = False
    coef_vector = coef_vector.copy()
    modification_vector = np.zeros_like(coef_vector)
    for idx in range(len(coef_vector)):
        if coef_vector[idx] == 1:
            modification_vector[idx] = 1
            if solved_for[idx]:
                for i in range(idx,num_vars):
                    coef_vector[i] ^= upper[idx,i]
            else:
                linearly_independent = True
                solved_for[idx] = 1
                upper[idx] = coef_vector
                lower[idx] = modification_vector
                break

    return linearly_independent, idx

class LUEqStore(IndexedEqStore):
    """Incremental LU decomposition store.

    Performs LU reduction on each inserted equation, rejecting linearly
    dependent information. ``solved_for`` tracks which columns have pivots
    (the rank), but does NOT contain actual variable values — run
    LU_Solver for back-substitution.
    """

    def _enforce_consistency(self):
        # set consistency flag:
        self.consistent = True

        if tuple() not in self.comb_to_idx:
            if self.dynamic:
                # For Dynamic stores, we can just update the known 
                # monomials to support the declaration that CONST(1)=1
                self._update_known_monomials(tuple())
            else:
                # For static stores without a const idx, we cant support that declaration:
                # thus, equation are consistent only when all inserted equations have
                # zero constant term. Any linear combination also has zero constant,
                # so 0=0 — never a contradiction. Insertion validates this at equation-
                # insertion time (see _update_known_monomials).
                
                # Because we can't support the declaration, there is nothing to do here
                return
                
        # Constant column present: mark CONST(1) = 1 as an axiom so that
        # _LU_reduction_consistent can detect contradictions via accumulated constants.
        if not self.solved_for[self.comb_to_idx[tuple()]]:
            self.solved_for[self.comb_to_idx[tuple()]] = 1
            self.num_eqs += 1
        else:
            raise ValueError("Haven't implemented the case where const already has a definition.")

    def __init__(
        self,
        comb_to_idx: dict[tuple[int,...], int] | None = None,
        consistent: bool = False
    ):
        super().__init__(comb_to_idx, consistent=False)
        self.filtering = True

        if self.dynamic:
            self.lower_matrix = np.eye(256, dtype = 'uint8')
            self.upper_matrix = np.eye(256, dtype = 'uint8')
            self.solved_for = np.zeros([256], dtype = 'uint8')
            self._extra_constants = np.zeros([256], dtype = 'uint8')
        else:
            self.lower_matrix = np.eye(self.num_vars, dtype = 'uint8')
            self.upper_matrix = np.eye(self.num_vars, dtype = 'uint8')
            self.solved_for = np.zeros([self.num_vars], dtype = 'uint8')
            self._extra_constants = np.zeros([self.num_vars], dtype = 'uint8')

        if consistent:
            self._enforce_consistency()

    def _expand_storage(self):
        """Double the capacity of L, U, solved_for, and _extra_constants."""
        if self.num_vars < self.lower_matrix.shape[1]:
            return

        new_size = self.num_vars * 2
        new_upper_matrix = np.eye(new_size, dtype=np.uint8)
        new_lower_matrix = np.eye(new_size, dtype=np.uint8)
        new_solved_for = np.zeros([new_size], dtype=np.uint8)
        new_extra_consts = np.zeros([new_size], dtype=np.uint8)

        new_upper_matrix[:self.num_vars, :self.num_vars] = self.upper_matrix
        new_lower_matrix[:self.num_vars, :self.num_vars] = self.lower_matrix
        new_solved_for[:self.num_vars] = self.solved_for
        new_extra_consts[:self.num_vars] = self._extra_constants

        self.upper_matrix = new_upper_matrix
        self.lower_matrix = new_lower_matrix
        self.solved_for = new_solved_for
        self._extra_constants = new_extra_consts

    def insert_equation(
        self,
        equation: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]],
        identifier: Any = None,
        translate_ANF: bool = True
    ):
        """Insert an equation into the LU store.

        The equation will be decomposed, and checked for linear dependence (filtering out
        informatuion which is redundant, allowing addition but not scaling). Returns whether
        or not the equation was successfully added or not (i.e. whether or not it was
        linearly independent). Has variants for both consistent and inconsistent systems.

        :param equation: The equation to insert, as a BooleanFunction(or numpy array, if the
            equation store is static)
        :type equation: BooleanFunction | np.ndarray[tuple[int],np.dtype[np.uint8]]
        :param identifier: Additional metadata about an equation such as a clock time,
            or cube used to generate the equation, defaults to None
        :type identifier: Any, optional
        :param translate_ANF: Whether or not to translate to ANF internally, defaults to True
            for safety and flexibility, but can be disabled for a small performance boost if
            you know your equations will be in ANF already.
        :type translate_ANF: bool, optional
        :raises ValueError: If passed bad input (i.e. wrong type/shape, etc.) for the equation
        :raises ValueError: If the store is consistent and you insert an inconsistent equation
        :raises TypeError: If the given equation is not of the right type.
        :return: Whether or not the equation was inserted correctly (i.e. whether or not it was
            linearly independent)
        :rtype: bool
        """
        if type(equation) == np.ndarray:
            if self.dynamic:
                raise ValueError(
                    "Cannot insert an ndarray into a dynamic equation store. " +
                    "An index mapping must be provided in order to enable this feature."
                )
            if equation.shape != (self.num_vars,):
                raise ValueError(
                    f"Equation store has {self.num_vars} variables, " +
                    f"but the provided ndarray has shape {equation.shape}."
                )
            coef_vector = equation

        elif isinstance(equation, BooleanFunction):
            for comb in extract_monomials(equation, translate_ANF):
                if comb not in self.comb_to_idx:
                    self._update_known_monomials(comb)

            coef_vector = np.zeros(self.lower_matrix.shape[1], dtype=np.uint8)
            coef_vector[:self.num_vars] = boolean_function_to_coef_vector(
                equation, self.comb_to_idx, self.num_vars, translate_ANF
            )
        else:
            raise TypeError(
                f"Equation must be a BooleanFunction or ndarray, not {type(equation)}."
            )

        if self.consistent and tuple() in self.comb_to_idx:
            linearly_independent, insertion_idx = _LU_reduction_consistent(
                self.upper_matrix,
                self.lower_matrix,
                self._extra_constants,
                self.solved_for,
                coef_vector,
                self.num_vars,
                self.comb_to_idx[tuple()]
            )
        else:
            linearly_independent, insertion_idx =  _LU_reduction(
                self.upper_matrix,
                self.lower_matrix,
                self.solved_for,
                coef_vector,
                self.num_vars,
            )

        if linearly_independent:
            self.equation_ids[insertion_idx] = identifier
            self.num_eqs += 1

        return linearly_independent

    @property
    def rank(self):
        return self.num_eqs

    @property
    def is_determined(self):
        return self.rank == self.num_vars
