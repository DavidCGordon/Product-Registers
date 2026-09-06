from typing import Any, Iterable
import numpy as np

from PyPR import FeedbackRegister
from PyPR.FeedbackFunctions import FeedbackFunction, Fibonacci, CMPR
from PyPR.BooleanLogic import BooleanFunction, BooleanANF, AND, XOR, VAR, CONST

from PyPR.Tools.RootCounting.MonomialProfile import TermSet, MonomialProfile
from PyPR.Tools.RootCounting.JordanSet import JordanSet
from PyPR.Tools.RootCounting.RootExpression import RootExpression

from random import randint, sample


class CrossJoin(FeedbackFunction):
    """A crossjoin nonlinear feedback shift register.

    Implements the crossjoin construction of Dubrova, which begins with a
    Fibonacci LFSR and introduces nonlinear (AND) terms that each appear
    at two distinct shifted positions in the feedback. The paired
    placement ensures that the nonlinear contributions cancel in a
    specific algebraic sense, preserving the period of the underlying
    LFSR while increasing linear complexity.

    The parameter tau delimits the nonlinear region: bits at index tau and
    above may carry nonlinear terms, while bits below tau remain purely
    linear. Currently only supports ANF (algebraic normal form) terms.

    :ivar size: The number of bits in the register.
    :vartype size: int
    :ivar primitive_polynomial: The coefficient list of the base LFSR's
        primitive polynomial, with index i holding the coefficient of x^i.
    :vartype primitive_polynomial: list[int]
    :ivar tau: The index of the lowest bit that may carry nonlinear terms.
    :vartype tau: int
    """

    primitive_polynomial: list[int]
    tau: int

    def __init__(self,
        size: int,
        primitive_poly: str | list[int]
    ) -> None:
        """Construct a crossjoin register from a Fibonacci LFSR base.

        Builds the linear Fibonacci shift from the given primitive
        polynomial. Nonlinear terms are not added until
        :meth:`generateNonlinearity` is called.

        :param size: The number of bits in the register.
        :type size: int
        :param primitive_poly: The primitive polynomial of the base LFSR.
            May be given as a coefficient list of length n+1 (with index i
            holding the coefficient of x^i) or as a Koopman hex string.
        :type primitive_poly: str | list[int]
        """
        # convert koopman string into polynomial:
        if isinstance(primitive_poly, str):
            primitive_poly = [int(x) for x in format(int(primitive_poly,16), f"0>{size}b")] + [1]
            
        self.primitive_polynomial = primitive_poly

        # same as fibonacci ANF generation:
        top_fn = [[size-idx] for (idx, t) in enumerate(primitive_poly) if t == 1][::-1][:-1]
        fn_list = [[[(i+1)%size]] for i in range(size-1)] + [top_fn]
        self.fn_list: list[BooleanFunction] = [
            XOR(BooleanFunction.from_ANF(bitFn),)
            for bitFn in fn_list
        ]

        self.size = size
        self.tau = self.size-1

    def shiftTerms(self,
        terms: list[Any],
        idxA: int,
        idxB: int
    ) -> None:
        """Move nonlinear terms from one bit position to another.

        Each term is shifted by (idxB - idxA) positions: all VAR indices
        are adjusted accordingly. The term is removed from position idxA's
        nonlinear node and added to position idxB's.

        :param terms: The AND terms to move.
        :type terms: list[BooleanFunction]
        :param idxA: The source bit position.
        :type idxA: int
        :param idxB: The destination bit position.
        :type idxB: int
        :raises ValueError: If a term contains non-VAR leaves, or if the
            shift would produce negative indices.
        """
        for term in terms:
            valid = True
            for var in term.args:
                if type(var) != VAR:
                    raise ValueError("types other than VAR are not currently supported")
                valid &= (var.index >= (idxA-idxB))
            if not valid:
                raise ValueError("Invalid shift attempted for term: " + str(term))
            
            newTerm = AND(*(VAR(var.index - idxA + idxB) for var in term.args))

            self.fn_list[idxA].args[-1].remove_arguments(term)
            self.fn_list[idxB].args[-1].add_arguments(newTerm)

    def getMinDestination(self, term: Any) -> int:
        """The lowest bit index to which this term can be shifted.

        :param term: An AND term.
        :type term: BooleanFunction
        :return: The minimum valid destination index.
        :rtype: int
        """
        return max((self.size - 1) - min(value.index for value in term.args), self.tau)

    def getMaxDestination(self, term: Any) -> int:
        """The highest bit index to which this term can be shifted.

        :param term: An AND term.
        :type term: BooleanFunction
        :return: The maximum valid destination index.
        :rtype: int
        """
        return min((self.size + self.tau) - (max(value.index for value in term.args)+1), self.size - 1)


    def addNonLinearTerm(self, maxAnds: int) -> None:
        """Add a random nonlinear AND term at two valid shifted positions.

        Randomly generates a product of 2 to `maxAnds` variables, then
        places it at two randomly chosen positions within the valid shift
        range. This paired placement is the core of the crossjoin
        construction.

        :param maxAnds: The maximum number of variables in the AND term.
        :type maxAnds: int
        """
        minDest = maxDest = 0

        while not (minDest < maxDest):
            numTaps = randint(2,maxAnds)
            newTerm = AND(*(VAR(i) for i in sample(range(1, self.size), numTaps)))
            maxDest = self.getMaxDestination(newTerm)
            minDest = self.getMinDestination(newTerm)

        idx1,idx2 = sample(range(minDest,maxDest+1),2)

        # add first copy
        self.fn_list[self.size - 1].args[-1].add_arguments(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx1)

        # add second copy
        self.fn_list[self.size - 1].args[-1].add_arguments(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx2)


    def generateNonlinearity(self,
        maxAnds: int = 4,
        tapDensity: float = .75
    ) -> None:
        """Populate the register with random nonlinear crossjoin terms.

        Adds paired AND terms until every bit below tau is referenced by
        at least one nonlinear term. The parameter `tapDensity` controls
        how much of the register is designated as the nonlinear region
        (tau = tapDensity * size).

        :param maxAnds: The maximum number of variables per AND term.
        :type maxAnds: int
        :param tapDensity: The fraction of bits in the nonlinear region.
        :type tapDensity: float
        """
        self.tau = min(self.tau,int(tapDensity * self.size))

        # add a set of nodes for nonlinear terms:
        for bit in range(self.size):
            self.fn_list[bit].add_arguments(XOR())

        # shift any needed linear terms to tau (always valid)
        for term in self.fn_list[self.size-1].args[0].args:
            if term.args[0].index >= self.tau: #type: ignore
                self.fn_list[self.size-1].args[1].add_arguments(term)
                self.fn_list[self.tau].args[1].add_arguments(term.shift_indices(self.tau-self.size+1))
                
        # main loop:
        tapped = set()
        while len(tapped) < self.tau:
            self.addNonLinearTerm(maxAnds)
            
            # for all nonlinear terms (tau & up)
            # determine which bits are tapped:
            for fn in self.fn_list[self.tau:]:
                for term in fn.args[-1].args:
                    tapped |= {val.index for val in term.args} #type: ignore

        # strip off any empty nodes
        for bit in range(self.size):
            nonlinear_terms = self.fn_list[bit].args[-1]
            if len(nonlinear_terms.args) == 0:
                self.fn_list[bit].remove_arguments(nonlinear_terms) 

        return

    @property
    def linear_feedback(self) -> list[Any]:
        """The linear (LFSR) portion of each bit's feedback function.

        :return: A per-bit list of the linear feedback component.
        :rtype: list[BooleanFunction]
        """
        return [f.args[0] for f in self.fn_list]

    @property
    def monomial_feedback(self) -> list[Any]:
        """The nonlinear portion of each bit's feedback function.

        For bits with no nonlinear terms, returns ``CONST(0)``.

        :return: A per-bit list of the nonlinear feedback component.
        :rtype: list[BooleanFunction]
        """
        output = []
        for f in self.fn_list:
            if len(f.args) > 1:
                output.append(XOR(*f.args[1:]))
            else:
                output.append(CONST(0))
        return output

    def compensation_list(self) -> list[BooleanFunction]:
        """Compute compensation filters for each bit.

        Each filter, when applied to the base Fibonacci LFSR's state,
        produces the same output as the corresponding crossjoin bit.
        This decomposes the crossjoin into a linear LFSR plus a set of
        nonlinear filter functions, which is useful for algebraic analysis.

        :return: A per-bit list of compensation filter functions.
        :rtype: list[BooleanFunction]
        """
        comp_list = []
        curr_fn = CONST(0)
        for fn in self.monomial_feedback[::-1]:
            # needs to be a boolean fn here to use shift_indicies
            curr_fn = curr_fn.shift_indices(-1)

            # moving in/out of BooleanANF to get nicer cancellation
            # slightly inefficient, but this is still relatively fast
            # and is much cleaner to read/write
            curr_fn = BooleanANF.from_BooleanFunction(curr_fn)
            curr_fn ^= BooleanANF.from_BooleanFunction(fn)
            curr_fn = curr_fn.to_BooleanFunction()

            comp_list.append(curr_fn.shift_indices(-1))
        
        # list is buit in reverse, so reverse when returning:
        return comp_list[::-1]
        
    def root_expressions(self) -> list[RootExpression]:
        """Compute a root expression for each bit of the register.

        Derived from the compensation filters: the degree of the largest
        AND term in each filter determines how many roots from the base
        LFSR's field appear in the exponential representation, bounding
        the linear complexity.

        :return: A per-bit list of root expressions.
        :rtype: list[RootExpression]
        """
        REs = []
        comp_list = self.compensation_list()

        for bit in range(self.size):
            term_lengths = [len(term.args) for term in comp_list[bit].args if type(term) != CONST]
            count = max(term_lengths, default = 1)
            count = min(count,self.size)
            REs.append(
                RootExpression({
                    (self.size,): set([JordanSet({self.size: count}, {1})])
                })
            )
        return REs

    def monomial_profiles(self) -> list[MonomialProfile]:
        """Compute a monomial profile for each bit of the register.

        Analogous to :meth:`root_expressions` but tracking the monomial
        structure rather than the root structure.

        :return: A per-bit list of monomial profiles.
        :rtype: list[MonomialProfile]
        """
        MPs = []
        comp_list = self.compensation_list()

        for bit in range(self.size):
            term_lengths = [len(term.args) for term in comp_list[bit].args if type(term) != CONST]
            count = max(term_lengths, default = 1)
            count = min(count,self.size)
            MPs.append(
                MonomialProfile(term_list=[TermSet({0:self.size},{0:count})])
            )
        return MPs

    @property
    def blocks(self) -> list[list[int]]:
        """The block decomposition of the register.

        A crossjoin has a single block spanning all bits.

        :return: A single-element list containing all bit indices.
        :rtype: list[list[int]]
        """
        return [list(range(self.size))]

    def filter_generator(self) -> tuple[Fibonacci, tuple[BooleanFunction,...]]:
        """Decompose this crossjoin into a base LFSR and filter functions.

        Returns the underlying Fibonacci LFSR and a tuple of per-bit filter
        functions. Applying filter[i] to the LFSR state produces the same
        output as bit i of the crossjoin.

        :return: A (base_lfsr, filters) pair.
        :rtype: tuple[Fibonacci, tuple[BoleanFunction,...]]
        """
        feedback_fn = Fibonacci(self.size, self.primitive_polynomial)
        comp_list = self.compensation_list()
        filter_fn = tuple([XOR(VAR(bit),comp_list[bit]) for bit in range(self.size)])

        return (feedback_fn,filter_fn)

    def convert_state(self,
        state: list[int] | np.ndarray[tuple[int],np.dtype[np.uint8]]
    ) -> list[int]:
        """Convert a crossjoin state to the equivalent base LFSR state.

        Applies the compensation filters to map a crossjoin state into the
        state of the underlying Fibonacci LFSR that would produce the same
        future output sequence.

        :param state: The crossjoin state vector.
        :type state: list[int]
        :return: The equivalent Fibonacci LFSR state.
        :rtype: list[int]
        """
        comp_list = self.compensation_list()
        return [
            (state[bit] ^ comp_list[bit].eval(state))
            for bit in range(self.size)
        ]