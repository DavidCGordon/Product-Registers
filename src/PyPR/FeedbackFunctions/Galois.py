from typing import Self, TYPE_CHECKING

from PyPR.BooleanLogic import BooleanANF, BooleanFunction, CONST
from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

if TYPE_CHECKING:
    from PyPR.FeedbackRegister import FeedbackRegister

class Galois(FeedbackFunction):
    """A linear feedback shift register in Galois configuration.

    The state shifts down (values progress toward bits with lower indices) and
    the bit departing position 0 is XORed back into every position whose
    corresponding polynomial coefficient is nonzero, distributing the feedback
    across several positions for a shorter critical path than the equivalent
    Fibonacci layout. When the register is built from a primitive polynomial
    of degree n the nonzero states form a single cycle of length 2^n - 1.

    The input polynomial follows the dual interpretation: it is the
    polynomial that convolves the output sequence to zero (matching the
    convention returned by Berlekamp-Massey), and the output sequence at any
    single bit satisfies its reverse. Equivalently, the register state is
    identified with an element of GF(2)[x] / P and the downward shift acts
    as multiplication by x^{-1}, which is what makes the input polynomial
    the dual rather than the primal. See
    `docs/conventions/Polynomial Conventions.md` for the full primal/dual
    framework. As a practical consequence,
    `Galois(*Galois.fromSeq(seq))` reproduces `seq` with no manual
    reversal at the user boundary.

    :ivar size: The number of bits in the register, equal to the degree of
        the polynomial.
    :vartype size: int
    :ivar primitive_polynomial: The coefficient list of the polynomial, of
        length n+1, with `primitive_polynomial[i]` the coefficient of x^i.
    :vartype primitive_polynomial: list[int]
    :ivar is_inverted: True if the register is currently configured to clock
        in the time-reversed direction; False for the forward direction.
    :vartype is_inverted: bool
    """

    primitive_polynomial: list[int]
    is_inverted: bool

    def __init__(self,
        size: int,
        primitive_polynomial: str | list[int]
    ) -> None:
        """Construct a Galois LFSR of the given size with the specified
        primitive polynomial.

        The polynomial is interpreted under the dual convention (see the
        class docstring): the output sequence at any single bit satisfies
        its reverse, not the input itself.

        :param size: The number of bits in the register, equal to the degree
            of the polynomial.
        :type size: int
        :param primitive_polynomial: The primitive polynomial of degree n.
            May be given either as a coefficient list of length n+1 (with
            index i holding the coefficient of x^i) or as a Koopman hex
            string. The Koopman format encodes the polynomial's lower n
            coefficients in hex; the leading 1 at degree n is implicit.
        :type primitive_polynomial: str | list[int]
        """
        self.size = size

        #convert koopman string into polynomial:
        if isinstance(primitive_polynomial, str):
            primitive_polynomial = [int(x) for x in format(int(primitive_polynomial,16), f"0>{size}b")] + [1]
        self.primitive_polynomial = primitive_polynomial

        #calculate function / taps
        self._fn_from_poly(primitive_polynomial)
        self.is_inverted = False

    #helper methods for anf construction
    def _fn_from_poly(self,
        polynomial: list[int]
    ) -> None:
        """Populate `fn_list` with the forward-time Galois feedback functions
        for the given polynomial.

        Builds the shift toward index 0 (each bit i takes the value of bit
        i+1 in the previous state) and adds the bit-0 contribution to every
        position whose polynomial coefficient is nonzero, realizing the
        Galois tap pattern derived from Dubrova's shift representation.

        :param polynomial: The coefficient list of the polynomial whose
            taps should be installed.
        :type polynomial: list[int]
        """
        # build the shift (towards zero):
        newFn = [[[i+1]] for i in range(self.size-1)] + [[]]

        # build tap set fn from polynomial
        # the shift by 1
        for i in range(self.size):
            if polynomial[i+1]:
                newFn[i] += [[0]]

        self.fn_list = [BooleanFunction.from_ANF(bitFn) for bitFn in newFn]

        # handle empty top case:
        if not polynomial[-1]:
            self.fn_list[-1].add_arguments(CONST(0))

    def _inverted_from_poly(self,
        polynomial: list[int]
    ) -> None:
        """Populate `fn_list` with the time-reversed Galois feedback functions
        for the given polynomial.

        Builds the shift toward index n-1 (each bit i takes the value of
        bit i-1 in the previous state) with taps derived so that the
        resulting register clocks the inverse of the forward update map.

        :param polynomial: The coefficient list of the polynomial whose
            inverse-time taps should be installed.
        :type polynomial: list[int]
        """
        newFn = [[]] + [[[i-1]] for i in range(1,self.size)]
        for i in range(self.size):
            if polynomial[i]:
                newFn[i] += [[self.size-1]]
        self.fn_list = [BooleanFunction.from_ANF(bitFn) for bitFn in newFn]

    def invert(self) -> None:
        """Toggle the register between its forward and time-reversed
        configurations.

        The forward configuration realizes the linear map associated with P;
        the inverted configuration realizes the inverse map, so that
        clocking the inverted register undoes a single clock of the original.
        Both configurations traverse the same orbit but in opposite directions.
        """
        #remake current anf based on the is_inverted attribute
        if not self.is_inverted:
            self._inverted_from_poly(self.primitive_polynomial)
        else:
            self._fn_from_poly(self.primitive_polynomial)

    @classmethod
    def fromSeq(cls,
        seq: list[int]
    ) -> tuple[list[int], "Galois"]:
        """Recover the Galois LFSR which generates a given binary sequence.

        Applies Berlekamp-Massey to recover the minimal characteristic
        polynomial of the sequence, then back-solves for the initial state
        of the Galois realization. The returned `(seed, register)` pair
        satisfies: clocking `register` from `seed` reproduces `seq` in the 
        values of bit 0 over time.

        :param seq: A prefix of a binary sequence, of length at least 2L,
            where L is the linear complexity of the sequence.
        :type seq: list[int]
        :return: The initial state and the recovered Galois LFSR.
        :rtype: tuple[list[int], Galois]
        """
        #run berlekamp massey to determine primitive polynomial
        L, c = berlekamp_massey(seq)

        # calculate inital state:
        s = []
        for i in range(L):
            s_i = 0
            for j in range(i+1):
                s_i ^= (seq[i-j] & c[j])
            s.append(s_i)

        #return Galois LFSR parameters
        return s, Galois(L, c[:L+1].tolist())

    @classmethod
    def fromReg(cls,
        F: "FeedbackRegister",
        bit: int = 0,
        numIters: int | None = None
    ) -> tuple[list[int], "Galois"]:
        """Recover the Galois LFSR which reproduces a single bit of a running
        register.

        Observes the specified bit of `F` over `numIters` clock cycles to
        obtain a binary sequence, then applies `fromSeq` to recover the
        Galois LFSR realizing it. This is useful for converting an arbitrary
        (possibly nonlinear) register into a linear surrogate that matches
        a chosen output bit.

        :param F: The source register to observe.
        :type F: FeedbackRegister
        :param bit: The index of the bit to observe.
        :type bit: int
        :param numIters: The number of clock cycles to observe. Defaults to
            2*F.size + 4, which is enough for Berlekamp-Massey to recover
            any linear recurrence of complexity at most F.size.
        :type numIters: int | None
        :return: The initial state and the recovered Galois LFSR.
        :rtype: tuple[list[int], Galois]
        """
        if not numIters:
            numIters = 2*F.size + 4
        seq = [state[bit] for state in F.run(numIters)]
        return Galois.fromSeq(seq)

