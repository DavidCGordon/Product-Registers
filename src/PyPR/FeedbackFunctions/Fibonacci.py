from typing import TYPE_CHECKING

from PyPR.BooleanLogic import BooleanANF, BooleanFunction, VAR
from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey
from PyPR.Tools.RegisterSynthesis.nlfsrSynthesis import BM_NL

from functools import cached_property

import numpy as np

if TYPE_CHECKING:
    from PyPR.FeedbackRegister import FeedbackRegister

class Fibonacci(FeedbackFunction):
    """A linear (or nonlinear) feedback shift register in Fibonacci configuration.

    The state shifts down (values progress toward bits with lower indices)
    and the bit at position n-1 is recomputed each cycle as the XOR of every
    position whose corresponding polynomial coefficient is nonzero. When the
    register is built from a primitive polynomial of degree n, the nonzero
    states form a single cycle of length 2^n - 1.

    The input polynomial follows the dual interpretation: it is the
    polynomial that convolves the output sequence to zero (matching the
    convention returned by Berlekamp-Massey), and the output sequence at any
    single bit satisfies its reverse. This is the disambiguation case where
    primal/dual terminology pays off: see
    `docs/conventions/Polynomial Conventions.md` for the full primal/dual framework
    and worked-out matrix examples. As a practical consequence,
    `Fibonacci(*Fibonacci.fromSeq(seq))` reproduces `seq` with no manual
    reversal at the user boundary.

    Fibonacci is mathematically equivalent to Galois with the same input
    polynomial — both realize the same linear recurrence — but Fibonacci
    concentrates the XOR into a single wide gate at position n-1, giving a
    longer critical path but a tap pattern that mirrors the polynomial
    coefficients more directly. The class also supports a nonlinear
    Fibonacci-style register via `fromSeq(..., nonlinear=True)`, in which
    the top-bit feedback is an arbitrary boolean function recovered by
    nonlinear Berlekamp-Massey.

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

    def _from_poly(self,
        size: int,
        primitive_polynomial: list[int]
    ) -> list[BooleanFunction]:
        """Build the Fibonacci feedback functions for a given polynomial.

        The result is a `size`-long list of boolean functions: for each
        i < n-1, position i takes the value of position i+1 in the previous
        state (the downward shift); position n-1 receives the XOR of every
        previous-state position whose corresponding polynomial coefficient
        is nonzero (excluding the leading-degree coefficient, which is the
        identity contribution and would tap the bit being computed).

        :param size: The number of bits in the register.
        :type size: int
        :param primitive_polynomial: The polynomial as a coefficient list of
            length n+1.
        :type primitive_polynomial: list[int]
        :return: The list of bit-update boolean functions, one per position.
        :rtype: list[BooleanFunction]
        """
        # Create the top function (reciprocal polynomial):
        # Example: [1,0,1,1] -> [0,2,3] -> [[3],[1],[0]]
        topFn = [[size-idx] for (idx, t) in enumerate(primitive_polynomial) if t == 1]
        topFn = [topFn[1:]]

        # Cosmetic Changes to order:
        topFn = topFn[::-1]

        #create the shift for all other bits
        shiftFn = [[[i+1]] for i in range(size-1)]

        return [BooleanFunction.from_ANF(bitFn) for bitFn in (shiftFn + topFn)]

    def __init__(self,
        size: int,
        primitive_polynomial: str | list[int]
    ) -> None:
        """Construct a Fibonacci LFSR of the given size with the specified
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
        #convert koopman strings
        if isinstance(primitive_polynomial, str):
            primitive_polynomial = [int(x) for x in format(int(primitive_polynomial,16), f"0>{size}b")] + [1]

        #assign attributes
        self.primitive_polynomial = primitive_polynomial
        self.size = size
        self.is_inverted = False

        self.fn_list = self._from_poly(size,primitive_polynomial)

    def invert(self) -> None:
        """Toggle the register between its forward and time-reversed
        configurations.

        The forward configuration realizes the linear map associated with the
        stored polynomial; the inverted configuration realizes its inverse,
        so that clocking the inverted register undoes a single clock of the
        original. Both configurations traverse the same orbit but in opposite
        directions. Internally, time-reversal swaps the stored polynomial for
        its reciprocal and re-runs the Fibonacci construction.

        .. note::
            This method currently has a known issue: the bit-relabeling step
            (`self.flip()`) needed to align indices after polynomial reversal
            is unimplemented (see the TODO in the source). Until that is
            fixed, calling `invert()` will produce a register whose tap
            pattern matches the inverse polynomial but whose bit indices are
            not aligned to the original layout.
        """
        if not self.is_inverted:
            #reverse taps:
            self.fn_list = self._from_poly(self.size,self.primitive_polynomial[::-1])
            # flip bit labelling:
            # TODO: FIX THIS LINE-> self.flip()
            self.is_inverted = True
        else:
            #remake anf:
            self.fn_list = self._from_poly(self.size,self.primitive_polynomial)
            self.is_inverted = False

    @classmethod
    def fromSeq(cls,
        seq: list[int],
        nonlinear: bool = False
    ) -> tuple[list[int], "Fibonacci"]:
        """Recover the Fibonacci register which generates a given binary sequence.

        In the linear case, applies Berlekamp-Massey to recover the dual
        polynomial annihilating the sequence and constructs a Fibonacci LFSR
        from it directly (no reversal, since the dual is exactly what
        `Fibonacci` expects). In the nonlinear case, applies the nonlinear
        Berlekamp-Massey variant to recover the smallest boolean function
        satisfying the sequence as a feedback rule and installs it as the
        top-bit update of a fresh Fibonacci register.

        The returned `(seed, register)` pair satisfies: clocking `register`
        from `seed` reproduces `seq` on bit 0.

        :param seq: A prefix of a binary sequence. For the linear case, must
            be of length at least 2L where L is the linear complexity; for
            the nonlinear case, length requirements depend on the recovered
            function's complexity.
        :type seq: list[int]
        :param nonlinear: If True, recover a nonlinear feedback rule using
            `BM_NL`; otherwise recover a linear one using Berlekamp-Massey.
        :type nonlinear: bool
        :return: The initial state and the recovered Fibonacci register.
        :rtype: tuple[list[int], Fibonacci]
        """
        if not nonlinear:
            #run berlekamp massey to determine primitive polynomial
            size, poly = berlekamp_massey(seq)
            fn = Fibonacci(size, poly[:size+1].tolist())

        else:
            size, f = BM_NL(seq)
            fn = Fibonacci(size, [])
            fn[size-1].add_arguments(f)

        #return Fibonacci LFSR parameters
        init_state = seq[:size]
        return init_state, fn

    @classmethod
    def fromReg(cls,
        F: "FeedbackRegister",
        bit: int = 0,
        numIters: int | None = None,
        nonlinear: bool = False
    ) -> tuple[list[int], "Fibonacci"]:
        """Recover the Fibonacci register which reproduces a single bit of a
        running register.

        Observes the specified bit of `F` over `numIters` clock cycles to
        obtain a binary sequence, then applies `fromSeq` to recover the
        Fibonacci register realizing it. Useful for converting an arbitrary
        register into a (possibly nonlinear) Fibonacci surrogate that matches
        a chosen output bit.

        :param F: The source register to observe.
        :type F: FeedbackRegister
        :param bit: The index of the bit to observe.
        :type bit: int
        :param numIters: The number of clock cycles to observe. Defaults to
            2*F.size + 4, which is enough for linear Berlekamp-Massey to
            recover any linear recurrence of complexity at most F.size.
        :type numIters: int | None
        :param nonlinear: Forwarded to `fromSeq`; if True, attempt to recover
            a nonlinear feedback rule.
        :type nonlinear: bool
        :return: The initial state and the recovered Fibonacci register.
        :rtype: tuple[list[int], Fibonacci]
        """
        if not numIters:
            numIters = 2*F.size + 4

        seq = [state[bit] for state in F.run(numIters)]
        return Fibonacci.fromSeq(seq, nonlinear = nonlinear)

    @cached_property
    def update_matrix(self) -> np.ndarray:
        """The GF(2) update matrix of the linear part of the register.

        The matrix M satisfies `next_state = M @ current_state` (over GF(2))
        on the linear portion of the feedback functions: entry `M[i, j]` is 1
        iff bit i's update function references bit j as a `VAR` leaf. Any
        nonlinear gates in `fn_list` are ignored.

        .. todo::
            Verify the precise relationship between the characteristic
            polynomial of M and the stored primitive polynomial for this
            register's specific layout. The expected relation (under the
            dual interpretation of the input) is that the characteristic
            polynomial of M equals the reverse of the stored polynomial,
            but this should be checked against the shift-down construction
            before being committed to the docstring.

        :return: The size x size GF(2) update matrix.
        :rtype: numpy.ndarray
        """
        matrix = np.zeros([self.size,self.size], dtype = 'uint8')

        for inpt in range(self.size):
            for outpt in range(self.size):
                #iterate through the VAR objects in the linear function portion
                for leaf in self.fn_list[outpt].inputs():
                    if isinstance(leaf,VAR) and leaf.index == inpt:
                        matrix[outpt][inpt] = 1

        return matrix