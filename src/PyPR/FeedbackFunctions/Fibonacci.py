from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from PyPR.BooleanLogic import VAR, BooleanANF, BooleanFunction

from PyPR.FeedbackFunctions import FeedbackFunction

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey
from PyPR.Tools.RegisterSynthesis.nlfsrSynthesis import BM_NL

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

    def _inverse_feedback(self) -> list[BooleanFunction]:
        """Build the feedback list realizing the inverse of the current map.

        Works for any feedback, linear or not, provided the register is a
        bijection. Writing the feedback function on the departing bit's
        opposite end as

            f = s[d] * A(other bits)  XOR  B(other bits)

        where `d` is the index the shift discards, the forward map is

            s'[i] = s[i +/- 1]        (the shift, which is trivially invertible)
            s'[e] = f(s)              (`e` the end the feedback writes to)

        Undoing the shift recovers every bit except s[d], and the remaining
        equation s'[e] = s[d]*A XOR B can be solved for s[d] exactly when
        A == 1 identically, giving s[d] = s'[e] XOR B. That is the
        "XOR it back out" condition: the discarded bit must enter the feedback
        linearly and on its own. When A is 0 the two states differing only in
        bit d collide, and when A is a non-constant function the map is
        non-injective on the states where A vanishes -- in both cases no
        inverse exists to build.

        The direction is read from `is_inverted`: a forward register shifts
        toward index 0 and discards s[0], so its inverse shifts toward n-1;
        an already-inverted register is the mirror image.

        :return: The feedback functions of the inverse register.
        :rtype: list[BooleanFunction]
        :raises ValueError: If the register is not a bijection, i.e. the
            discarded bit does not appear in the feedback as a lone linear term.
        """
        n = self.size

        if not self.is_inverted:
            # shift toward 0: bit 0 is discarded, the feedback lands on bit n-1.
            # Undoing it shifts toward n-1, and B's variables move down one index.
            discarded, feedback_at, offset = 0, n - 1, -1
        else:
            # mirror image: bit n-1 is discarded, the feedback lands on bit 0.
            discarded, feedback_at, offset = n - 1, 0, +1

        terms = [set(t) for t in BooleanANF.from_BooleanFunction(self.fn_list[feedback_at])]
        coefficient = [t - {discarded} for t in terms if discarded in t]
        remainder = [t for t in terms if discarded not in t]

        if coefficient != [set()]:
            raise ValueError(
                f"This register is not invertible: bit {discarded} is discarded by "
                f"the shift, so it can only be recovered if it appears in the "
                f"feedback as a lone linear term. Its coefficient in the ANF is "
                f"{'0 (the bit is absent)' if not coefficient else 'not the constant 1'}"
                f", so distinct states collide and no inverse map exists."
            )

        # The recovered bit: the feedback's own value, XOR the rest of the ANF
        # re-indexed through the undone shift. Every other bit just reads the
        # neighbour the undone shift brought it from.
        recovered: list[list[int]] = (
            [[feedback_at]] + [sorted(v + offset for v in t) for t in remainder]
        )
        new_fns: list[list[list[int]]] = [
            recovered if i == discarded else [[i + offset]]
            for i in range(n)
        ]

        return [BooleanFunction.from_ANF(fn) for fn in new_fns]

    def invert(self) -> None:
        """Toggle the register between its forward and time-reversed
        configurations.

        The forward configuration realizes the register's update map; the
        inverted configuration realizes its inverse, so that clocking the
        inverted register undoes a single clock of the original. Calling this
        twice restores the original configuration.

        This works for any invertible Fibonacci-shaped register, linear or not,
        by one of two routes.

        The work is done by :meth:`_inverse_feedback`, which solves the feedback
        for the bit the shift discards. For a linear register that is the same
        answer as the textbook route -- reverse the polynomial, then reverse the
        bit labelling -- and the equivalence is worth spelling out, because it
        is a nice illustration of why the primal/dual distinction is kept
        explicit:

        1. **Reversing the polynomial reverses time.** `Fibonacci` takes its
           input as the *dual* of the output sequence: P convolves S to zero.
           Reversal acts on both sides of that relationship at once -- P is
           dual with S exactly when the reciprocal P-tilde is dual with the
           reversed sequence S-tilde (the four-way equivalence in
           `docs/conventions/Polynomial Conventions.md`, which follows from
           reversal inverting roots on either side). So `_from_poly` applied to
           the reciprocal builds the register generating the time-reversed
           sequence.

        2. **Reversing the bit labelling realigns the window.** That is not yet
           the inverse *map*, because the two registers hold the same sliding
           window in opposite index order. A shift-down Fibonacci taking its
           output from bit 0 has bit i holding s[t+i], so its state is the
           window (s[t], ..., s[t+n-1]); the reciprocal register's state is
           that window read backwards, with bit j holding s[t+n-1-j]. The two
           agree under i -> n-1-i, which is exactly what `flip` applies. Only
           after that relabelling does the reciprocal register's transition act
           as the inverse on the *same* state vector.

        Skipping step 2 leaves a register that is not the inverse of anything:
        it clocks a valid LFSR, just one whose bits are labelled inconsistently
        with the original, so composing the two returns the starting state only
        for the handful of states fixed by the relabelling.

        `_inverse_feedback` reaches the same register without going through the
        polynomial at all, which is why it is the only route used here: it also
        covers the registers that have no polynomial to reverse, such as the
        nonlinear ones from `fromSeq(..., nonlinear=True)`. A test asserts that
        the two agree wherever both apply.

        :raises ValueError: If the register is not a bijection, so no inverse
            exists to build. Writing the feedback as f = s[d]*A XOR B for the
            discarded bit s[d], inversion needs A == 1 -- the discarded bit must
            enter the feedback linearly and on its own, so it can be XORed back
            out. A linear register always satisfies this. The registers `BM_NL`
            recovers generally do not: it optimizes for reproducing a sequence,
            not for bijectivity, and the resulting feedback carries the
            discarded bit inside larger products.
        """
        self.fn_list = self._inverse_feedback()
        self.is_inverted = not self.is_inverted

        # update_matrix is derived from fn_list, which was just rebuilt
        if 'update_matrix' in self.__dict__: del self.update_matrix

    @classmethod
    def fromSeq(cls,
        seq: list[int],
        nonlinear: bool = False,
        bijective: bool = False
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
        :param bijective: If True, the recovered register's update map is a
            bijection, so it has an inverse (`invert` works on it), preperiod 0
            from every seed, and a state graph that is a disjoint union of
            cycles. Applies to both the linear and nonlinear paths.

            Not "a well-defined period": every orbit of a map on a finite set is
            eventually periodic, so the period is always defined. What
            bijectivity adds is *pure* periodicity -- no transient.

            The costs are asymmetric. On a random sequence, or one believed to
            come from a bijective source, the flag is nearly free -- a fraction
            of a bit of register length on average -- and what it buys is mainly
            the *retention* of a property the unconstrained search would
            otherwise discard to save that bit (on uniform random inputs, a
            sampling experiment put that at roughly a third of linear fits). On a
            sequence whose transient genuinely never recurs there is
            no bijective register to find, and the search degenerates silently:
            the returned length grows with the amount of data supplied rather
            than converging. Hence the default is False -- the assertion has to
            be earned. See `docs/architecture/Bijective Register Synthesis.md`.
        :type bijective: bool
        :return: The initial state and the recovered Fibonacci register.
        :rtype: tuple[list[int], Fibonacci]
        :raises ValueError: If `bijective` is True and no bijective register
            shorter than `len(seq)` generates the sequence. Raised by the
            underlying constrained search -- `_berlekamp_massey_bijective` on
            the linear path, `_BM_NL_bijective` on the nonlinear one.
        """
        if not nonlinear:
            #run berlekamp massey to determine primitive polynomial
            size, poly = berlekamp_massey(seq, bijective = bijective)
            fn = Fibonacci(size, poly[:size+1].tolist())

        else:
            size, f = BM_NL(seq, bijective = bijective)
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
        iff bit i's update function references bit j as a `VAR` leaf. For the
        Fibonacci layout that is the downward shift together with the single
        wide XOR feeding bit n-1, so the polynomial occupies that row -- the
        mirror of :attr:`Galois.update_matrix`, where it occupies column 0.
        Any nonlinear gates in `fn_list` are ignored.

        The characteristic polynomial of M is the reverse of the stored
        polynomial. The stored polynomial is the dual (it convolves the output
        sequence to zero), whereas the characteristic polynomial of a
        transition map is the primal, whose roots are the roots of the state
        sequence; primal and dual are reverses of each other. See
        `docs/conventions/Polynomial Conventions.md`. Fibonacci and Galois
        built from the same polynomial realize the same recurrence through
        similar matrices, so they share this characteristic polynomial even
        though the matrices themselves differ.

        This is a cached property computed from `fn_list`; `invert` discards
        it, since it rebuilds that list.

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
