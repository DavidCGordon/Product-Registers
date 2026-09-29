from math import ceil, floor, log2
from typing import TYPE_CHECKING

from PyPR.BooleanLogic import AND, VAR, XOR, BooleanFunction

from PyPR.FeedbackFunctions import FeedbackFunction

from PyPR.Tools.RegisterSynthesis.fcsrSynthesis import BM_FCSR

if TYPE_CHECKING:
    from PyPR.FeedbackRegister import FeedbackRegister

# Chosen layout:
# carry feeding into the adder ([2*i]) = [2*i + 1]
# value feeding into the adder ([2*i]) = [2*(i+1)]

class FCSR(FeedbackFunction):
    """A Feedback with Carry Shift Register (FCSR).

    An FCSR is the 2-adic analogue of a linear feedback shift register:
    where an LFSR realizes a linear recurrence over GF(2) (equivalently,
    a rational element of the formal power series ring GF(2)[[x]]), an
    FCSR realizes a rational element of the 2-adic integers. The output
    sequence of an FCSR is eventually periodic and determined by the
    2-adic fraction p/q, where q is the connection integer and p is fixed
    by the initial state.

    The register interleaves value bits and carry bits: position 2i holds
    a value bit and position 2i+1 holds the carry that feeds into the
    adder at position 2i. The total register width is 2d - 1, where d is
    the 2-adic complexity (there are d value positions and d - 1 carry
    positions).

    The state encodes the numerator, and how it is read depends on whether
    the feedback reaches the top value cell. Write a for the value word and
    c for the carry word. When q >= 2^d - 1 the loop is closed, every weight
    is positive, and a register in that state emits the expansion of
    -(a + 2c) / q, so the numerator is non-positive. When q < 2^d - 1 the
    loop is open: the top value cell holds its own value, which forces its
    weight to be -2^(d-1) rather than +2^(d-1). The reading is then
    a + 2c - 2^d * v_{d-1}, which is the value word taken as a
    two's-complement integer at width d, plus 2c -- the carry word is not
    part of that pattern -- so the numerator may have either sign.
    FCSR(2, 1) is open, and its state [0, 0, 1] has a = 2, c = 0, giving a
    reading of 2 - 4 = -2 and hence numerator +2. See
    `docs/architecture/FCSR Implementation.md` for the layout, the tap
    derivation and the encoding, and
    `docs/theory/2-adic Integers and Rational Sequences.md` for the arithmetic.

    :ivar size: The total number of bits in the register (values + carries),
        equal to 2 * diadic_complexity - 1.
    :vartype size: int
    :ivar connection_int: The connection integer q (an odd positive integer),
        which plays the role that the connection polynomial plays in an LFSR:
        the tap positions are the set bits of (q + 1) / 2. Primitivity is an
        extra property q may or may not have; the constructor does not
        require it.
    :vartype connection_int: int
    """

    connection_int: int

    def __init__(self,
        diadic_complexity: int,
        q: int
    ) -> None:
        """Construct an FCSR with the given 2-adic complexity and connection
        integer.

        The connection integer q encodes the feedback taps via its binary
        representation: the tap positions are derived from (q + 1) / 2. At
        each tapped position, a full binary adder replaces the simple XOR
        of the base shift, introducing carry propagation.

        :param diadic_complexity: The 2-adic complexity of the register,
            equal to the number of value-bit positions. The total register
            width is 2 * diadic_complexity - 1.
        :type diadic_complexity: int
        :param q: The connection integer (an odd positive integer).
        :type q: int
        """
        taps = [int(x) for x in bin((q+1)//2)[2:][::-1]]
        #print(taps,len(taps))

        self.connection_int = q
        self.size = 2*diadic_complexity - 1
        # self.size // 2 gives number of non-carry bits

        self.fn_list = [BooleanFunction() for _ in range(self.size)]

        # place base connections:
        for i in range(diadic_complexity-1):
            self.fn_list[2*i] = XOR(VAR(2*(i+1)), VAR(2*i+1))
            self.fn_list[2*i + 1] = AND(VAR(2*(i+1)), VAR(2*i+1))
        self.fn_list[-1] = VAR(self.size-1)

        # overwrite functions for bits with feed in:
        for i, tap in enumerate(taps):
            if tap:
                #print(i)
                if i == self.size // 2:
                    self.fn_list[2*i] = VAR(0)
                else:
                    self.fn_list[2*i] = XOR(VAR(2*(i+1)), VAR(2*i + 1), VAR(0))
                    self.fn_list[2*i + 1] = XOR(
                        AND(VAR(2*(i+1)), VAR(2*i + 1)),
                        AND(VAR(2*(i+1)), VAR(0)),
                        AND(VAR(2*i + 1), VAR(0))
                    )

        #self.fn_list[-1] = VAR(0)

    @property
    def carries(self) -> list[BooleanFunction]:
        """The carry-bit feedback functions (one per adder position).

        :return: The d - 1 carry functions at odd-indexed positions.
        :rtype: list[BooleanFunction]
        """
        return [self.fn_list[2*i + 1] for i in range(self.size//2)]

    @property
    def values(self) -> list[BooleanFunction]:
        """The value-bit feedback functions for the paired positions.

        Returns the d - 1 value functions at even-indexed positions that
        have a corresponding carry bit. The final value position
        (index 2(d-1)) is not included.

        :return: The paired value functions at even-indexed positions.
        :rtype: list[BooleanFunction]
        """
        return [self.fn_list[2*i] for i in range(self.size//2)]

    @classmethod
    def fromSeq(cls,
        seq: list[int]
    ) -> tuple[list[int], "FCSR"]:
        """Recover the FCSR which generates a given binary sequence.

        Applies the 2-adic Berlekamp-Massey algorithm to recover the
        rational fraction p/q whose 2-adic expansion matches the sequence,
        then constructs an FCSR from the connection integer q and computes
        the initial state from p/q.

        :param seq: A prefix of a binary sequence long enough for the
            2-adic Berlekamp-Massey algorithm to converge.
        :type seq: list[int]
        :return: The initial state and the recovered FCSR.
        :rtype: tuple[list[int], FCSR]
        """
        #run berlekamp massey to determine primitive polynomial
        size, num, den = BM_FCSR(seq)
        size, init_state = FCSR.state_from_frac(num,den)
        return init_state, FCSR(size, den)

    @classmethod
    def fromReg(cls,
        F: "FeedbackRegister",
        bit: int = 0,
        numIters: int | None = None
    ) -> tuple[list[int], "FCSR"]:
        """Recover the FCSR which reproduces a single bit of a running
        register.

        Observes the specified bit of `F` over `numIters` clock cycles to
        obtain a binary sequence, then applies `fromSeq` to recover the
        FCSR realizing it.

        :param F: The source register to observe.
        :type F: FeedbackRegister
        :param bit: The index of the bit to observe.
        :type bit: int
        :param numIters: The number of clock cycles to observe. Defaults to
            2*F.size + 4.
        :type numIters: int | None
        :return: The initial state and the recovered FCSR.
        :rtype: tuple[list[int], FCSR]
        """
        if not numIters:
            numIters = 2*F.size + 4

        seq = [state[bit] for state in F.run(numIters)]
        return FCSR.fromSeq(seq)

    @classmethod
    def state_from_frac(cls,
        num: int,
        den: int
    ) -> tuple[int, list[int]]:
        """Convert a 2-adic rational fraction to an FCSR initial state.

        Given a fraction p/q (where q is the connection integer), computes
        the register size and the interleaved value/carry state vector that
        realizes that fraction. The fraction is intentionally not simplified:
        a fraction and its reduction generate the same sequence but name
        different registers, since q is the connection integer and reducing
        it changes the taps. Keeping p/q as given is what lets `fromSeq`
        hand the same den to this method and to the constructor, so that
        state and register are built for one q.

        For non-positive numerators, each (value, carry) pair contributes
        value + 2 * carry, and |num| is split by thirds: the quotient
        |num| // 3 goes to each of the value word a and the carry word c,
        with the remainder allocated to a (remainder 1) or c (remainder 2).
        The split is a choice among those with a + 2c = |num|, and the one
        made here is the one guaranteeing a < 2^(size-1) and c < 2^(size-1):
        the value word leaves the top cell -- the open loop's sign bit --
        clear, so the state reads the same in either regime, and the carry
        word fits its size - 1 cells. See
        `docs/architecture/FCSR Implementation.md` §Encoding a fraction for
        the residue argument and for a split that fails.

        A positive numerator takes the other route. Its expansion terminates,
        which a register can only produce with its feedback left open, so the
        size is chosen to put den below 2^size - 1 and the value word is set to
        2^size - num -- the two's complement pattern of -num at that width, its
        top cell being the sign bit the open loop reads as -2^(size-1). The
        carry word is empty. 0/1 and 1/1 are special-cased: the first needs no
        register, and the second is the one fraction where the logarithmic
        sizing would close the loop.

        :param num: The numerator p of the 2-adic fraction.
        :type num: int
        :param den: The denominator q (connection integer, positive).
        :type den: int
        :return: The 2-adic complexity and the interleaved initial state.
        :rtype: tuple[int, list[int]]
        """
        # fraction is not simplified: reducing it changes the denominator,
        # hence the taps, hence which register the state belongs to
        # but this does assume no negative denominators

        # handle 0/1 edge case (undefined log)
        if den == 1 and num == 0:
            return (1,[0])

        # to handle the 1/1 edge case, we need an extra bit to
        # turn off the feedback (to get the all zeros state)
        if den == 1 and num == 1:
            return (2,[1,0,1])

        if  num > 0:
            size = 1 + ceil(log2(max(den,abs(num))))
            values = 2**size - num
            carries = 0

        else:
            size = max(
                ceil(log2(abs(num)/3 + 1)) + 1,
                ceil(log2(den))
            )

            values = carries = floor(abs(num) / 3)
            if floor(abs(num)) % 3 == 1:
                values += 1
            elif floor(abs(num)) % 3 == 2:
                carries += 1
            else:
                pass

        # convert a,c into a state.
        out_state = [0 for i in range(2*size-1)]
        for i,bit in enumerate([int(x) for x in bin(values)[2:][::-1]]):
            out_state[2*i] = bit
        for i,bit in enumerate([int(x) for x in bin(carries)[2:][::-1]]):
            out_state[2*i+1] = bit

        return size, out_state

