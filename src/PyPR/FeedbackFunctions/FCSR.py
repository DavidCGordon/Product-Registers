from typing import TYPE_CHECKING

from PyPR.BooleanLogic import BooleanFunction, XOR, AND, VAR, CONST
from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.Tools.RegisterSynthesis.fcsrSynthesis import BM_FCSR

from math import ceil, floor, log2

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

    :ivar size: The total number of bits in the register (values + carries),
        equal to 2 * diadic_complexity - 1.
    :vartype size: int
    :ivar connection_int: The connection integer q (an odd positive integer),
        which plays the role that the primitive polynomial plays in an LFSR.
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
        realizes that fraction. The fraction is intentionally not simplified,
        so that the resulting state is valid for a larger FCSR if desired.

        For negative numerators, each (value, carry) pair contributes
        value + 2 * carry to the total, so the natural partition of
        |num| distributes by thirds: the quotient |num| // 3 is split
        between values and carries, with the remainder allocated to
        whichever component absorbs the extra unit.

        :param num: The numerator p of the 2-adic fraction.
        :type num: int
        :param den: The denominator q (connection integer, positive).
        :type den: int
        :return: The 2-adic complexity and the interleaved initial state.
        :rtype: tuple[int, list[int]]
        """
        # fraction is not simplified in order to
        # create valid states for larger FCSRs
        # but this does assume no negative denominators

        # handle 0/1 and 1/1 edge case (undefined log)
        if den == 1 and num in (0,1):
            return (1,[num])

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

