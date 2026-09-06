from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.BooleanLogic import AND, XOR, CONST, VAR

# Linear Complexity and Monomial estimation
from PyPR.Tools.RootCounting.MonomialProfile import TermSet, MonomialProfile
from PyPR.Tools.RootCounting.JordanSet import JordanSet
from PyPR.Tools.RootCounting.RootExpression import RootExpression
from PyPR.Tools.RootCounting.Combinatorics import powerset

# Other libs
import time

# Single bit MPR
_M1 = MPR(1,[1,1],[0,1])

class TFunction(CMPR):
    """A triangular feedback function (T-Function)

    A T-function (as defined by Klimov-Shamir) is a state-update map
    F: GF(2)^n -> GF(2)^n in which each bit's update depends only on bits at
    earlier positions in the dependency order. The triangular dependency
    structure makes T-functions a convenient substrate for proving full-period
    guarantees inductively, starting from the bit with no dependencies and
    propagating outward. The class of T-functions is broad.

    The default no-chaining construction realizes a binary counter that increments
    by 1 modulo 2^n on each clock, producing a single orbit of length 2^n. The
    counter is laid out with bit n-1 as the least significant bit — opposite
    the usual numeric convention. This choice is deliberate: it aligns with
    the high-to-low chaining flow used throughout PyPR's CMPR layouts (where
    nonlinear chaining propagates from higher index blocks into lower ones),
    so that the T-function's triangular dependency mirrors a CMPR's
    chaining structure rather than the conventional little-endian carry
    direction.

    :ivar induction_order: A traversal order over bit indices (top-down)
        suitable for inductive reasoning that exploits the triangular
        dependency structure: bit i is reached only after every bit j > i
        has already been processed.
    :vartype induction_order: list[int]
    """

    induction_order: list[int]

    def __init__(self,
        size: int
    ) -> None:
        """Construct an n-bit T-function realizing the default binary counter.

        :param size: The number of bits in the register (and the bit-width of
            the resulting binary counter).
        :type size: int
        """
        super().__init__([_M1.__copy__() for i in range(size)])
        for i in range(self.size-2, -1,-1):
            self.fn_list[i].args[0].add_arguments(
                AND(*(VAR(j) for j in range(self.size-1,i,-1)))
            )

        self.fn_list[-1].add_arguments(CONST(True))
        self.induction_order = list(range(self.size-1,-1,-1))

    def monomial_profiles(self, verbose: bool = False, force_default: bool = False) -> list[MonomialProfile]:
        mps = []
        for i in range(self.size):
            mps.append(MonomialProfile([
                TermSet({i:1 for i in index_set},{i:1 for i in index_set})
                for index_set in powerset(range(i))
            ] + [TermSet({i:1},{i:1})]))
        return mps[::-1]
        
    def root_expressions(self, locked_list: list[int] | None = None, verbose: bool = False, force_default: bool = False) -> list[RootExpression]:
        #return super().root_expressions(locked_list, verbose, force_default)
        res = []
        for i in range(self.size):
            roots = JordanSet({},set(range(1,2**i+2)))
            res.append(RootExpression({tuple():set([roots])}))
        return res[::-1]