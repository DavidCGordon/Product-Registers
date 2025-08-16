from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.BooleanLogic import AND, XOR, CONST, VAR

# Linear Complexity and Monomial estimation
from PyPR.Tools.RootCounting.MonomialProfile import TermSet,MonomialProfile
from PyPR.Tools.RootCounting.JordanSet import JordanSet
from PyPR.Tools.RootCounting.RootExpression import RootExpression

# Other libs
import time

# Single bit MPR
_M1 = MPR(1,[1,1],[0,1])

class TFunction(CMPR):
    def __init__(self, size):
        super().__init__([_M1.__copy__() for i in range(size)])
        for i in range(self.size-2, -1,-1):
            self.fn_list[i].args[0].add_arguments(
                AND(*(VAR(j) for j in range(self.size-1,i,-1)))
            )

        self.fn_list[-1].add_arguments(CONST(True))
        self.induction_order = list(range(self.size-1,-1,-1))
