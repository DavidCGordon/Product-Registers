from ProductRegisters.BooleanLogic import XOR, AND, VAR
from ProductRegisters.FeedbackFunctions import FeedbackFunction


class FCSR(FeedbackFunction):
    def __init__(self,q):
        taps = [int(x) for x in bin((q+1)//2)[2:][::-1]]

        self.connection_int = q

        self.size = len(taps)*2 - 1
        # self.size // 2 gives number of non-carry bits

        self.fn_list = [None for _ in range(self.size)]
        # carry feeding into the adder ([2*i]) = [2*i + 1]
        # value feeding into the adder ([2*i]) = [2*(i+1)]

        for i, tap in enumerate(taps[:-1]):
            if tap:
                self.fn_list[2*i] = XOR(VAR(2*(i+1)), VAR(2*i + 1), VAR(0))
                self.fn_list[2*i + 1] = XOR(
                    AND(VAR(2*(i+1)), VAR(2*i + 1)),
                    AND(VAR(2*(i+1)), VAR(0)),
                    AND(VAR(2*i + 1), VAR(0))
                )
            else:
                self.fn_list[2*i] = XOR(VAR(2*(i+1)), VAR(2*i+1))
                self.fn_list[2*i + 1] = AND(VAR(2*(i+1)), VAR(2*i+1))
        self.fn_list[-1] = VAR(0)

    @property
    def carries(self):
        return [self.fn_list[2*i + 1] for i in range(self.size//2)]

    @property
    def values(self):
        return [self.fn_list[2*i] for i in range(self.size//2)]

    @classmethod
    def state_from_numerator(self,size,n):
        out_state = [0 for i in range(size)]

        remainder = n % 2
        a = (n-1) // 2
        c = (n-a) // 2

        a = [int(x) for x in bin(a)[2:][::-1]]
        c = [int(x) for x in bin(c)[2:][::-1]]
        c.append(0) # match lengths, as a < 2c

        for i in range(len(a)):
            out_state[2*i]     = a[i]
            out_state[2*i + 1] = c[i]
        out_state[0] = remainder

        return out_state

