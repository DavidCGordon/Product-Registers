from PyPR.BooleanLogic import BooleanFunction, XOR, AND, VAR, CONST
from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.Tools.RegisterSynthesis.fcsrSynthesis import BM_FCSR

from math import ceil, floor, log2

# Chosen layout:
# carry feeding into the adder ([2*i]) = [2*i + 1]
# value feeding into the adder ([2*i]) = [2*(i+1)]

class FCSR(FeedbackFunction):
    def __init__(self, diadic_complexity, q):
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
    def carries(self):
        return [self.fn_list[2*i + 1] for i in range(self.size//2)]

    @property
    def values(self):
        return [self.fn_list[2*i] for i in range(self.size//2)]

    @classmethod
    def fromSeq(cls,seq):
        #run berlekamp massey to determine primitive polynomial
        size, num, den = BM_FCSR(seq)
        size, init_state = FCSR.state_from_frac(num,den)
        return init_state, FCSR(size, den)

    @classmethod 
    def fromReg(cls, F, bit = 0, numIters = None):
        if not numIters:
            numIters = 2*F.size + 4
            
        seq = [state[bit] for state in F.run(numIters)]
        return FCSR.fromSeq(seq)

    @classmethod
    def state_from_frac(cls,num,den):
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

