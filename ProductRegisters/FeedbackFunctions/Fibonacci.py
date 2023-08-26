from ProductRegisters.BooleanLogic import ANF_spec_repr, BooleanFunction
from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

class Fibonacci(FeedbackFunction):

    def _from_poly(self,size,primitive_polynomial):
        #Create the top function:
        # Example: [1,0,1,1] -> [0,2,3] -> [[3],[1],[0]]
        topFn = [[size-idx] for (idx, t) in enumerate(primitive_polynomial) if t == 1]
        # Example: [[3],[1],[0]] -> [[0],[1],[3]]
        topFn = topFn[::-1]
        # Example: [[0],[1],[3]] -> [[[0],[1]]]
        topFn = [topFn[:-1]]

        #create the shift for all other bits
        shiftFn = [[[i+1]] for i in range(size-1)]

        return [BooleanFunction.construct_ANF(bitFn) for bitFn in (shiftFn + topFn)]

    def __init__(self, size, primitive_polynomial):
        #convert koopman strings
        if type(primitive_polynomial) == str:
            primitive_polynomial = [int(x) for x in format(int(primitive_polynomial,16), f"0>{size}b")] + [1]
        
        #assign attributes
        self.primitive_polynomial = primitive_polynomial
        self.size = size        
        self.is_inverted = False

        self.fn_list = self._from_poly(size,primitive_polynomial)

    def invert(self):
        if not self.is_inverted:
            #reverse taps:
            self.fn_list = self._from_poly(self.size,self.primitive_polynomial[::-1])
            #flip bit labelling:
            # TODO: FIX THIS LINE-> self.flip()
            self.is_inverted = True
        else:
            #remake anf:
            self.fn_list = self._from_poly(self.size,self.primitive_polynomial)
            self.is_inverted = False
            
    @classmethod
    def fromSeq(self,seq):
        #run berlekamp massey to determine primitive polynomial
        L, c = berlekamp_massey(seq)
        
        #return Fibonacci LFSR parameters
        init_state = seq[:L]
        return init_state, Fibonacci(L, c[:L+1])

    @classmethod 
    def fromReg(self, F, bit):
        numIters = 2*F.size + 4
        seq = [state[bit] for state in F.run(numIters)]
        return Fibonacci.fromSeq(seq)