from ProductRegisters.BooleanLogic import ANF_spec_repr, BooleanFunction, CONST
from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

class Galois(FeedbackFunction):
        
    def __init__(self, size, primitive_polynomial):
        self.size = size

        #convert koopman string into polynomial:
        if type(primitive_polynomial) == str:
            primitive_polynomial = [int(x) for x in format(int(primitive_polynomial,16), f"0>{size}b")] + [1]
        self.primitive_polynomial = primitive_polynomial

        #calculate function / taps
        self._fn_from_poly(primitive_polynomial)
        self.is_inverted = False

    #helper methods for anf construction
    def _fn_from_poly(self, polynomial):
        # build the shift (towards zero):
        newFn = [[[i+1]] for i in range(self.size-1)] + [[]]

        # build tap set fn from polynomial
        # this is derived from dubrova's shifting        
        for i in range(self.size):
            if polynomial[i+1]:
                newFn[i] += [[0]]

        self.fn_list = [BooleanFunction.construct_ANF(bitFn) for bitFn in newFn]
    
        # handle empty top case:
        if not polynomial[-1]:
            self.fn_list[-1].add_arguments(CONST(0))

    def _inverted_from_poly(self, polynomial):
        newFn = [[]] + [[[i-1]] for i in range(1,self.size)]
        for i in range(self.size):
            if polynomial[i]:
                newFn[i] += [[self.size-1]]
        self.fn_list = [BooleanFunction.construct_ANF(bitFn) for bitFn in newFn]
        
    def invert(self):
        #remake current anf based on the is_inverted attribute
        if not self.is_inverted:
            self._inverted_from_poly(self.primitive_polynomial)
        else:
            self._fn_from_poly(self.primitive_polynomial)

    @classmethod
    def fromSeq(self,seq):
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
        return s, Galois(L, c[:L+1])

    @classmethod 
    def fromReg(self, F, bit):
        numIters = 2*F.size + 4
        seq = [state[bit] for state in F.run(numIters)]
        return Galois.fromSeq(seq)

