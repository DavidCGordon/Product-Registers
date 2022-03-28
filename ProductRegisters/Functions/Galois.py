from BitVector import BitVector

from ProductRegisters import ANF
from ProductRegisters.Functions import FeedbackFunction
from ProductRegisters.Tools.BerlekampMassey import berlekamp_massey

class Galois(FeedbackFunction):
        
    def __init__(self, size, primitive_polynomial):
        self.size = size

        #convert koopman string into polynomial:
        if type(primitive_polynomial) == str:
            primitive_polynomial = list(BitVector( \
                intVal = int(primitive_polynomial, 16), size = size)) + [1]
        self.primitive_polynomial = primitive_polynomial

        #calculate function / taps
        self._anf_from_poly(primitive_polynomial)
        self.is_inverted = False

    #helper methods for anf construction
    def _anf_from_poly(self, polynomial):
        #build tap set fn from polynomial
        newFn = [[[i+1]] for i in range(self.size-1)] + [[]]
        for i in range(self.size):
            if polynomial[i+1]:
                newFn[i] += [[0]]
        self.anf = [ANF(bitFn) for bitFn in newFn]
    
    def _inverted_from_poly(self, polynomial):
        newFn = [[]] + [[[i-1]] for i in range(1,self.size)]
        for i in range(self.size):
            if polynomial[i]:
                newFn[i] += [[self.size-1]]
        self.anf = [ANF(bitFn) for bitFn in newFn]
        
    def invert(self):
        #remake current anf based on the is_inverted attribute
        if not self.is_inverted:
            self._inverted_from_poly(self.primitive_polynomial)
        else:
            self._anf_from_poly(self.primitive_polynomial)

    @classmethod
    def fromSeq(self,seq):
        #run berlekamp massey to determine primitive polynomial
        L, c = berlekamp_massey(seq)

        #calculate inital state:
        s = []
        for i in range(L):
            s_i = 0
            for j in range(i+1):
                s_i ^= (seq[i-j] & c[j])
            s.append(s_i)
        s = BitVector(bitlist = s[::-1])

        #return Galois LFSR parameters 
        return s, Galois(L, c[:L+1])

    @classmethod 
    def fromReg(self, F, bit):
        numIters = 2*F.size + 4
        seq = [state[bit] for state in F.run(numIters)]
        return Galois.fromSeq(seq)

