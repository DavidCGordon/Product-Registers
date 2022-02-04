from BitVector import BitVector

from ProductRegisters.Functions import FeedbackFunction
from ProductRegisters import ANF

from copy import deepcopy

#inversion = opposite polynomial + bit flip
#alternatively
#reciprocal polynomial = inversion + bit flip

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
        #remove current taps:
        if not self.is_inverted:
            self._inverted_from_poly(self.primitive_polynomial)
        else:
            self._anf_from_poly(self.primitive_polynomial)


    #returns an equivalent Galois LFSR ANF and initial state:
    @classmethod
    def fromSeq(self,seq):
        N = len(seq)
        #c = current connection polynomial
        c = [1] + [0 for i in range(N-1)]
        #b = prev. connection polynomial
        b = [1] + [0 for i in range(N-1)]

        #L = len(LFSR frame) = upper bound on deg(C)
        L = 0
        #m is index of last change
        m = -1
        
        #n = bit we are correcting.
        for n in range(N):
            #calculate discrepancy from LFSR frame
            d = 0
            for i in range(L+1):
                d ^= (c[i] & seq[n-i])

            #handle discrepancy if needed
            if d != 0:
                
                #store copy of C
                temp = c[:]

                #c = c - x^(n-m)b
                shift = n-m
                for i in range (shift, N):
                    c[i] ^= b[i - shift]

                #if 2L <= n, then the polynomial is unique
                #it's safe to update
                if 2*L <= n:
                    L = n + 1 - L
                    b = temp
                    m = n
        
        #calculate inital state:
        s = []
        for i in range(L):
            s_i = 0
            for j in range(i+1):
                s_i ^= (seq[i-j] & c[j])
            s.append(s_i)
        s = BitVector(bitlist = s[::-1])
        return s, Galois(L, c[:L+1])

    @classmethod 
    def fromReg(self, F, bit):
        numIters = 2*F.size + 4
        seq = [state[bit] for state in F.run(numIters)]
        return Galois.fromSeq(seq)

