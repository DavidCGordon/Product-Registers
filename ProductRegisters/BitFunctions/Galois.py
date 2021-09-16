from BitVector import BitVector
from .BitFunction import BitFunction
from copy import deepcopy


class Galois(BitFunction):
    #construction
    def _taps_from_poly(self, p):
        #build tap set fn from polynomial
        self.taps = []
        for i in range(self.size):
            if p[i+1]:
                self.taps.append(self.size-1-i)
                
        for tap in self.taps:
            self.fn[tap] += [[self._tap_bit]]
            
    def __init__(self, size, polynomial):
        self.size = size

        #internal use
        self._tap_bit = 0
        self._polynomial = polynomial
        
        #calculate function
        self.fn = [[[i+1]] for i in range(self.size-1)] + [[]]
        self._taps_from_poly(self._polynomial)
        

    def invert(self):
        #remove current taps:
        for tap in self.taps:
            self.fn[tap].remove([self._tap_bit])

        #flip the non-tap bits
        self.flip()
        
        #flip the polynomial and set the tap_bit
        self._polynomial = self._polynomial[::-1]
        self._tap_bit = self.size-1-self._tap_bit
        
        #add back new taps:
        self._taps_from_poly(self._polynomial)
            

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
