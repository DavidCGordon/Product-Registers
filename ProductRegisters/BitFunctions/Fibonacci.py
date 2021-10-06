from BitVector import BitVector
from .BitFunction import BitFunction
from copy import deepcopy

class Fibonacci(BitFunction):
    def __init__(self, size, primitive_poly):
        #convert koopman string into polynomial:
        if type(primitive_poly) == str:
            primitive_poly = list(BitVector( \
                intVal = int(primitive_poly, 16), size = size)) + [1]
        self.primitive_polynomial = primitive_poly

        #[1,0,1,1] -> [0,2,3] -> [3,1,0] -> [0,1,3] -> [0,1]
        top_fn = [[size-idx] for (idx, t) in enumerate(primitive_poly) if t == 1][::-1][:-1]
        self.fn = [[[i+1]] for i in range(size-1)] + [top_fn]
        
        self.size = size        
        self._inverted = False
        
    def invert(self):
        if not self._inverted:
            #reverse taps:
            self.fn[self.size-1] = [[0]] + \
                [[self.size-i[0]] for i in self.fn[-1][1:]][::-1]
            #flip bit labelling:
            self.flip()
            self._inverted = True
        else:
            #flip bit labelling:
            self.flip()
            #reverse taps:
            self.fn[self.size-1] = [[0]] + \
            [[self.size-i[0]] for i in self.fn[-1][1:]][::-1]
            self._inverted = False
            
    #use berlekamp-massey algorithm to generate a Fibonacci LFSR
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
                    
        #return (initial_state of register,fibonacci bitFn)
        init_state = BitVector(bitlist = seq[:L][::-1])
        return init_state, Fibonacci(L, c[:L+1])
            
