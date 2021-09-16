from BitVector import BitVector
from .BitFunction import BitFunction
from copy import deepcopy

class Fibonacci(BitFunction):
    def __init__(self, size, topFunc):
        if type(topFunc) == str:
            taps = list(BitVector(intVal = int(topFunc, 16), size = size))
            topFunc = [[(size-idx)%size] for (idx, t) in enumerate(taps) if t == 1]

        self.fn = [[[(i+1)%size]] for i in range(size-1)] + [topFunc]
        self.size = size

    def invert(self):
        #flip_taps
        self.fn[self.size-1] = [[0]] + \
            [[self.size-i[0]] for i in self.fn[self.size-1][1:]][::-1]

        #flip_bit_labelling:
        self.flip()
        
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
        fn = Fibonacci(L, [[L-i] for i in range(1,L+1) if c[i]][::-1])
        init_state = int(BitVector(bitlist = seq[:L][::-1]))
        return init_state, fn
            
