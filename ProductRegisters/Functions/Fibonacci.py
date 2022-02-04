from BitVector import BitVector

from ProductRegisters.Functions import FeedbackFunction
from ProductRegisters import ANF

from copy import deepcopy

class Fibonacci(FeedbackFunction):

    def _anf_from_poly(self,size,primitive_polynomial):
        #Create the top function:
        # Example: [1,0,1,1] -> [0,2,3] -> [[3],[1],[0]]
        topFn = [[size-idx] for (idx, t) in enumerate(primitive_polynomial) if t == 1]
        # Example: [[3],[1],[0]] -> [[0],[1],[3]]
        topFn = topFn[::-1]
        # Example: [[0],[1],[3]] -> [[[0],[1]]]
        topFn = [topFn[:-1]]

        #create the shift for all other bits
        shiftFn = [[[i+1]] for i in range(size-1)]

        return [ANF(bitFn) for bitFn in (shiftFn + topFn)]

    def __init__(self, size, primitive_polynomial):
        #convert koopman strings
        if type(primitive_polynomial) == str:
            primitive_polynomial = list(BitVector( \
                intVal = int(primitive_polynomial, 16), size = size)) + [1]
        
        #assign attributes
        self.primitive_polynomial = primitive_polynomial
        self.size = size        
        self.is_inverted = False

        self.anf = self._anf_from_poly(size,primitive_polynomial)

    def invert(self):
        if not self.is_inverted:
            #reverse taps:
            self.anf = self._anf_from_poly(self.size,self.primitive_polynomial[::-1])
            #flip bit labelling:
            self.flip()
            self.is_inverted = True
        else:
            #remake anf:
            self.anf = self._anf_from_poly(self.size,self.primitive_polynomial)
            self.is_inverted = False
            
    #use berlekamp-massey algorithm to generate a Fibonacci LFSR
    @classmethod
    def fromSeq(self,seq):
        #N = total number of bits to process
        N = len(seq)
        #c = current connection polynomial guess
        c = [1] + [0 for i in range(N-1)]
        #b = prev. connection polynomial guess
        b = [1] + [0 for i in range(N-1)]

        #L = len(LFSR frame) = upper bound on deg(C)
        L = 0
        #m is the index of last change
        m = -1
        
        #n = bit we are correcting.
        for n in range(N):

            #calculate discrepancy from LFSR frame
            d = 0
            for i in range(L+1):
                d ^= (c[i] & seq[n-i])

            #handle discrepancy (if needed)
            if d != 0:
                
                #store copy of C
                temp = c[:]

                #c = c - (x**(n-m) * b)
                shift = n-m
                for i in range (shift, N):
                    c[i] ^= b[i - shift]

                #if 2L <= n, then the polynomial is unique
                #it's safe to update the linear complexity.
                if 2*L <= n:
                    L = n + 1 - L
                    b = temp
                    m = n
                    
        #return (initial_state of register,fibonacci bitFn)
        init_state = BitVector(bitlist = seq[:L][::-1])
        return init_state, Fibonacci(L, c[:L+1])

    @classmethod 
    def fromReg(self,F, bit):
        numIters = 2*F.size + 4
        seq = [state[bit] for state in F.run(numIters)]
        return Fibonacci.fromSeq(seq)
