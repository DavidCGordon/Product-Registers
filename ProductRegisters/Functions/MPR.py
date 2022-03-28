from BitVector import BitVector

from ProductRegisters.Functions import FeedbackFunction
from ProductRegisters.ANF import ANF

from ProductRegisters import FeedbackRegister
from ProductRegisters.Tools.BerlekampMassey import berlekamp_massey

from functools import cached_property

class MPR(FeedbackFunction):
    def __init__(self, size, primitive_poly, update_poly = None):
        self.size = size

        # Format Update & Primitive polyomials: -----------------------------

        # convert U to a list:
        if not update_poly:
            update_poly = [0,1] + [0]*(size-2)
        elif type(update_poly) == int:
            update_poly = list(BitVector(intVal = update_poly, size = size))[::-1]

        # covert update to powers ([1,0,0,1,1] -> [0,3,4])
        self.update_polynomial = update_poly
        update_powers = [idx for (idx, t) in enumerate(update_poly) if t == 1]


        # P can be either polynomial list or koopman hex string:
        #   -koopman format note: binary interpretation is missing final 1
        #   -example: "12" -> "10010" -> "100101" = (1 + x^3 + x^5) -> [0,3,5]
        if type(primitive_poly) == str:
            primitive_poly = list(BitVector(intVal = int(primitive_poly, 16), size = size))+[1]

        self.primitive_polynomial = primitive_poly
        primitive_powers = [(idx) for (idx, t) in enumerate(primitive_poly) if t == 1]

        
        #represent multiplication in GF(2^n): ----------------------------

        #the anf also includes n-1 "hypothetical bits" for higher powers
        #anf[:size] are real bits, anf [size:] are hypothetical bits
        functions = [[] for i in range(2*size - 1)]
        
        # Multiply by update polynomial U
        for idx in range(size):
            for power in update_powers:
                functions[idx+power].append([idx])

        #convert to ANF objects (easiest here):
        functions = [ANF(func) for func in functions]

        # Mod by primitive polynomial P
        for idx in range(2*size-2, size-1, -1):
            for power in primitive_powers[:-1]:

                #Xor the ANFs
                functions[idx - size + power] += functions[idx]

        #return real bits
        self.anf = functions[:size]

    @cached_property
    def minimal_polynomial(self):
        #create a feedback register to determine minimal polynomial
        seq = []
        testReg = FeedbackRegister(2**self.size-1,self)
        for state in testReg.run((self.size+1)*2):
            seq.append(state[0])
        _, m = berlekamp_massey(seq)
        return m