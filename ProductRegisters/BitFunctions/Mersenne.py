from BitVector import BitVector
from .BitFunction import BitFunction

class Mersenne(BitFunction):
    def __init__(self, size, primitive_poly, update_poly = None):
        self.size = size

        #Format Update & Primitive polyomials: -----------------------------

        #convert U to a list:
        if not update_poly:
            update_poly = [0,1] + [0]*(size-2)
        elif type(update_poly) == int:
            update_poly = list(BitVector(intVal = update_poly, size = size))[::-1]

        #covert update to powers ([1,0,0,1,1] -> [0,3,4])
        self.update_polynomial = update_poly
        update_powers = [idx for (idx, t) in enumerate(update_poly) if t == 1]
        

        #P can be either polynomial list or koopman hex string:
        #   -koopman format note: binary interpretation is missing final 1
        #   -example: "12" -> "10010" -> "100101" = (1 + x^3 + x^5) -> [0,3,5]
        if type(primitive_poly) == str:
            primitive_poly = list(BitVector(intVal = int(primitive_poly, 16), size = size))+[1]

        self.primitive_polynomial = primitive_poly
        primitive_powers = [(idx) for (idx, t) in enumerate(primitive_poly) if t == 1]


        #represent multiplication in GF(2^n): ----------------------------

        #the anf also includes n-1 "hypothetical bits" for higher powers
        #anf[:size] is real bits, anf [size:] is hypothetical
        anf = [[] for i in range(2*size - 1)]
        
        #Shift U copies (multiply by update polynomial)
        for idx in range(size):
            for power in update_powers:
                anf[idx+power].append([idx])

        #Shift P copies (mod by prime polynomial)
        for idx in range(2*size-2, size-1, -1):
            for power in primitive_powers[:-1]:

                #list xor
                for term in anf[idx]:
                    if term in anf[idx - size + power]:
                        anf[idx-size+power].remove(term)
                    else:
                        anf[idx-size+power].append(term)

        #return real bits
        self.fn = [sorted(bitFn) for bitFn in anf[:size]]
