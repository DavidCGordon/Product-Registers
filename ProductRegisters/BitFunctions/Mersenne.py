from BitVector import BitVector
from .BitFunction import BitFunction

class Mersenne(BitFunction):
    def __init__(self, size, hexStr, update = None):
        self.size = size
        # here U denotes update polynomial, P denotes primitive polynomial 

        #Format U & P:-----------------------------------------
        #U format: [1,0,0,1,1] = x^4 + x^3 + 1
        if not update: update = [0,1] + [0]*(size-2)
        #covert update: [1,0,0,1,1] -> [0,3,4] -> powers of x
        update = [idx for (idx, t) in enumerate(update) if t == 1]

        #get P from koopman string:
        #P format: "12" -> "10010" (koopman format)-> "100101" (x^5 + x^2 + 1) -> [3,5] 
        primePolynomial = list(BitVector(intVal = int(hexStr, 16), size = size))[1:] + [1]
        primePoly = [(idx+1) for (idx, t) in enumerate(primePolynomial) if t == 1]

        #represent multiplication in GF(2^n): -----------------
        #Shift U copies (multiply by update polynomial)
        bitValues = [[] for _ in range(2*size + 1)]
        for i in range(size):
            for u in update:
                bitValues[i+u].append([i])
        
        #Shift P copies (mod by prime polynomial)
        for i in range(2*size-2, size-1, -1):
            for p in primePoly:
                for x in bitValues[i]:
                    #adding mod 2
                    if x in bitValues[i-p]:
                        bitValues[i-p].remove(x)
                    else:
                        bitValues[i-p].append(x)

        self.fn = [sorted(bitFunc) for bitFunc in bitValues[:size]]