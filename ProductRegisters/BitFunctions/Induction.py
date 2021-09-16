from BitVector import BitVector
from .BitFunction import BitFunction
from random import randint, sample

class Induction(BitFunction):
    def __init__(self, size, randomize = True):
        self.fn = [ [[i], list(range(i))] for i in range(size) ]
        self.fn[0][1] = True
        self.size = size

        if randomize: self.GenerateNonlinearity()
    
    def GenerateNonlinearity(self):
        for bitIdx in range(3,self.size):
            numTerms = randint(1,2**(bitIdx//2))

            addedTerms = 0
            while addedTerms < numTerms:
                termSize = randint(1,bitIdx-1)
                newTerm = sorted(sample(range(bitIdx), termSize))

                if newTerm not in self.fn[bitIdx]:
                    self.fn[bitIdx] += [newTerm]
                    addedTerms += 1