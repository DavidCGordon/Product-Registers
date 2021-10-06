from BitVector import BitVector
from .BitFunction import BitFunction
from random import randint, sample, shuffle

class Induction(BitFunction):
    def __init__(self, size, nonlinear = True, permute = False):
        self.fn = [ [[i], list(range(i))] for i in range(size) ]
        self.fn[0][1] = True
        self.size = size
        self.induction_order = list(range(self.size))

        if nonlinear: self.GenerateNonlinearity()
        if permute:
            shuffle(self.induction_order)
            newAnf = []
            for bitFn in self.fn:
                newFn = []
                for term in bitFn:
                    if type(term) == bool:
                        newFn.append(term)
                    else:
                        newTerm = []
                        for var in term:
                            newTerm.append(self.induction_order[var])
                        newFn.append(newTerm)
                newAnf.append(newFn)
            for i in range(self.size):
                self.fn[self.induction_order[i]] = newAnf[i]
                
                            
    #better sampling eventually, but this works for now
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
