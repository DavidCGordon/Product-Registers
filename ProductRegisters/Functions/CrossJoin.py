from BitVector import BitVector

from ProductRegisters.Functions import FeedbackFunction
from ProductRegisters import ANF

from random import randint, sample

class CrossJoin(FeedbackFunction):
    def __init__(self, size, primitive_poly, \
                 nonlinear = True, maxAnds = 4, tapDensity = .75):
        #convert koopman string into polynomial:
        if type(primitive_poly) == str:
            primitive_poly = list(BitVector(intVal = int(primitive_poly, 16), size = size)) + [1]
        self.primitive_polynomial = primitive_poly

        #same as fibonacci ANF generation:
        top_fn = [[size-idx] for (idx, t) in enumerate(primitive_poly) if t == 1][::-1][:-1]
        self.anf = [[[(i+1)%size]] for i in range(size-1)] + [top_fn]
        self.anf = [ANF(bitFn) for bitFn in self.anf]

        self.size = size
        self.tau = self.size-1

        if nonlinear: self.generateNonlinearity(maxAnds, tapDensity)

    def shiftTerms(self, terms, idxA, idxB):
        for term in terms:
            valid = True
            for var in term:
                valid &= (var >= (idxA-idxB))
            if not valid:
                raise ValueError("Invalid shift attempted for term: " + str(term))
            
            newTerm = frozenset([(i - idxA + idxB) for i in term])

            self.anf[idxA].remove(term)
            self.anf[idxB] += ANF([newTerm])

    def getMinDestination(self, term):
        return max((self.size - 1) - min(term), self.tau)

    def getMaxDestination(self, term):
        return min((self.size + self.tau) - (max(term) + 1), self.size - 1)

    def addNonLinearTerm(self,maxAnds):
        minDest = maxDest = 0
        while not (minDest < maxDest):
            numTaps = randint(2,maxAnds)
            newTerm = frozenset(sample(range(1, self.size), numTaps))
            maxDest = self.getMaxDestination(newTerm)
            minDest = self.getMinDestination(newTerm)

        idx1,idx2 = sample(range(minDest,maxDest+1),2)

        #add first copy
        self.anf[self.size - 1].add(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx1)

        #add second copy
        self.anf[self.size - 1].add(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx2)

        #fix edge case caused by set merging
        if idx1 == self.size-1:
             self.anf[self.size-1] += ANF([newTerm])


    def generateNonlinearity(self, maxAnds, tapDensity):
        self.tau = int(tapDensity * self.size)
        tapped = BitVector(intVal = 0, size = self.size)

        #shift linear terms:
        self.shiftTerms(
            [term for term in self.anf[self.size-1] if tuple(term)[0] > self.tau],
             self.size-1, self.tau
        )

        #may be slow, but is easier to read.
        tapped = frozenset()
        while len(tapped) < self.tau:
            self.addNonLinearTerm(maxAnds)
            for fn in self.anf[:self.tau]:
                for term in fn:
                    tapped |= term
        return
