from BitVector import BitVector
from .BitFunction import BitFunction
from random import randint, sample

class CrossJoin(BitFunction):
    def __init__(self, size, hexStr, maxAnds = 3, density = .75):
        taps = list(BitVector(intVal = int(hexStr, 16), size = size))
        LFSRFunc = [[(size-idx)%size] for (idx, t) in enumerate(taps) if t == 1]

        self.fn = [[[(i+1)%size]] for i in range(size-1)] + [LFSRFunc]
        self.size = size
        self.tau = self.size-1

        self.generateNonlinearity(maxAnds, density)

    def __str__(self):
        outstr = ""
        for i in range(len(self.fn)-1,-1,-1):
            outstr += ("Bit " + str(i)
                + ": " + str(self.fn[i][0][0]) + " + ("
                + (" + ".join([str(term) for term in self.fn[i]][1:])) + ")\n")
        return outstr

    def shiftTerms(self, terms, idxA, idxB):
        alteredTerms = []
        for term in terms:
            valid = True
            for i in term:
                valid &= (i >= (idxA-idxB))
            if valid:
                altTerm = [(i - idxA + idxB) for i in term]
                alteredTerms.append(altTerm)
                self.fn[idxA].remove(term)
                self.fn[idxB] += [altTerm]
            else:
                raise ValueError("Invalid shift attempted for term: " + str(term))
        return [idx for term in alteredTerms for idx in term]

    def getMinDestination(self, term):
        return max((self.size - 1) - min(term), self.tau)

    def getMaxDestination(self, term):
        return min((self.size + self.tau) - (max(term) + 1), self.size - 1)

    def addNonLinearTerm(self,maxAnds):
        minDest = maxDest = 0
        while not (minDest < maxDest):
            numTaps = randint(2,maxAnds)
            newTerm = sample(range(1, self.size), numTaps)
            maxDest = self.getMaxDestination(newTerm)
            minDest = self.getMinDestination(newTerm)

        idx1,idx2 = sample(range(minDest,maxDest+1),2)

        self.fn[self.size - 1] += [newTerm]
        self.fn[self.size - 1] += [newTerm[:]]
        
        newTaps = []
        newTaps.append(self.shiftTerms([newTerm], self.size-1, idx1))
        newTaps.append(self.shiftTerms([newTerm], self.size-1, idx2))
        return [tap for sublist in newTaps for tap in sublist]

    def generateNonlinearity(self, maxAnds, density):
        self.tau = int(density * self.size)
        tapped = BitVector(intVal = 0, size = self.size)

        self.shiftTerms(
            [term for term in self.fn[self.size-1] if term[0] > self.tau],
             self.size-1, self.tau
        )

        while tapped.count_bits() < self.tau:
            newTaps = self.addNonLinearTerm(maxAnds)
            for tap in newTaps:
                tapped[tap] = 1
        return