from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.BooleanLogic import BooleanFunction, ANF_spec_repr

from random import randint, sample

# The CrossJoin bitfunction implements the concepts outlined in Elena Dubrova's papers
# (https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6290394) and 
# (https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5290281) 
# to create a scalable sequence generator

class CrossJoin(FeedbackFunction):
    def __init__(self, size, primitive_poly, \
                 nonlinear = True, maxAnds = 4, tapDensity = .75):
                 
        #convert koopman string into polynomial:
        if type(primitive_poly) == str:
            primitive_poly = [int(x) for x in format(int(primitive_poly,16), f"0>{size}b")] + [1]
            #primitive_poly = list(BitVector(intVal = int(primitive_poly, 16), size = size)) + [1]
            
        self.primitive_polynomial = primitive_poly

        #same as fibonacci ANF generation:
        top_fn = [[size-idx] for (idx, t) in enumerate(primitive_poly) if t == 1][::-1][:-1]
        self.fn_list = [[[(i+1)%size]] for i in range(size-1)] + [top_fn]
        self.fn_list = [BooleanFunction.construct_ANF(bitFn) for bitFn in self.fn_list]

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

            self.fn_list[idxA].remove(term)
            self.fn_list[idxB] += ANF_spec_repr([newTerm])

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
        self.fn_list[self.size - 1].add(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx1)

        #add second copy
        self.fn_list[self.size - 1].add(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx2)

        #fix edge case caused by cancelations
        if idx1 == self.size-1:
             self.fn_list[self.size-1] += ANF_spec_repr([newTerm])


    def generateNonlinearity(self, maxAnds, tapDensity):
        self.tau = int(tapDensity * self.size)

        #shift linear terms:
        self.shiftTerms(
            [term for term in self.fn_list[self.size-1] if tuple(term)[0] > self.tau],
             self.size-1, self.tau
        )

        #may be slow, but is easier to read.
        tapped = set()
        while len(tapped) < self.tau:
            self.addNonLinearTerm(maxAnds)
            
            for fn in self.fn_list[:self.tau]:
                for term in fn:
                    tapped |= term
        return
