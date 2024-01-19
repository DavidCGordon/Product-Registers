from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.BooleanLogic import BooleanFunction, ANF_spec_repr, AND, XOR, VAR

from random import randint, sample

# The CrossJoin bitfunction implements the concepts outlined in Elena Dubrova's papers
# (https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6290394) and 
# (https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5290281) 
# to create a scalable sequence generator


# For Now: ONLY SUPPORTS ANF TERMS:
class CrossJoin(FeedbackFunction):
    def __init__(self, size, primitive_poly):
                 
        # convert koopman string into polynomial:
        if type(primitive_poly) == str:
            primitive_poly = [int(x) for x in format(int(primitive_poly,16), f"0>{size}b")] + [1]
            # primitive_poly = list(BitVector(intVal = int(primitive_poly, 16), size = size)) + [1]
            
        self.primitive_polynomial = primitive_poly

        # same as fibonacci ANF generation:
        top_fn = [[size-idx] for (idx, t) in enumerate(primitive_poly) if t == 1][::-1][:-1]
        self.fn_list = [[[(i+1)%size]] for i in range(size-1)] + [top_fn]
        self.fn_list = [
        XOR(
            BooleanFunction.construct_ANF(bitFn),
        ) for bitFn in self.fn_list
        ]

        self.size = size
        self.tau = self.size-1

    def shiftTerms(self, terms, idxA, idxB):
        for term in terms:
            valid = True
            for var in term.args:
                if type(var) != VAR:
                    raise ValueError("types other than VAR are not currently supported")
                valid &= (var.index >= (idxA-idxB))
            if not valid:
                raise ValueError("Invalid shift attempted for term: " + str(term))
            
            newTerm = AND(*(VAR(var.index - idxA + idxB) for var in term.args))

            self.fn_list[idxA].args[-1].remove_arguments(term)
            self.fn_list[idxB].args[-1].add_arguments(newTerm)

    def getMinDestination(self, term):
        return max((self.size - 1) - min(value.index for value in term.args), self.tau)

    def getMaxDestination(self, term):
        return min((self.size + self.tau) - (max(value.index for value in term.args)+1), self.size - 1)

    def addNonLinearTerm(self,maxAnds):
        minDest = maxDest = 0


        while not (minDest < maxDest):
            numTaps = randint(2,maxAnds)
            newTerm = AND(*(VAR(i) for i in sample(range(1, self.size), numTaps)))
            maxDest = self.getMaxDestination(newTerm)
            minDest = self.getMinDestination(newTerm)

        idx1,idx2 = sample(range(minDest,maxDest+1),2)

        # add first copy
        self.fn_list[self.size - 1].args[-1].add_arguments(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx1)

        # add second copy
        self.fn_list[self.size - 1].args[-1].add_arguments(newTerm)
        self.shiftTerms([newTerm], self.size-1, idx2)


    def generateNonlinearity(self, maxAnds = 4, tapDensity = .75):
        self.tau = min(self.tau,int(tapDensity * self.size))

        # add a set of nodes for nonlinear terms:
        for bit in range(self.size):
            self.fn_list[bit].add_arguments(XOR())

        # shift any needed linear terms
        for term in self.fn_list[self.size-1].args[0].args:
            if term.args[0].index >= self.tau:
                self.fn_list[self.size-1].args[0].remove_arguments(term)
                self.fn_list[self.tau].args[0].add_arguments(term.shift_indices(self.tau-self.size+1))
                
        # main loop:
        tapped = set()
        while len(tapped) < self.tau:
            self.addNonLinearTerm(maxAnds)
            
            # for all nonlinear terms above tau
            for fn in self.fn_list[self.tau:]:
                for term in fn.args[-1].args:
                    tapped |= {val.index for val in term.args}

        # strip off any empty nodes
        for bit in range(self.size):
            nonlinear_terms = self.fn_list[bit].args[-1]
            if len(nonlinear_terms.args) == 0:
                self.fn_list[bit].remove_arguments(nonlinear_terms)

        return
