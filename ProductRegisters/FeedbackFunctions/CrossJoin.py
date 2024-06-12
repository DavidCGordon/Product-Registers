from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.BooleanLogic import BooleanFunction, ANF_spec_repr, AND, XOR, VAR, CONST

from ProductRegisters.Tools.RootCounting.MonomialProfile import TermSet,MonomialProfile
from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet
from ProductRegisters.Tools.RootCounting.RootExpression import RootExpression

from random import randint, sample

# The CrossJoin bitfunction implements the concepts outlined in Elena Dubrova's papers
# (https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6290394) and 
# (https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5290281) 
# to create a scalable sequence generator


# For Now: ONLY SUPPORTS ANF TERMS:
class CrossJoin(FeedbackFunction):
    def __init__(self, size, primitive_poly):
        """
        Constructor for the CrossJoin family of function families.

        Args:
          size: 
            the number of bits in the feedback register
          primitive_poly:
            the primitive polynomial of the base LFSR. This can be either in koopman string format,
            or given as a list with 0/1 entries, where p[i] is the coeffient of x^i in the primitive polynomial.
        Returns:
            None
        """
        # convert koopman string into polynomial:
        if type(primitive_poly) == str:
            primitive_poly = [int(x) for x in format(int(primitive_poly,16), f"0>{size}b")] + [1]
            
        self.primitive_polynomial = primitive_poly

        # same as fibonacci ANF generation:
        top_fn = [[size-idx] for (idx, t) in enumerate(primitive_poly) if t == 1][::-1][:-1]
        self.fn_list = [[[(i+1)%size]] for i in range(size-1)] + [top_fn]
        self.fn_list = [
            XOR(BooleanFunction.construct_ANF(bitFn),)
            for bitFn in self.fn_list
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

    @property
    def linear_feedback(self):
        return [f.args[0] for f in self.fn_list]

    @property
    def monomial_feedback(self):
        output = []
        for f in self.fn_list:
            if len(f.args) > 1:
                output.append(XOR(*f.args[1:]))
            else:
                output.append(CONST(0))
        return output
    
    def compensation_list(self):
        comp_list = []
        curr_fn = CONST(0)
        for fn in self.monomial_feedback[::-1]:
            curr_fn = curr_fn.shift_indices(-1)

            # moving in/out of ANF_spec_repr to get nicer cancellation
            # slightly inefficient, but this is still relatively fast
            # and is much cleaner to read/write
            curr_fn = ANF_spec_repr.from_BooleanFunction(curr_fn)
            curr_fn ^= ANF_spec_repr.from_BooleanFunction(fn)
            curr_fn = curr_fn.to_BooleanFunction()

            comp_list.append(curr_fn)
        
        # list is buit in reverse, so reverse when returning:
        return comp_list[::-1]
        

    def root_expressions(self):
        REs = []
        compensation_list = self.compensation_list()
        
        for bit in range(self.size):
            term_lengths = [len(term.args) for term in compensation_list[bit].args if type(term) != CONST]
            count = max(term_lengths, default = 1)
            count = min(count,self.size)
            REs.append(RootExpression([JordanSet({self.size:count},1)]))
        return REs

    def filter_generator(self):
        pass