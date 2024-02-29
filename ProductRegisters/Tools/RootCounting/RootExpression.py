from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet
from ProductRegisters.Tools.RootCounting.OverlappingRectangle import rectangle_solve
from ProductRegisters.Tools.RootCounting.EvaluationFunctions import pessimistic_expected_value
from ProductRegisters.Tools.RootCounting.Combinatorics import binsum, powerset
from ProductRegisters.Tools.RootCounting.PartialOrders import \
    isSubspace, isEmbeddedSet, isEmbeddedSubset, isExactSubset, maximalElements

from ProductRegisters.BooleanLogic import CONST, VAR

from itertools import product
import time

class RootExpression:
    # wrapper for JordanSets, with extra functionality:
    # self.terms is a set of JordanSets.
    def __init__(self,term_list = None):
        if term_list is None:
            self.terms = set()
        else:
            self.terms = set(term_list)
    
    def extend(self, extension_jordan_set):
        max_mult = 0
        for term in self.terms:
            if isEmbeddedSubset(extension_jordan_set, term):
                max_mult = max(max_mult,term.m)
        return self + RootExpression([JordanSet(extension_jordan_set.roots, max_mult + 1)])

    def __str__(self):
        return " + ".join(str(term) for term in self.terms)
    
    def __copy__(self):
        return RootExpression([jordanset.__copy__() for jordanset in self.terms])


    def __xor__(self, other): return self.__add__(other)
    def __add__(self, other):
        print(f'Addition Started: {len(self.terms)} x {len(other.terms)} Terms')
        start_time = time.time()

        #clean out redundant subsets and merge.
        new_anf = maximalElements(
            leq_ordering=isExactSubset, 
            inputs=[self.terms, other.terms]
        )

        print(f'Multiplication Finished: {len(new_anf)} Terms - Time: {time.time()-start_time}')
        return RootExpression(new_anf)

    def __and__(self, other): return self.__mul__(other)
    def __mul__(self, other):
        print(f'Multiplication Started: {len(self.terms)} Terms x {len(other.terms)}')
        start_time = time.time()

        new_term_sets = []
        for a, b in product(self.terms, other.terms):
            new_term_sets.append(a * b)

        #clean out redundant subsets.
        new_anf = maximalElements(
            leq_ordering = isExactSubset,
            inputs = new_term_sets
        )

        print(f'Multiplication Finished: {len(new_anf)} Terms - Time: {time.time()-start_time}')
        return RootExpression(new_anf)

    @classmethod
    def logical_one(self): return RootExpression([JordanSet({}, 1)])

    @classmethod
    def logical_zero(self): return RootExpression([])

    def __invert__(self): return self ^ RootExpression.logical_one()


    """
    Evaluation
    """

    def upper(self, locked_list = None):
        # initialize values
        linear_complexity = 0
        basis_table = {}

        # Assume any embedded sets are present:
        #  - this means no cleaning is necessary
        #  - but each embedded term might needs to be added to multiple entries in the basis table.
        #       - because the roots span across several subfields.
        for term in self.terms: 
            full = []           # - All bases where b == c (contains the 1 coset)
            partial = []        # - All bases where b != c (no embedded sets)
            
            # Add each (b,c) pair to the right list:
            for b,c in term.roots.items():
                if b == c: full.append((b,binsum(b,c-1)))
                else: partial.append((b,binsum(b,c)))

            # Use the full list to distribute embedded sets to appropriate entries
            #   - every subset of bases in the full list represents a valid subfield
            for modifier in powerset(full):
                pairs = sorted(list(modifier) + partial)
                basis = tuple([x[0] for x in pairs])
                counts = tuple([term.m] + [x[1] for x in pairs])

                if basis in basis_table:
                    basis_table[basis].append(counts)
                else:
                    basis_table[basis] = [counts]

        # evaluate the basis table using hyperrec algorithm
        for x in basis_table.values():
            linear_complexity += rectangle_solve(x)
        return linear_complexity


    #calculate lower bound (ignoring bases in locked list)
    def lower(self, locked_list = None, safety_factor = None):
        # initalize values:
        linear_complexity = 0        
        basis_table = {}

        # Evaluate only the maximum coset? 
        cleaned_terms = maximalElements(
            leq_ordering = isEmbeddedSet,
            inputs = [self.terms]
        )

        # add each term into the basis table.
        for term in cleaned_terms:
            pairs = sorted([
                (b,pessimistic_expected_value(b,c))
                for b,c in term.roots.items()
            ])

            basis = tuple([x[0] for x in pairs])
            counts = tuple([term.m] + [x[1] for x in pairs])

            if basis in basis_table:
                basis_table[basis].append(counts)
            else:
                basis_table[basis] = [counts]
            
        #solve basis table using hyperrec algorithm:
        for rectangle_list in basis_table.values():
            linear_complexity += rectangle_solve(rectangle_list)
        return linear_complexity


    """    
    ALTERNATE MULTIPLICATION (test speeds/profile as future work).

    def __mul__(self, other):
        new_root_expression = RootExpression()
        for term_a, term_b in product(self.terms, other.terms):
            new_root_expresion += (term_a * term_b)
        return new_root_expression
    """