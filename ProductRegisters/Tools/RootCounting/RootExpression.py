from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet
from ProductRegisters.Tools.RootCounting.OverlappingRectangle import rectangle_solve
from ProductRegisters.Tools.RootCounting.EvaluationFunctions import pessimistic_expected_value
from ProductRegisters.Tools.RootCounting.Combinatorics import binsum, powerset


from ProductRegisters.BooleanLogic import CONST, VAR

from itertools import product
import time

def isEmbeddedSubset(a,b):
    a_keyset = set(a.roots.keys())
    b_keyset = set(b.roots.keys())

    # there should be a strict subset relationship on keys:
    if not a_keyset.issubset(b_keyset):
        return False

    # all shared keys should have A <= B .
    for k in b_keyset.intersection(a_keyset):
        if a.roots[k] > b.roots[k]:
            return False

    # all non shared keys should be maximum (for each pair (b,c), b should equal c).
    for k in (b_keyset - a_keyset):
        if b.roots[k] != k:
            return False

    # if all conditions pass, return true
    return True

def isRootSubset(a,b):
    # all counts in A must be <= B to be a subset.
    for k in a.roots.keys():
        if a.roots[k] > b.roots[k]:
            return False
    return True


class RootExpression:
    def __init__(self,root_table):
        self.root_table = root_table
    
    def extend(self, extension_jordan_set):
        max_mult = 0
        for field_tuple in self.root_table:
            for term in self.root_table[field_tuple]:
                if isEmbeddedSubset(extension_jordan_set, term):
                    max_mult = max(max_mult,*term.mults)
    
        field_tuple = tuple(sorted(extension_jordan_set.roots.keys()))
        return self + RootExpression({
            field_tuple: set([JordanSet(extension_jordan_set.roots, set([max_mult + 1]))])
        })

    def __str__(self):
        outstr = "{\n"
        for field_tuple in self.root_table:
            outstr += f"{field_tuple}: " + " + ".join(str(term) for term in self.root_table[field_tuple])
            outstr += "\n"
        outstr += "}"
        return outstr
    
    def __copy__(self):
        new_table = {}
        for field_tuple in self.root_table:
            new_table[field_tuple] = set()
            for JSet in self.root_table[field_tuple]:
                new_table[field_tuple].add(JSet.__copy__())
        return RootExpression(new_table)


    
    def __xor__(self, other): return self.__add__(other)
    def __add__(self, other):
        # clean out redundant subsets and merge.

        # add selfs terms
        output = self.__copy__()
        for field_tuple in other.root_table:
            # if there is no conflict, just write directly:
            if field_tuple not in self.root_table:
                output.root_table[field_tuple] = set()
                for JSet in other.root_table[field_tuple]:
                    output.root_table[field_tuple].add(JSet)
                continue
            
            # otherwise, split mults:
            # problem there could be smaller things, which were kept for mult reasons
            for JSet in other.root_table[field_tuple]:
                new_term = JSet.__copy__()

                # add in the new term:
                to_remove = []
                merge_location = None
                for already_added in output.root_table[field_tuple]:
                    smaller = isRootSubset(new_term, already_added)
                    bigger = isRootSubset(already_added,new_term)

                    if smaller and bigger:
                        merge_location = already_added
                    elif smaller:
                        new_term.mults -= already_added.mults
                        # can break early in this case
                        if not new_term.mults:
                            break
                    elif bigger:
                        already_added.mults -= new_term.mults
                        if not already_added.mults:
                            to_remove.append(already_added)
                
                # update the table entry accordingly:
                if new_term.mults and merge_location is not None:
                    merge_location.mults |= new_term.mults
                elif new_term.mults and merge_location is None:
                    output.root_table[field_tuple].add(new_term)
                for term in to_remove:
                    output.root_table[field_tuple].remove(term)

        #print("\n")
        return output


    def __and__(self, other): return self.__mul__(other)
    def __mul__(self, other):
        # clean out redundant subsets and merge.

        # add selfs terms
        output = RootExpression({})
        for field_tuple1 in self.root_table:
            for field_tuple2 in other.root_table:
                target_field_tuple = tuple(sorted((set(field_tuple1) | set(field_tuple2))))

                # add if needed
                if target_field_tuple not in output.root_table:
                    output.root_table[target_field_tuple] = set()

                for JSet1 in self.root_table[field_tuple1]:
                    for JSet2 in other.root_table[field_tuple2]:
                        prod = JSet1 * JSet2

                        # add in the new term:
                        to_remove = []
                        merge_location = None
                        for already_added in output.root_table[target_field_tuple]:
                            smaller = isRootSubset(prod, already_added)
                            bigger = isRootSubset(already_added,prod)

                            if smaller and bigger:
                                merge_location = already_added
                            elif smaller:
                                prod.mults -= already_added.mults
                                # can break early in this case
                                if not prod.mults:
                                    break
                            elif bigger:
                                already_added.mults -= prod.mults
                                if not already_added.mults:
                                    to_remove.append(already_added)
                        
                        # update the table entry accordingly:
                        if prod.mults and merge_location is not None:
                            merge_location.mults |= prod.mults
                        elif prod.mults and merge_location is None:
                            output.root_table[target_field_tuple].add(prod)
                        for term in to_remove:
                            output.root_table[target_field_tuple].remove(term)
        return output    

    @classmethod
    #def logical_one(self): return RootExpression({(): set([JordanSet({}, {1})])})
    def logical_one(self): return RootExpression({(): set([JordanSet({}, {1})])})

    @classmethod
    def logical_zero(self): return RootExpression({})

    def __invert__(self): return self ^ RootExpression.logical_one()


    """
    Evaluation
    """

    def upper(self):
        # initialize values
        linear_complexity = 0
        basis_table = {}

        # Assume any embedded sets are present:
        #  - this means no cleaning is necessary
        #  - but each embedded term might needs to be added to multiple entries in the basis table.
        #       - because the roots span across several subfields.
        for field_tuple in self.root_table:
            for term in self.root_table[field_tuple]: 
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
                    counts = tuple([max(term.mults)] + [x[1] for x in pairs])

                    if basis in basis_table:
                        basis_table[basis].append(counts)
                    else:
                        basis_table[basis] = [counts]

        # evaluate the basis table using hyperrec algorithm
        for x in basis_table.values():
            linear_complexity += rectangle_solve(x)
        return linear_complexity


    #calculate lower bound (ignoring bases in locked list)
    def lower(self):
        # initalize values:
        linear_complexity = 0        
        basis_table = {}

        # Evaluate only the maximum coset? 
        # cleaned_terms = maximalElements(
        #     leq_ordering = isEmbeddedSet,
        #     inputs = [self.terms]
        # )

        # add each term into the basis table.
        for field_tuple in self.root_table:
            for term in self.root_table[field_tuple]:
                pairs = sorted([
                    (b,pessimistic_expected_value(b,c))
                    for b,c in term.roots.items()
                ])

                basis = tuple([x[0] for x in pairs])
                counts = tuple([max(term.mults)] + [x[1] for x in pairs])

                if basis in basis_table:
                    basis_table[basis].append(counts)
                else:
                    basis_table[basis] = [counts]
                
        #solve basis table using hyperrec algorithm:
        for rectangle_list in basis_table.items():
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