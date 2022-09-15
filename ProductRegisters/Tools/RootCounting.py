from itertools import product
from itertools import groupby

from functools import cached_property

#class only used internally for CMPR algorithm at the moment

class RootExpression:
    # ANF FORMAT:
    # XOR: xors are handled as a list of terms
    # AND: terms are handled as a dict of basis->count
    # GF(2): GF(2) values are ignored, and are handled outside.

    def __init__(self,anf = None):
        if anf is None:
            anf = []
        self.anf = anf
    
    def __str__(self):
        " + ".join(
            ("<" + ", ".join(f"{k}:{v}" for k,v in term.items()) + ">")
            for term in self.anf
        )

    def __add__(self, other):
        #clean out redundant subsets and merge.
        new_anf = maximalElements(
            leq_ordering=isSubset, 
            input_lists=[self.anf, other.anf]
        )

        return RootExpression(new_anf)

    # handle foiling multiplication
    def __mul__(self, other):
        #foil
        new_term_lists = []
        for term_a, term_b in product(self.anf, other.anf):
            new_term_lists.append([merge_terms(term_a,term_b)])

        #clean out redundant subsets.
        new_anf = maximalElements(
            leq_ordering = isSubset,
            input_lists = new_term_lists
        )

        return RootExpression(new_anf)


    #calculate upper bound (ignoring bases in locked list)
    def upper(self, locked_list = None):

        #clean/remove embedded SUBsets:
        cleaned_terms = maximalElements(
            leq_ordering = isEmbeddedSubset,
            input_lists = [self.anf]
        )

        #create basis table:
        basis_table = createBasisTable(
            term_list = cleaned_terms,
            evaluation_function = optimistic_evaluation,
            locked_list = locked_list
        )

        #evaluate basis table using hyperrec algorithm:
        linear_complexity = 0
        for x in basis_table.values():
            linear_complexity += rectangle_solve(x)
        return linear_complexity


    #calculate lower bound (ignoring bases in locked list)
    def lower(self, locked_list = None, safety_factor = None):
        
        #clean/remove embedded sets:
        cleaned_terms = maximalElements(
            leq_ordering = isEmbeddedSet,
            input_lists = [self.anf]
        )

        #create basis table:
        basis_table = createBasisTable(
            term_list = cleaned_terms,
            evaluation_function = pessimistic_expected_value,
            locked_list = locked_list,
        )

        #evaluate basis table using hyperrec algorithm:
        linear_complexity = 0
        for x in basis_table.values():
            linear_complexity += rectangle_solve(x)
        return linear_complexity

        
#convert pairs to numerical values for the rectangle algorithm
def createBasisTable(term_list,evaluation_function,locked_list):
    basis_table = {}
    for term in term_list:
        converted = [] # stores binsum(basis,count) for the term
        basis = []     # stores the bases for the term

        # sort the term for consistent ordering of dimensions.
        ordered = sorted(term.items(),key = lambda x:x[0], reverse = True)

        locked = False

        for b,c in ordered:
            
            # don't add terms whose basis is locked:
            if locked_list and b in locked_list:
                locked = True
                break

            # different flags for different computations:
            converted.append(evaluation_function(b,c))
            basis.append(b)

        # if this term was not locked, append it to the list.
        if not locked:
            basis = tuple(basis)
            if basis in basis_table:
                basis_table[basis].append(tuple(converted))
            else:
                basis_table[basis] = [tuple(converted)]
    return basis_table


#different statistics to change the way the bounds are given:

def optimistic_evaluation(b,c):
    return binsum(b,c)

def pessimistic_evaluation(b,c):
    if b == c: 
        return binsum(b,c-1)
    else:
        return binsum(b,c)

def pessimistic_expected_value(b,c):
    if b == c: 
        return (binsum(b,c-1) * (2**b-1)**2) / (2**(2*b))
    else:
        return (binsum(b,c) * (2**b-1)**2) / (2**(2*b))

""" 
WORK IN PROGRESS DO NOT USE:
 
def confidence_interval_lower(b,c,base_function, z = 2.58):
    mean = base_function(b,c) * (2**b-1)**2) / (2**(2*b))
    deviation = sqrt(((2**b-1) / 2**(2*b)) / (base_function(b,c) / b))
    return mean - z*deviation

def confidence_interval_upper(b,c,base_function, z = 2.58):
    mean = base_function(b,c) * (2**b-1)**2) / (2**(2*b))
    deviation = sqrt(((2**b-1) / 2**(2*b)) / (base_function(b,c) / b))
    return mean + z*deviation
"""


# Union (AND) of terms, combining counts as needed
# returns a new dict (new term).
def merge_terms(*terms):
    new_term = {}
    
    # sum all counts in the same basis
    # uses a table for efficient lookup
    for term in terms:
        for basis,count in term.items():
            
            if basis in new_term:
                new_term[basis] += count
            else:
                new_term[basis] = count

    # reduce the sums so no value > key
    for basis in new_term:
        new_term[basis] = min(basis,new_term[basis])
        
    return new_term

#HELPERS FOR EVALUATION:
def choose(n,k):
    prod = 1
    for i in range(k):
        prod *= (n-i)/(k-i)
    return round(prod)

def binsum(n,d):
    tot = 0
    for k in range(1,d+1):
        tot += choose(n,k)
    return tot

"""
TERM RELATIONSHIPS ------------------------------------------------------------------------------------------

Math Note:
all three methods "is____" methods define partial orders;
this makes the O(n^2) single pass "kick out" algorithm correct
    - it is also likely optimal, unfortunately.
"""

#leq_ordering(A,B) returns true if A <= B
def maximalElements(leq_ordering, input_lists):

    maximal_list = []
    for input_list in input_lists:
        for term in input_list:
            maximal = True

            for already_added in maximal_list:
                if leq_ordering(already_added, term):
                    maximal_list.remove(already_added)
                    break

                elif leq_ordering(term,already_added):
                    maximal = False
                    break
                
            if maximal:
                maximal_list.append(term)

    return maximal_list


# If the roots in A are a subset of those in B
def isSubset(a,b):
    if a.keys() != b.keys():
        return False

    # all counts in A must be < B to be a subset.
    for k in a.keys():
        if a[k] > b[k]:
            return False

    # if all conditions pass, return true
    return True


# If A might be embedded in B (causing a degeneracy)
# we remove these to avoid duplicates (we only want to count the largest one)
def isEmbeddedSet(a,b):
    a_keyset = set(a.keys())
    b_keyset = set(b.keys())

    # there should be a strict subset relationship on keys:
    if not a_keyset.issubset(b_keyset):
        return False

    # all shared keys should have equal values.
    for k in b_keyset.intersection(a_keyset):
        if a[k] != b[k]:
            return False

    # all non shared keys should be maximum (for each pair (b,c), b should equal c).
    for k in (b_keyset - a_keyset):
        if b[k] != k:
            return False

    # if all conditions pass, return true
    return True

#  If set A is a subset of something that might be embedded in B:
#  we remove these when we assume B will be included.
def isEmbeddedSubset(a,b):
    a_keyset = set(a.keys())
    b_keyset = set(b.keys())

    # there should be a strict subset relationship on keys:
    if not a_keyset.issubset(b_keyset):
        return False

    # all shared keys should have A <= B .
    for k in b_keyset.intersection(a_keyset):
        if a[k] > b[k]:
            return False

    # all non shared keys should be maximum (for each pair (b,c), b should equal c).
    for k in (b_keyset - a_keyset):
        if b[k] != k:
            return False

    # if all conditions pass, return true
    return True




"""
HYPERREC-ALG FOR FASTER EVALUATION --------------------------------------------------------------------------

API:
solve
"""

#group tuples by first dimension, and add an empty "0" layer
def preprocess(rectangle_list):
    output = sorted(rectangle_list, key = lambda x: x[0], reverse=True)
    output = [(k, [x[1:] for x in g]) for k, g in groupby(output, key = (lambda x: x[0]))]
    output.append((0,[]))
    return output

# check if tuple 1 is <= than tuple 2 in every index
# this indicates rectangle A is contained entirely in rectangle B
def rect_sset(rect_A,rect_B):
    return all(a <= b for a,b in zip(rect_A,rect_B))

#insert tuples from a into b, if they are still needed.
def insert_tuples(upper_layer, lower_layer):
    still_needed = []
    for upper_rectangle in upper_layer:
        needed = True

        # if any rectangle in the lower layer needs includes this one:
        # it is not needed for further area calculations
        for lower_rectangle in lower_layer:
            if rect_sset(upper_rectangle,lower_rectangle):
                needed = False
                break
        if needed:
            still_needed.append(upper_rectangle)

    #add the rectangles which are still needed:
    lower_layer += still_needed

#recursive call to solve (grows poorly with dim, but MUCH better with
#number of rectangles: significantly faster/more useable.
def _solve_rec(dimension,rectangle_list):
    # base case:
    if dimension == 1:
        return max(x[0] for x in rectangle_list)
    
    total_area = 0
    processed_list = preprocess(rectangle_list)

    for layer in range(len(processed_list)-1):

        # use the current dimension to get the width of the layer
        layer_width = processed_list[layer][0] - processed_list[layer+1][0]
        # solve the subproblem 1 dimension down to get the cross-sectional area
        cross_section_area = _solve_rec(dimension-1, processed_list[layer][1])
        # multiply to get the volume of the cross-section and add to total area
        total_area += layer_width * cross_section_area

        # insert any still needed rectangles into the next layer:
        # this ensures the cross-sectional area remains correct.
        insert_tuples(processed_list[layer][1], processed_list[layer+1][1])

    return total_area

#wrapper to preprocess dimension for ease of use.
def rectangle_solve(rectangle_list):
    dimension = len(rectangle_list[0])
    return _solve_rec(dimension,rectangle_list)