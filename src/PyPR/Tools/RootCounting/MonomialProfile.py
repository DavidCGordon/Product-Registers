from PyPR.Tools.RootCounting.Combinatorics import choose, binsum
from PyPR.Tools.RootCounting.OverlappingRectangle import rectangle_solve
from PyPR.Tools.RootCounting.PartialOrders import maximalElements
from PyPR.BooleanLogic import AND, XOR, CONST, VAR
from itertools import product, combinations, tee, cycle

from math import comb

# A list of Blocks, and the corresponding weight
# Completely ignores constants (I hope that works)
class TermSet:
    def __init__(self,totals,counts):
        self.totals = totals
        self.counts = counts 
    
    def __copy__(self):
        return TermSet(
            {k:v for k,v, in self.totals.items()},
            {k:v for k,v, in self.counts.items()}
        )

    def __mul__(self,other):
        output = TermSet({},{})
        for block_id in self.totals:
            if block_id in output.totals:
                output.counts[block_id] += self.counts[block_id]
            else:
                output.counts[block_id] = self.counts[block_id]
                output.totals[block_id] = self.totals[block_id]
                
        for block_id in other.totals:
            if block_id in output.totals:
                output.counts[block_id] += other.counts[block_id]
            else:
                output.counts[block_id] = other.counts[block_id]
                output.totals[block_id] = other.totals[block_id]
        
        # reduce the counts in each multiplication
        for block_id in output.totals:
            output.counts[block_id] = min(output.counts[block_id],output.totals[block_id])

        return output

    def __str__(self):
        order = sorted(self.totals.keys(), key = lambda x: self.totals[x], reverse=True)
        return  "<" +  ", ".join(f"{block_id}:{self.counts[block_id]}/{self.totals[block_id]}" for block_id in order) + ">"










def isMonomialSubset(a,b):
    # used in MonomialProfile.__mul__
    # Don't exclude terms with smaller bases:
    # Because a termset represents terms with at least 1 variable (1 up to count)
    # but not zero, you can't have a term with a larger basis consume the smaller
    # ones without losing track of those monomials.
    if a.totals.keys() != b.totals.keys():
        return False
    
    # all counts in A must be < B to be a subset.
    for block_id in a.totals:
        compare_value = b.counts[block_id] if block_id in b.counts else 0
        if a.counts[block_id] > compare_value:
            return False
    return True



class MonomialProfile:
    # wrapper for TermSets, with extra functionality:
    # self.terms is a set of TermSets.
    def __init__(self,term_list = None):
        if term_list is None:
            self.terms = set()
        else:
            self.terms = set(term_list)

    @classmethod             
    def from_merged(cls, fn_list, blocks):
        total_len = sum(len(block) for block in blocks)
        bitmap = [
            MonomialProfile.logical_zero() 
            for i in range(total_len)
        ]
        
        for block_id in range(len(blocks)):
            for bit in blocks[block_id]:
                bitmap[bit] = MonomialProfile([TermSet(
                    totals={block_id:len(blocks[block_id])},
                    counts={block_id:1}
                )])

        total_fn = XOR(*fn_list)
        total_fn = total_fn.remap_constants([
            (0, MonomialProfile.logical_zero()),
            (1, MonomialProfile.logical_one())
        ])

        return total_fn.eval_ANF(bitmap)
    
    def to_BooleanFunction(self):
        output = XOR()
        for term in self.terms:
            if len(term.counts) == 0:
                output.add_arguments(CONST(1))
                continue

            term_fn = AND()    
            for block,count in term.counts.items():
                term_fn.add_arguments(*([VAR(block)] * count))
            output.add_arguments(term_fn)
        return output
    

    

    def __str__(self):
        termlist = sorted(
            list(self.terms),
            key = lambda x: (
                tuple(sorted(zip(x.totals.values(),x.counts.values())))
            )
        )

        return " + ".join(str(term) for term in termlist)
    
    def __copy__(self):
        return MonomialProfile([termset.__copy__() for termset in self.terms])


    def __xor__(self, other): return self.__add__(other)
    def __add__(self, other):
        #clean out redundant subsets and merge.
        new_terms = maximalElements(
            leq_ordering=isMonomialSubset, 
            inputs=[self.terms, other.terms]
        )

        return MonomialProfile(new_terms)

    def __and__(self, other): return self.__mul__(other)
    def __mul__(self, other):

        new_terms = [a*b for a, b in product(self.terms, other.terms)]

        #clean out redundant subsets.
        new_terms = maximalElements(
            leq_ordering = isMonomialSubset,
            inputs = [new_terms]
        )

        return MonomialProfile(new_terms)

    # When multiplying by Logical One, you should leave the result untouched
    # When adding with Logical One, you should add an indicator term 
    # These effects are accomplished by the MonomialProfile with an Empty TermSet
    @classmethod
    def logical_one(cls): return MonomialProfile([TermSet({},{})])

    # When multiplying by Logical Zero, you should cancel out terms
    # When adding with Logical zero, you should leave the result untouched
    # These effects are accomplished by the empty MonomialProfile
    @classmethod
    def logical_zero(cls): return MonomialProfile([])

    # The Monomial Profile of an inverted function is the same
    def __invert__(self): return self ^ MonomialProfile.logical_one()


    def upper(self):
        # initialize
        num_monomials = 0
        basis_table = {}

        # build basis table
        for termset in self.terms:
            basis = tuple(sorted((termset.totals.keys())))
            values = tuple([binsum(termset.totals[id],termset.counts[id]) for id in basis])
            
            # handle empty monomial profile:
            if basis == ():
                basis_table[basis] = [(1,)]
                continue

            if basis in basis_table:
                basis_table[basis].append(values)
            else:
                basis_table[basis] = [values]

        # evaluate the basis table using hyperrec algorithm
        for basis, rectangle_list in basis_table.items():
            num_monomials += rectangle_solve(rectangle_list)
        return num_monomials








    # for cube attacks
    def get_cube_candidates(self):
        candidates =  []
        for term_set in self.terms:
            for block_id in term_set.totals:
                # create the candidate:
                modified_set = term_set.__copy__()
                modified_set.counts[block_id] -= 1

                # test if the candidate is useful
                useful = True
                for other in self.terms:
                    # dont compare to the termset this one was derived from
                    if other == term_set:
                        continue

                    # all counts in A must be <= B to be a subset.
                    is_subset = True
                    for t in modified_set.totals:
                        compare_value = other.counts[t] if t in other.counts else 0
                        if modified_set.counts[t] > compare_value:
                            is_subset = False
                            break

                    if is_subset:
                        useful = False
                        break

                if useful:
                    num_cubes = 1
                    for i in modified_set.totals:
                        num_cubes *= choose(
                            modified_set.totals[i],
                            modified_set.counts[i]
                        )

                    # for every bit which is in the block, but not in the term we are testing
                    # there is a coin flip on whether its full monomial (i.e. test monomial * bit)
                    # appears. This is the chance that all of those monomials fail to appear.
                    cube_success_rate = 2**(-(
                        term_set.totals[block_id]-term_set.counts[block_id]
                    ))

                    candidates.append((
                        modified_set,
                        block_id,
                        num_cubes,
                        cube_success_rate
                    ))
                        
        return candidates









    def get_monomials(self,complete_subsets=False):
        if complete_subsets:
            return self._get_monomials_complete()
        else:
            return self._get_monomials_exact()

    def _get_monomials_exact(self):
        # construct most general totals matrix:
        total_dim = 0
        for term in self.terms:
            total_dim = max([total_dim, *term.totals.keys()])
        total_dim += 1
        
        totals = [0 for i in range(total_dim)]

        for term in self.terms:
            for k,v in term.totals.items():
                totals[k] = v

        # build basis table
        basis_table = {}
        for termset in self.terms:
            basis = tuple(sorted((termset.totals.keys())))
            values = tuple([termset.counts[id] for id in basis])
            
            # handle empty monomial profile:
            if basis == ():
                continue

            if basis in basis_table:
                basis_table[basis].append(values)
            else:
                basis_table[basis] = [values]

        # perform rollover for each basis:
        for basis, rects in basis_table.items():
            curr_vec = [1 for i in range(len(basis))]
            curr_vec[0] = 0

            # sort rectangles into correct order for iteration
            rects = sorted(rects, key = lambda x: x[0], reverse=True)
            for d in range(1,len(basis)):
                rects = sorted(rects, key = lambda x: x[d])

            rect_idx = 0
            while rect_idx < len(rects):
                # increment degree
                curr_vec[0] += 1


                # rollover loop:
                rollover_idx = 0
                rollover_copy = [x for x in curr_vec]
                # loop until indices are not too large:
                while any((
                    (rollover_copy[i] > rects[rect_idx][i])
                    for i in range(len(rects[rect_idx]))
                )):
                    # normal rollover for everything but last place:
                    if rollover_idx < len(basis)-1:
                        rollover_copy[rollover_idx] = 1
                        rollover_copy[rollover_idx + 1] += 1
                        rollover_idx += 1

                    # if necessary swap to next rect & reset rollover attempt:
                    elif rollover_idx == len(basis)-1:
                        rollover_copy = [x for x in curr_vec]
                        rollover_idx = 0
                        rect_idx += 1

                        if rect_idx == len(rects):
                            break

                # copy rollover back into curr_vec and begin the combinatorics:
                curr_vec = rollover_copy
                if rect_idx < len(rects):
                    # set up basic combinations objects
                    comb_iters = [combinations(range(totals[i]), 0) for i in range(total_dim)]

                    # insert the unique ones for this degree combination
                    for i in range(len(basis)):
                        comb_iters[basis[i]] = combinations(range(totals[basis[i]]), curr_vec[i]) 
                    
                    # combine into a product iterator and yield
                    total_iter = iproduct(*comb_iters)
                    for item in total_iter:
                        yield item

    def _get_monomials_complete(self):
        # convert counts to rectangle list
        rects = [rect(term) for term in self.terms]
        
        # construct most general totals matrix:
        dim = max([0] + [len(r) for r in rects])
        totals = [0 for i in range(dim)]
        for term in self.terms:
            for k,v in term.totals.items():
                totals[k] = v

        # rounds of stable sorting to get rectangles in correct order
        rects = sorted(rects, key = safe_get(0), reverse=True)
        for d in range(1,dim):
            rects = sorted(rects, key = safe_get(d))

        # single rollover loop
        rect_idx = 0
        curr_vec = [0 for i in range(dim)]
        while rect_idx < len(rects):

            # increment degree
            curr_vec[0] += 1

            
            # rollover loop:
            rollover_idx = 0
            rollover_copy = [x for x in curr_vec]

            # as long as any index is too large:
            while any((
                (rollover_copy[i] > rects[rect_idx][i])
                for i in range(len(rects[rect_idx]))
            )):
                # normal rollover for everything but last place:
                if rollover_idx < len(rects[rect_idx])-1:
                    rollover_copy[rollover_idx] = 0
                    rollover_copy[rollover_idx + 1] += 1
                    rollover_idx += 1

                # if necessary swap to next rect & reset rollover attempt:
                elif rollover_idx == len(rects[rect_idx])-1:
                    rollover_copy = [x for x in curr_vec]
                    rollover_idx = 0
                    rect_idx += 1

                    if rect_idx == len(rects):
                        break

            # copy rollover back into curr_vec
            curr_vec = rollover_copy

            if rect_idx < len(rects):
                comb_iters = [combinations(range(totals[i]), curr_vec[i]) for i in range(dim)]
                total_iter = iproduct(*comb_iters)

                for item in total_iter:
                    yield item

# lazy product implementation for faster skipping of unusable sets :)
# attribution: https://discuss.python.org/t/a-product-function-which-supports-large-infinite-iterables/5753
def iproduct(*iterables, repeat=1):
    iterables = [item for row in zip(*(tee(iterable, repeat) for iterable in iterables)) for item in row]
    N = len(iterables)
    saved = [[] for _ in range(N)]  # All the items that we have seen of each iterable.
    exhausted = set()               # The set of indices of iterables that have been exhausted.
    for i in cycle(range(N)):
        if i in exhausted:  # Just to avoid repeatedly hitting that exception.
            continue
        try:
            item = next(iterables[i])
            yield from product(*saved[:i], [item], *saved[i+1:])  # Finite product.
            saved[i].append(item)
        except StopIteration:
            exhausted.add(i)
            if not saved[i] or len(exhausted) == N:  # Product is empty or all iterables exhausted.
                return
    yield ()  # There are no iterables.

def safe_get(d):
    def f(x):
        if d < len(x):
            return x[d]
        else:
            return 0
    return f

def rect(term):
    out = [0 for i in range(1+max(
        # add 0 to handle empty terms
        list(term.totals.keys()) + [0]
    ))]

    for k,c, in term.counts.items():
        out[k] = c
    return tuple(out)