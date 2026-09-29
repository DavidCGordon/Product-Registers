from itertools import combinations, cycle, product, tee

from PyPR.BooleanLogic import AND, CONST, VAR, XOR

from PyPR.Tools.RootCounting.Combinatorics import binsum, choose
from PyPR.Tools.RootCounting.OverlappingRectangle import rectangle_solve
from PyPR.Tools.RootCounting.PartialOrders import maximalElements


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
            self.terms,
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
            basis = tuple(sorted(termset.totals.keys()))
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








    # for cube attacks:
    def get_cube_candidates(self):
        """Cube profiles for a cube attack, each with the degree its superpolys can reach.

        A candidate is a term of this profile with one variable removed from one
        block: a cube I with that per-block count. Summing the output over I leaves
        the superpoly, whose monomials come from the terms whose monomials can
        contain T_I -- those with at least the candidate's count in every block --
        with the candidate's variables removed. A term exceeding the candidate by
        e variables in total therefore contributes superpoly monomials of degree
        at most e, over the blocks where it exceeds it.

        So each candidate carries a degree bound, the largest excess over any term
        containing it, and target blocks, those where some containing term has
        variables to spare. Degree 1 is a linear superpoly. A candidate that only
        its own profile contains has a constant superpoly and says nothing about
        the state, so it is not returned; nor is a candidate with no variables,
        whose superpoly is the output itself.

        :return: ``(candidate, target_blocks, num_cubes, degree)`` for each
            candidate, where ``num_cubes`` counts the index sets with the
            candidate's per-block counts.
        :rtype: list[tuple[TermSet, tuple[int, ...], int, int]]
        """
        candidates = []
        already_added = set()
        for term_set in self.terms:
            for block_id in term_set.totals:
                # create the candidate and check if it's already been processed
                candidate = term_set.__copy__()
                candidate.counts[block_id] -= 1
                # a cube of no variables sums nothing: its "superpoly" is the
                # output itself, which is an algebraic attack's system, not a cube's
                if sum(candidate.counts.values()) == 0:
                    continue
                # keyed by block, not by the sorted sizes and counts: those
                # collide for distinct candidates such as <0:1/7, 1:2/5> and
                # <0:2/7, 1:1/5>, and the second was silently dropped
                already_added_key = tuple(sorted(
                    (i, candidate.totals[i], candidate.counts[i]) for i in candidate.totals
                ))

                if already_added_key in already_added:
                    continue
                already_added.add(already_added_key)

                # the superpoly's degree bound and the blocks its monomials use
                degree = 0
                targets = set()
                for other in self.terms:

                    # compare over every block either side mentions: a containing
                    # term with variables in a block the candidate doesn't list
                    # still puts those variables in the superpoly
                    diffs = {}
                    for i in set(candidate.counts) | set(other.counts):
                        other_count = other.counts[i] if i in other.counts else 0
                        candidate_count = candidate.counts[i] if i in candidate.counts else 0
                        diffs[i] = other_count - candidate_count

                    if any(x < 0 for x in diffs.values()):
                        continue # This set "sticks out" past the other term and is not a subset

                    # zero excess is the candidate's own profile: a constant
                    excess = sum(diffs.values())
                    if excess > 0:
                        degree = max(degree, excess)
                        for i, x in diffs.items():
                            if x > 0:
                                targets.add(i)

                if degree == 0:
                    continue

                num_cubes = 1
                for i in candidate.totals:
                    num_cubes *= choose(
                        candidate.totals[i],
                        candidate.counts[i]
                    )

                candidates.append((
                    candidate,
                    tuple(sorted(targets)),
                    num_cubes,
                    degree
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
            basis = tuple(sorted(termset.totals.keys()))
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
                while any(
                    (rollover_copy[i] > rects[rect_idx][i])
                    for i in range(len(rects[rect_idx]))
                ):
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

        # rollover loop which dynamically switches between the active term/rect
        rect_idx = 0
        curr_vec = [0 for i in range(dim)]

        # counteract the first incrementto start first yield with all zeroes
        # this allows the method to yield the constant term:
        curr_vec[0] -= 1

        while rect_idx < len(rects):
            # increment degree
            curr_vec[0] += 1

            # rollover loop:
            rollover_idx = 0
            rollover_copy = [x for x in curr_vec]

            # as long as any index is too large for the current rectangle:
            while any(
                (rollover_copy[i] > rects[rect_idx][i])
                for i in range(len(rects[rect_idx]))
            ):
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
