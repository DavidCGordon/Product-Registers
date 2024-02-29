from ProductRegisters.Tools.RootCounting.Combinatorics import choose
from ProductRegisters.Tools.RootCounting.PartialOrders import maximalElements
from ProductRegisters.BooleanLogic import AND, XOR, CONST, VAR
from itertools import product




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



def isSubset(a,b):
    # all counts in A must be < B to be a subset.
    for block_id in a.totals:
        compare_value = b.counts[block_id] if block_id in b.counts else 0
        if a.counts[block_id] > compare_value:
            return False

    # if all conditions pass, return true
    return True


class MonomialProfile:
    # wrapper for JordanSets, with extra functionality:
    # self.terms is a set of JordanSets.
    def __init__(self,term_list = None):
        if term_list is None:
            self.terms = set()
        else:
            self.terms = set(term_list)

    @classmethod             
    def from_merged(self, fn_list, blocks):
        bitmap = [None for i in range(sum(len(block) for block in blocks))]
        for block_id in range(len(blocks)):
            for bit in blocks[block_id]:
                bitmap[bit] = MonomialProfile([TermSet(
                    totals={block_id:len(blocks[block_id])},
                    counts={block_id:1}
                )])

        total_fn = XOR(*fn_list)
        total_fn = total_fn.remap_constants({
            0: MonomialProfile.logical_zero(),
            1: MonomialProfile.logical_one()
        })

        return total_fn.eval_ANF(bitmap)
    

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
            leq_ordering=isSubset, 
            inputs=[self.terms, other.terms]
        )

        return MonomialProfile(new_terms)

    def __and__(self, other): return self.__mul__(other)
    def __mul__(self, other):

        new_terms = [a*b for a, b in product(self.terms, other.terms)]

        #clean out redundant subsets.
        new_terms = maximalElements(
            leq_ordering = isSubset,
            inputs = [new_terms]
        )

        return MonomialProfile(new_terms)

    # When multiplying by Logical One, you should leave the result untouched
    # When adding with Logical Zero, you should add an indicator term 
    # These effects are accomplished by the MonomialProfile with an Empty TermSet
    @classmethod
    def logical_one(self): return MonomialProfile([TermSet({},{})])

    # When multiplying by Logical Zero, you should cancel out terms
    # When adding with Logical zero, you should leave the result untouched
    # These effects are accomplished by the empty MonomialProfile
    @classmethod
    def logical_zero(self): return MonomialProfile([])

    # The Monomial Profile of an inverted function is the same
    def __invert__(self): return self ^ MonomialProfile.logical_one()


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

                    if isSubset(modified_set,other):
                        useful = False
                        break

                if useful:
                    num_cubes = 1
                    for i in modified_set.totals:
                        num_cubes *= choose(
                            modified_set.totals[i],
                            modified_set.counts[i]
                        )

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



    
                

    