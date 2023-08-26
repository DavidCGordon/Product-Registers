# leq_ordering(A,B) returns true if A <= B
# input sets is an iterable of iterables.
def maximalElements(leq_ordering, inputs):
    maximal_set = set()
    for container in inputs:
        for element in container:
            maximal = True

            for already_added in maximal_set:
                if leq_ordering(already_added, element):
                    maximal_set.remove(already_added)
                    break

                elif leq_ordering(element, already_added):
                    maximal = False
                    break
                
            if maximal:
                maximal_set.add(element)
    return maximal_set

# In general:
# Exact => Same multiplicities
# Subset => All non-full bases have count_a < count_b
# Embedded => Allow for additional full bases in the larger term.


# If the roots in A are a subset of those in B, with the same multiplicity, then A is redundant
# For general cleaning, we can't necessarily remove different multiplicies due needing them for future jordan calculations
# We remove these frequently to declutter root expressions
def isExactSubset(a,b):
    if a.m != b.m:
        return False

    if a.roots.keys() != b.roots.keys():
        return False

    # all counts in A must be < B to be a subset.
    for k in a.roots.keys():
        if a.roots[k] > b.roots[k]:
            return False

    # if all conditions pass, return true
    return True

# If the roots in A are a subset of those in B, with same or lower multiplicity
# We can finally remove these at evaluation to avoid double counting.
def isSubset(a,b):
    if a.m > b.m:
        return False

    if a.roots.keys() != b.roots.keys():
        return False

    # all counts in A must be < B to be a subset.
    for k in a.roots.keys():
        if a.roots[k] > b.roots[k]:
            return False

    # if all conditions pass, return true
    return True

def isSubspace(a,b):
    a_keyset = set(a.roots.keys())
    b_keyset = set(b.roots.keys())

    # for sets in python this operator indicates a strict subset
    return a_keyset < b_keyset










# If a copy of A might be embedded in B, with the same multiplicity
# This is handled in the evaluation method (e.g. <7,7> = 126 or 127)
# We remove these frequently, to declutter the root expressions
def isEmbeddedSet(a,b):
    if a.m != b.m:
        return False

    a_keyset = set(a.roots.keys())
    b_keyset = set(b.roots.keys())
    
    # there should be a strict subset relationship on keys:
    if not a_keyset.issubset(b_keyset):
        return False

    # all shared keys should have equal values.
    for k in b_keyset.intersection(a_keyset):
        if a.roots[k] != b.roots[k]:
            return False

    # all non shared keys should be maximum (for each pair (b,c), b should equal c).
    for k in (b_keyset - a_keyset):
        if b.roots[k] != k:
            return False

    # if all conditions pass, return true
    return True

# Formally: If set A is a subset of something that might be embedded in B
# Alternatively: If all the roots are inside the area optimistically covered by B.
# we remove these when we optimistically assume that we will get everything in B
# if we aren't optimistic, we might still want the smaller one. 

# unclear if we can ever actually get these?
def isEmbeddedSubset(a,b):
    if a.m > b.m:
        return False
    
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

