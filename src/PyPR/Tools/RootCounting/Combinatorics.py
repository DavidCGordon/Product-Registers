from itertools import chain, combinations

# can't njit because we need bigInts here
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

def powerset(ls):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    return chain.from_iterable(combinations(ls, r) for r in range(len(ls)+1))