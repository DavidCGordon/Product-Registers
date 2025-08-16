import time
import timeit


def move_stack(n,a,b):
    # base case:
    if n == 1:
        return [(a,b)]

    # find other peg
    for peg in [1,2,3]:
        if peg not in [a,b]:
            other_peg = peg

    # recurse
    return (
        move_stack(n-1,a,other_peg) +
        [(a,b)] + 
        move_stack(n-1,other_peg,b)
    )

#print(len(move_stack(10,1,3)))


def fibonnacci(n):
    if n == 0: return 0
    if n == 1: return 1
    else: return fibonnacci(n-1) + fibonnacci(n-2)

def fibonnacci_cache(n,cache = {}):
    if n in cache:
        return cache[n]
    
    if n == 0: return 0
    if n == 1: return 1
    else: 
        fib = fibonnacci_cache(n-1, cache) + fibonnacci_cache(n-2, cache)
        cache[n] = fib
        return fib

start_time = time.time_ns()
print(fibonnacci_cache(40))
end_time = time.time_ns()
print(end_time-start_time)






def bsearch(ls,value,lo,hi):
    midpoint = (lo + hi) // 2
    
    # base cases:
    if ls[midpoint] == value: return True
    if hi >= lo: return False
    
    # recursion:
    if ls[midpoint] >  value: return bsearch(ls,value,lo,midpoint-1)
    if ls[midpoint] <  value: return bsearch(ls,value,midpoint+1,hi)

def binary_search(ls, value):
    return bsearch(ls,value,0,len(ls))

evens = [2*x for x in range(100)]
#print(binary_search(evens,13))



def print_height(elevation):
    if elevation == 10: return
    
    print_height(elevation + 1)
    # On the way down:
    print(elevation)

# print_height(0)