from itertools import islice

def periodTest(test, show = False):
    matches = []
    for idx, val in enumerate(islice(test,2**(test.size)+1)):
        if show: print(idx, val)
        if idx == 0:
            start = val
        elif val == start:
            matches.append(idx)
    return matches

def isLinear(bitFn, allowAfine = False):
    valid = True
    for ls in bitFn.fn:
        for l in ls:
            if type(l) == bool:
                valid &= allowAfine
            elif len(l) != 1:
                valid = False
    return valid

def berlekamp_massey(bits):
    N = len(bits)
    c = [1] + [0 for _ in range(N-1)]
    c[0] = 1
    b = [1] + [0 for _ in range(N-1)]
    b[0] = 1

    L = 0
    m = -1
    n = 0

    while n < N:
        d = bits[n]
        for i in range(1,L+1):
            d ^= (c[i] & bits[n-i])
        
        if d != 0:
            t = c[:]
            
            for i in range (n-m, N):
                c[i] ^= b[i - n + m]

            if L <= (n/2):
                L = n + 1 - L
                m = n
                b = t
        n += 1
    return c[:L]

def getOffset(bitSeq1, bitSeq2):
    N = len(bitSeq1)
    n = 0 #the index being tested
    i = 0 #the number correct in a row
    while i < N:
        if n > N+1: return None #prevents infinite loop in case of no match

        if bitSeq1[(n+i)%N] != bitSeq2[(i)%N]:
            n += 1
            i = 0
        else:
            i += 1
    return n