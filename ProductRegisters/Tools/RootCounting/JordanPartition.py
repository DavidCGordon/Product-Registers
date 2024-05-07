from math import comb
# What a crazy algorithm: https://arxiv.org/pdf/math/0612437.pdf
# note this should be symmetric.

# counts multiplicity of the prime instead of computing the combinations
# this is slower at small scales, but scales better and is more accurate/stable
def Dp(i,s,t,p):
    if i >= s: return False

    count = 0

    for term in range(i+1):
        for decrement in range(t-1-i):
            num_term = (s+t-2-2*i) + (term - decrement)
            den_term = (t-1-i) + (term - decrement)

            while num_term % p == 0:
                num_term //= p
                count += 1

            while den_term % p == 0:
                den_term //= p
                count -= 1
    return (count > 0)

# used to iterate through the Dp Array, as specified in paper.
def update(i,s,t,p):
    f = 0
    while (Dp(i + f,s,t,p) and (i + f) < s):
        f += 1
    
    length = s+t-1-2*i - f
    degree = f+1
    next_idx = i+f+1

    return length, degree, next_idx

# final method
def JP_solve(s,t,p):
    s,t = sorted([s,t])
    out = []
    idx = 0
    while idx < s:
        length, degree, idx = update(idx,s,t,p)
        out.append((length,degree))
    return out

