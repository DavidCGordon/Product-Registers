from math import ceil, copysign

# May be very optimiseable, but nicely leans on pythons integer
# implementations right now - 
# uses python bigInts to represent p_adic integer operations
# http://cs.engr.uky.edu/~klapper/pdf/fcsr.pdf primarily
# some modifications/ideas taken from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1056960

# http://ijns.jalaxy.com.tw/contents/ijns-v21-n1/ijns-2019-v21-n1-p1-6.pdf
# - 2-adic complexity of T-Functions (we can beat them)

from math import log2

def BM_FCSR(seq):
    k = 0
    while k < len(seq) and seq[k] == 0:
        k += 1

    # handle all zero case
    if k == len(seq):
        return 1,0,1
    
    m = k
    seq_approx = 2**k

    num_prev = 0
    den_prev = 1

    num_curr = 2**k
    den_curr = 1

    for i in range(k+1,len(seq)):
        seq_approx += int(seq[i])*(2**i)
        
        #discrepancy check
        if (seq_approx * den_curr - num_curr) % (2**(i+1)) != 0:
            scale = 2**(i-m)
            if phi(num_curr,den_curr) < scale * phi(num_prev,den_prev):
                temp_num = num_curr
                temp_den = den_curr

                #scale the current guess
                d = D(scale*num_prev, scale*den_prev, num_curr, den_curr)
                num_curr = d * num_curr + scale * num_prev
                den_curr = d * den_curr + scale * den_prev

                #update previous condition
                num_prev = temp_num
                den_prev = temp_den
                m = i
            else:
                #scale previous guess:
                d = D(num_curr, den_curr, scale*num_prev, scale*den_prev)
                num_curr = num_curr + d * 2**(i-m) * num_prev
                den_curr = den_curr + d * 2**(i-m) * den_prev

    # correct signs
    if den_curr < 0:
        den_curr *= -1
        num_curr *= -1

    return FCSR_size(num_curr,den_curr), num_curr, den_curr




#takes an iterator and returns an iterator
def BM_FCSR_iterator(seq, yield_rate):
    seq = enumerate(seq)

    k = 0
    while next(seq)[1] == 0:
        k += 1

    m = k
    seq_approx = 2**k

    num_prev = 0
    den_prev = 1

    num_curr = 2**k
    den_curr = 1

    for i, bit in seq:
        
        #if it's time, yield:
        if i % yield_rate == 0:
            if den_curr < 0:
                yield (
                    FCSR_size(num_curr * -1,den_curr * -1), 
                    num_curr * -1, 
                    den_curr * -1
                )
            else:
                yield (
                    FCSR_size(num_curr,den_curr), 
                    num_curr,
                    den_curr
                )


        seq_approx += int(bit)*(2**i)
        
        #discrepancy check
        if (seq_approx * den_curr - num_curr) % (2**(i+1)) != 0:
            scale = 2**(i-m)
            if phi(num_curr,den_curr) < scale * phi(num_prev,den_prev):
                temp_num = num_curr
                temp_den = den_curr

                #scale the current guess
                d = D(scale*num_prev, scale*den_prev, num_curr, den_curr)
                num_curr = d * num_curr + scale * num_prev
                den_curr = d * den_curr + scale * den_prev

                #update previous condition
                num_prev = temp_num
                den_prev = temp_den
                m = i
            else:
                #scale previous guess:
                d = D(num_curr, den_curr, scale*num_prev, scale*den_prev)
                num_curr = num_curr + d * 2**(i-m) * num_prev
                den_curr = den_curr + d * 2**(i-m) * den_prev

    #pseudo-return:
    # fix numerator and denominator signs:
    if den_curr < 0:
        yield (
            FCSR_size(num_curr * -1,den_curr * -1), 
            num_curr * -1, 
            den_curr * -1
        )
    else:
        yield (
            FCSR_size(num_curr,den_curr), 
            num_curr,
            den_curr
        )


"""
HELPER METHODS
"""
# Find the appropriate scale factor
# b1/b2 = g is the one that gets scaled.

# Note that the original paper has this implemented incorrectly
# use https://core.ac.uk/download/pdf/82034673.pdf instead

def D(a1,a2,b1,b2):
    options = []
    if b1 != b2:
        options += list(odd_round(-(a1-a2)/(b1-b2)))
    if b1 != -b2:
        options += list(odd_round(-(a1+a2)/(b1+b2)))
    d = min(options, key = lambda x: phi(x*b1+a1,x*b2+a2))
    return d

def odd_round(x):
    v = ceil(x) // 2 * 2 + 1
    return (v, v-2)

def phi(num,den): 
    return max(abs(num),abs(den))

# Determine the FCSR size for this numerator and denominator
def FCSR_size(num,den):
    if  num > 0:
        size = 1 + ceil(log2(max(abs(num),abs(den))))
    else:
        size = max(1 + ceil(log2(abs(num)/3 + 1)), ceil(log2(den)))
    return size


# alternate sequence generation to debug against
def fcsr_eval(n,d,lim = None):
    i = 0
    while (not lim) or (i < lim):
        i += 1
        a = (n % 2)
        n = (n-d*a) // 2
        yield a