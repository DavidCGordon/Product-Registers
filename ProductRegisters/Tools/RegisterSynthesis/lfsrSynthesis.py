def berlekamp_massey(seq):
    #N = total number of bits to process
    N = len(seq)

    # current connection polynomial guess
    curr_guess = [1] + [0 for i in range(N-1)]
    # prev. connection polynomial guess
    prev_guess = [1] + [0 for i in range(N-1)]

    # L = current linear complexity
    L = 0
    # m = index of last change
    m = -1
    
    #n = index of bit we are correcting.
    for n in range(N):

        # calculate discrepancy from LFSR frame
        d = 0
        for i in range(L+1):
            d ^= (curr_guess[i] & seq[n-i])

        #handle discrepancy (if needed)
        if d != 0:
            
            #store copy of current guess
            temp = curr_guess[:]

            #curr_guess = curr_guess - (x**(n-m) * prev_guess)
            shift = n-m
            for i in range (shift, N):
                curr_guess[i] ^= prev_guess[i - shift]

            #if 2L <= n, then the polynomial is unique
            #it's safe to update the linear complexity.
            if 2*L <= n:
                L = n + 1 - L
                prev_guess = temp
                m = n

    #return the linear complexity and connection polynomial
    return (L, curr_guess[:L+1])

def berlekamp_massey_iterator(seq, yield_rate = 1000):
    
    arr = []
    curr_guess = [1]
    prev_guess = [1]

    # L = current linear complexity
    linear_complexity = 0
    # m = index of last change
    m = -1

    for n, bit in enumerate(seq):

        arr.append(bit)
        curr_guess.append(0)
        prev_guess.append(0)

        if n % yield_rate == 0:
            yield (linear_complexity, curr_guess[:linear_complexity + 1])

        #calculate discrepancy from LFSR frame
        discrepancy = 0
        for i in range(linear_complexity + 1):
            discrepancy ^= (curr_guess[i] & arr[n-i])

        #handle discrepancy (if needed)
        if discrepancy:
            
            #store copy of current guess
            temp = curr_guess[:]

            #update current guess
            shift = n-m
            for i in range (shift, n+1):
                curr_guess[i] ^= prev_guess[i - shift]

            #update LC 
            if 2 * linear_complexity <= n:
                linear_complexity = (n + 1) - linear_complexity
                prev_guess = temp
                m = n

    while True:
        yield(linear_complexity, curr_guess[:linear_complexity + 1])