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

def BM_iterative(reg, bit_idx, block_size):

    seq = []
    N = 0

    #generate initial segment of sequence
    for state in reg.run(block_size):
        seq.append(state[bit_idx])
    N += block_size

    curr_guess = [1] + [0 for i in range(N-1)]
    prev_guess = [1] + [0 for i in range(N-1)]

    # L = current linear complexity
    L = 0
    # m = index of last change
    m = -1
    
    # n = index of bit we are correcting.
    n = -1
    while True:
        n += 1

        #if we have finished a block:
        if n == N:
            yield (L)
            for state in reg.run(block_size):
                seq.append(state[bit_idx])
                curr_guess.append(0)
                prev_guess.append(0)
            N += block_size

        #calculate discrepancy from LFSR frame
        d = 0
        for i in range(L+1):
            d ^= (curr_guess[i] & seq[n-i])

        #handle discrepancy (if needed)
        if d != 0:
            
            #store copy of current guess
            temp = curr_guess[:]

            #update current guess
            shift = n-m
            for i in range (shift, N):
                curr_guess[i] ^= prev_guess[i - shift]

            #update LC 
            if 2*L <= n:
                L = n + 1 - L
                prev_guess = temp
                m = n
