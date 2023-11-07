from numba import njit
import numpy as np

def berlekamp_massey(seq):
    N = len(seq)
    if type(seq) != np.ndarray:
        seq = np.asarray(seq, dtype='uint8')
    return _berlekamp_massey(N,seq)

@njit
def _berlekamp_massey(N,seq):
    # N = total number of bits to process
    # current connection polynomial guess
    curr_guess = np.zeros(N, dtype='uint8')
    curr_guess[0] = 1
    # prev. connection polynomial guess
    prev_guess = np.zeros(N, dtype='uint8')
    prev_guess[0] = 1

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
            temp = curr_guess.copy()

            #curr_guess = curr_guess - (x**(n-m) * prev_guess)
            shift = n-m
            for i in range(shift, N):
                curr_guess[i] ^= prev_guess[i - shift]

            #if 2L <= n, then the polynomial is unique
            #it's safe to update the linear complexity.
            if 2*L <= n:
                L = n + 1 - L
                prev_guess = temp
                m = n

    #return the linear complexity and connection polynomial
    return (L, curr_guess[:L+1])




    

@njit
def _bm_iterator_core(
    start_idx,yield_rate,
    arr,curr_guess,prev_guess,
    linear_complexity,last_update):

    # if it's time to resize (powers of 2)
    for n in range(start_idx, start_idx + yield_rate):

        #calculate discrepancy from LFSR frame
        discrepancy = 0
        for i in range(linear_complexity + 1):
            discrepancy ^= (curr_guess[i] & arr[n-i])

        #handle discrepancy (if needed)
        if discrepancy:

            #store copy of current guess
            temp = curr_guess.copy()

            #update current guess
            shift = n-last_update
            for i in range (shift, n+1):
                curr_guess[i] ^= prev_guess[i - shift]

            #update LC 
            if 2 * linear_complexity <= n:
                linear_complexity = (n + 1) - linear_complexity
                prev_guess = temp
                last_update = n
    return arr, curr_guess, prev_guess, linear_complexity, last_update

def berlekamp_massey_iterator(seq, yield_rate = 1000):
    arr_size = 2**10
    arr = np.zeros(arr_size, dtype='uint8')

    curr_guess = np.zeros(arr_size, dtype='uint8')
    curr_guess[0] = 1

    prev_guess = np.zeros(arr_size, dtype='uint8')
    prev_guess[0] = 1

    linear_complexity = 0
    last_update = -1

    NotEnded = True
    start_idx = 0

    while NotEnded:
        # if it's time to resize (powers of 2)
        while start_idx + yield_rate >= arr_size:
            new_arr = np.zeros(arr_size * 2, dtype='uint8')
            new_arr[:arr_size] = arr
            arr = new_arr

            new_curr_guess = np.zeros(arr_size * 2, dtype='uint8')
            new_curr_guess[:arr_size] = curr_guess
            curr_guess = new_curr_guess

            new_prev_guess = np.zeros(arr_size * 2, dtype='uint8')
            new_prev_guess[:arr_size] = prev_guess
            prev_guess = new_prev_guess

            arr_size *= 2

        # grow arr:
        new_chunk = [x for _,x in zip(range(yield_rate), seq)]
        arr[start_idx : start_idx+len(new_chunk)] = new_chunk
        if len(new_chunk) < yield_rate:
            NotEnded = False
        

        # update variables with JIT code
        arr,curr_guess,prev_guess,linear_complexity,last_update \
            = _bm_iterator_core(
                start_idx, len(new_chunk),
                arr,curr_guess,prev_guess,
                linear_complexity,last_update
            )

        # update the index and yield
        start_idx += yield_rate
        yield(linear_complexity, curr_guess[:linear_complexity + 1])
    raise StopIteration