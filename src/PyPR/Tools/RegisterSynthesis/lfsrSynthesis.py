from numba import njit
import numpy as np

def _berlekamp_massey_bijective(seq):
    """Shortest LFSR generating `seq` whose update map is a bijection.

    A Fibonacci register is a bijection exactly when the bit the shift discards
    enters the feedback linearly and alone, i.e. f(s) = s[0] XOR B(s[1..m-1]).
    Imposing that on the generation condition f(seq[t..t+m-1]) = seq[t+m] and
    moving the known terms together gives

        B(seq[t+1], ..., seq[t+m-1])  =  seq[t+m] XOR seq[t]
               window w_t, m-1 bits           target d_t

    For a *linear* register B must itself be linear, so this is a system of
    GF(2) equations in the m-1 unknown tap coefficients rather than the lookup
    table the nonlinear variant solves. A fit exists at length m exactly when
    that system is consistent, which elimination decides.

    .. warning::
        A bijective register emits a *purely* periodic sequence, because every
        state of a bijection lies on a cycle. If `seq` is not purely periodic -- if it has a genuine transient -- no
        bijective register reproduces it, and the search degenerates rather than
        failing. Writing rho for the least period and tau for the minimal
        preperiod, positions tau-1 and tau-1+rho carry equal windows but opposite
        targets, so every length up to N - tau - rho is inconsistent and the
        result is N - tau - rho + 1: slope 1 in len(seq). See `docs/architecture/Bijective Register Synthesis.md`.

    :param seq: A prefix of a binary sequence.
    :type seq: list[int]
    :return: The register length m and the degree-m connection polynomial, in
        the same coefficient convention `berlekamp_massey` returns.
    :rtype: tuple[int, numpy.ndarray]
    :raises ValueError: If no length below len(seq) admits a bijective fit.
    """
    seq = [int(b) for b in seq]

    for m in range(1, len(seq)):
        rows = [[seq[t + 1 + j] for j in range(m - 1)] for t in range(len(seq) - m)]
        targets = [seq[t + m] ^ seq[t] for t in range(len(seq) - m)]

        # Gaussian elimination over GF(2) on [rows | targets]
        width = m - 1
        augmented = [row[:] + [targets[i]] for i, row in enumerate(rows)]
        pivots, pivot_row = [], 0
        for col in range(width):
            candidate = next(
                (i for i in range(pivot_row, len(augmented)) if augmented[i][col]), None
            )
            if candidate is None:
                continue
            augmented[pivot_row], augmented[candidate] = \
                augmented[candidate], augmented[pivot_row]
            for i in range(len(augmented)):
                if i != pivot_row and augmented[i][col]:
                    augmented[i] = [x ^ y for x, y in zip(augmented[i], augmented[pivot_row])]
            pivots.append(col)
            pivot_row += 1

        # inconsistent iff some row reads 0 ... 0 | 1
        if any(all(v == 0 for v in row[:width]) and row[width] for row in augmented):
            continue

        taps = [0] * width
        for r, col in enumerate(pivots):
            taps[col] = augmented[r][width]

        # Coefficient convention: _from_poly taps bit (size - idx) for each
        # nonzero P[idx] with idx >= 1, and drops the lowest-index entry, which
        # is why P[0] is set. Bit 0 needs idx = m; window position j is state
        # bit j+1, hence idx = m-1-j.
        polynomial = np.zeros(m + 1, dtype='uint8')
        polynomial[0] = 1
        polynomial[m] = 1
        for j, bit in enumerate(taps):
            if bit:
                polynomial[m - 1 - j] = 1

        return m, polynomial

    raise ValueError(
        f"No bijective LFSR of length below {len(seq)} generates this sequence."
    )


def berlekamp_massey(seq, bijective = False):
    """Recover the shortest LFSR whose output matches `seq`.

    The polynomial returned satisfies the **dual** relationship with the
    sequence -- it convolves it to zero -- which is the convention `Fibonacci`
    and `Galois` take as constructor input, so no reversal is needed at the
    boundary. See `docs/conventions/Polynomial Conventions.md`.

    The register this produces is **not** always a bijection. It taps bit 0, and
    hence is invertible, only when the returned polynomial reaches full degree L.
    Nothing in the algorithm forces that tap, so a singular update map is a
    routine outcome on unstructured input rather than an edge case worth
    ignoring. Passing `bijective=True`
    restricts the search to the invertible ones -- but read the warning in
    :func:`_berlekamp_massey_bijective` first, because on a sequence with a
    transient that search degenerates rather than failing outright.

    :param seq: A prefix of a binary sequence.
    :type seq: list[int] | numpy.ndarray
    :param bijective: If True, return the shortest LFSR whose update map is a
        bijection. Defaults to False, the original behaviour.
    :type bijective: bool
    :return: The linear complexity and the connection polynomial.
    :rtype: tuple[int, numpy.ndarray]
    """
    if bijective:
        return _berlekamp_massey_bijective(seq)

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

def berlekamp_massey_iterator(seq, yield_rate = 1000, bijective = False):
    """Stream the Berlekamp-Massey fit, yielding the current result periodically.

    :param seq: A binary sequence or generator of one.
    :type seq: Iterable[int]
    :param yield_rate: Number of bits consumed between yields.
    :type yield_rate: int
    :param bijective: If True, yield the shortest fit whose update map is a
        bijection, matching `berlekamp_massey(..., bijective=True)` on the
        prefix consumed so far.
    :type bijective: bool
    :return: Successive (length, polynomial) pairs, the last covering all of seq.
    :rtype: Iterator[tuple[int, numpy.ndarray]]
    :raises ValueError: If bijective is True and some prefix admits no
        bijective fit.
    """
    if bijective:
        # The bijective fit is a rescan over m, not an incremental update, so
        # there is no streaming core to reuse: each yield reports the fit for
        # the prefix consumed so far. Cost per yield is that of the batch call.
        prefix = []
        for bit in seq:
            prefix.append(bit)
            if len(prefix) % yield_rate == 0:
                yield _berlekamp_massey_bijective(prefix)
        if len(prefix) % yield_rate:
            yield _berlekamp_massey_bijective(prefix)
        return

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