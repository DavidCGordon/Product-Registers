"""Synthesis of (nonlinear) feedback shift registers from an observed sequence.

Throughout this module a register of length m has two distinct maps attached to
it, and keeping them apart is what makes the index arithmetic below readable:

- the **feedback function** ``f: GF(2)^m -> GF(2)``, which produces one bit;
- the **update map** ``F: GF(2)^m -> GF(2)^m``, which produces the next state.

For a Fibonacci-configured register the second is built from the first:

    F(s)  =  (s[1], s[2], ..., s[m-1], f(s))

That is, the state shifts toward index 0 -- discarding s[0], which is the output
bit -- and f supplies the single new bit entering at index m-1.

**Sequences and states.** Seed the register with the first m bits of a sequence,
so that state s_t is the length-m frame starting at position t:

    s_t  =  (seq[t], seq[t+1], ..., seq[t+m-1])

Then F(s_t) = s_{t+1} holds for the shift coordinates automatically, and pins the
feedback on exactly one value: f(s_t) = seq[t+m]. So

    the register generates seq  <=>  f(seq[t], ..., seq[t+m-1]) = seq[t+m] for all t

and every claim below about "the sequence" is a claim about f evaluated on the
sliding frames of that sequence. Synthesis is the problem of finding the
shortest m for which some f satisfies that system, together with such an f.

See `docs/architecture/Bijective Register Synthesis.md` for the full development,
including the proof that the update map is a bijection exactly when the discarded
bit occurs in f as a lone linear term.
"""
from PyPR.BooleanLogic import XOR, AND, NOT, CONST, VAR


def KMP_table(seq):
    """The Knuth-Morris-Pratt prefix function of `seq`.

    Entry ``output[i]`` is the length of the longest proper prefix of
    ``seq[:i+1]`` that is also a suffix of it. `BM_NL` uses it on the *reversed*
    prefix of the observed sequence, where it answers a question about repeated
    windows -- see that function's docstring.

    :param seq: The sequence to build the table for.
    :type seq: list[int]
    :return: The prefix function, one entry per position of `seq`.
    :rtype: list[int]
    """
    output = [0]
    pref_len = 0
    for idx in range(1,len(seq)):
        while (pref_len > 0) and (seq[pref_len] != seq[idx]):
            pref_len = output[pref_len-1]
        
        if seq[pref_len] == seq[idx]:
            pref_len += 1

        output.append(pref_len)
    return output


def _BM_NL_bijective(seq):
    """Shortest Fibonacci register generating `seq` whose update map is a bijection.

    **When is the update map invertible?** Recovering s from F(s) is free in all
    but one coordinate: F(s) = (s[1], ..., s[m-1], f(s)) hands back s[1..m-1]
    directly, since they sit in positions 0..m-2 of the image. The only bit at
    risk is s[0], which the shift discards, and the only place it can still be
    recorded is the single new bit f(s).

    Over GF(2) each variable has degree at most 1, so the ANF of f splits
    uniquely on s[0]:

        f(s)  =  s[0] * A(s[1..m-1])  XOR  B(s[1..m-1])

    Write t = F(s). Since s[1..m-1] = t[0..m-2] are already known, A and B are
    known values, and the final coordinate reads t[m-1] = s[0]*A XOR B. When
    A = 1 this gives s[0] = t[m-1] XOR B, a unique preimage. When A = 0 the
    equation is t[m-1] = B, which does not mention s[0] at all: the states s and
    s XOR e_0 then share an image. Since A is a function of the remaining bits it
    can be 1 on some states and 0 on others, so F is injective exactly when
    A is identically 1 -- that is, when

        f(s)  =  s[0] XOR B(s[1..m-1])

    with s[0] occurring in no other monomial. A linear register satisfies this
    automatically whenever bit 0 is a tap, since every monomial is a singleton.

    **Turning that into a fitting problem.** Substituting the constrained form
    into the generation condition f(s_t) = seq[t+m] of the module docstring:

        seq[t+m]  =  seq[t] XOR B(seq[t+1], ..., seq[t+m-1])

    Both occurrences of seq on the right are known, so move them together:

        B(seq[t+1], ..., seq[t+m-1])  =  seq[t+m] XOR seq[t]
               window w[t], m-1 bits          target d[t]

    B is otherwise unconstrained -- it may be as nonlinear as it likes -- so the
    whole bijectivity requirement now sits in the shape of the equation rather
    than in B. What remains is a consistency question about a table: a suitable
    B exists for this m exactly when d[t] is a function of w[t], i.e. when no two
    positions share a window but demand different targets.

    Note w[t] is m-1 bits, not m: it is the state window with its low end
    s[0] = seq[t] removed, that term having moved to the left-hand side. The
    target seq[t+m] was never part of the state window -- it is the next bit.

    **Termination.** The search walks m upward and stops at the first consistent
    length. It always terminates for len(seq) >= 2: the loop runs to
    m = len(seq) - 1, where range(len(seq) - m) yields the single index t = 0, so
    the table holds one entry and cannot conflict.

    **Cost.** The bijective feedbacks are a subset of those the unconstrained
    search ranges over, so the length returned here is at least `BM_NL`'s, and
    in practice it is equal or slightly longer. The obstruction to a fit at a
    given m is a repeated (m-1)-window carrying disagreeing targets; there are
    about N^2/2 pairs of windows and 2^(m-1) windows to collide in, so repeated
    windows become scarce once 2^(m-1) outgrows the number of pairs -- which
    happens at roughly the length the unconstrained search needs anyway.

    :param seq: A prefix of a binary sequence.
    :type seq: list[int]
    :return: The register length m and a feedback function f with
        f(seq[t..t+m-1]) = seq[t+m] on every window, whose update map is a
        bijection.
    :rtype: tuple[int, BooleanFunction]
    :raises ValueError: If no length below len(seq) is consistent. For
        len(seq) >= 2 this cannot happen; for len(seq) == 1 the loop body never
        runs and the error is raised immediately.
    """
    for m in range(1, len(seq)):
        # d[t] must be determined by the window that survives the shift
        table = {}
        consistent = True
        for t in range(len(seq) - m):
            window = tuple(seq[t + 1:t + m])
            d = seq[t + m] ^ seq[t]

            if window in table and table[window] != d:
                # this window already demanded the other value, so no B can
                # satisfy both positions at this length
                consistent = False
                break
            table[window] = d

        if not consistent:
            continue

        # B = XOR of the minterms selecting the windows where d == 1. Window
        # position j is state bit j+1, since the window drops the discarded bit.
        minterms = []
        for window, d in table.items():
            if not d:
                continue
            if not window:
                minterms.append(CONST(1))
            else:
                minterms.append(AND(*[
                    VAR(j + 1) if bit else NOT(VAR(j + 1))
                    for j, bit in enumerate(window)
                ]))

        return m, XOR(VAR(0), *minterms)

    raise ValueError(f"No bijective register found for a sequence of length {len(seq)}")


def BM_NL(seq, bijective = False):
    """Recover a short Fibonacci register whose feedback generates `seq`.

    **The default search.** The register keeps a frame length m and a feedback
    function f accumulated as a XOR of minterms. At position n it predicts
    seq[n] by evaluating f on the current frame (seq[n-m], ..., seq[n-1]); on a
    mismatch it XORs in the minterm that is 1 at exactly that frame and 0
    everywhere else, flipping f's value at that one point and leaving the rest
    of the table intact.

    That correction is sound only if the frame determines the required output.
    If the same m-frame occurred earlier demanding a different next bit, flipping
    f at that frame repairs one position and breaks the other, and no amount of
    further correction converges. So m must be large enough that frames
    demanding different outputs are distinct, and the algorithm's other job is
    detecting when it is not and growing m.

    **What the KMP table measures.** Ambiguity means the current frame recurs
    earlier in the sequence, so the quantity needed is the length of the longest
    suffix of seq[:n] that also occurs ending at some earlier position. The
    prefix function supplies it after one reversal. In R = seq[:n] reversed, a
    prefix of length L is the most recent L bits of the sequence -- the current
    frame -- while a suffix of R[:p+1] of length L is an L-window ending earlier.
    The KMP prefix function of R compares exactly those, so

        s = max(KMP_table(seq[:n][::-1]))

    is that longest repeated suffix length. Worked example: for
    seq[:n] = [1,0,1,1,0,1] the table is [0,0,1,1,2,3], giving s = 3; the last
    three bits [1,0,1] also occur as seq[0..2], while the last four
    [1,1,0,1] occur nowhere earlier.

    That last observation is the general one: length s repeats and length s+1
    does not, so s+1 is the shortest frame separating the two occurrences. The
    code grows to exactly that when the current frame is among the ambiguous
    ones -- ``if s > m-1``, i.e. s >= m, set ``m = s+1`` -- and leaves m alone
    when s < m, where the frame is already unique.

    The counter k gates the measurement rather than changing it: it is reset on
    each growth and decremented every position, so the KMP check runs only once
    it falls below zero. This is the role the ``2L <= n`` test plays in linear
    Berlekamp-Massey -- both withhold a length increase until enough new symbols
    have arrived to justify one -- though the quantity being accumulated differs.

    **The bijective variant.** The register this returns reproduces the sequence
    but is generally *not* a bijection: nothing in the search constrains how the
    discarded bit enters f, and the recovered feedback carries it inside larger
    products. Such a register has no inverse (`Fibonacci.invert` rejects it), and
    its state graph has transients rather than being a disjoint union of cycles,
    so a seed can sit off the cycles and report a nonzero preperiod. The period
    itself is still well defined -- every orbit of a map on a finite set is
    eventually periodic; what fails is *pure* periodicity. Passing
    ``bijective=True`` restricts the search to the feedbacks that are invertible;
    see :func:`_BM_NL_bijective` for the characterization and the cost.

    :param seq: A prefix of a binary sequence.
    :type seq: list[int]
    :param bijective: If True, return the shortest register whose update map is
        a bijection, which may be longer. Defaults to False, the original
        behaviour.
    :type bijective: bool
    :return: The register length m and a feedback function f satisfying
        f(seq[t..t+m-1]) = seq[t+m] on every window of the observed sequence.
    :rtype: tuple[int, BooleanFunction]
    """
    if bijective:
        return _BM_NL_bijective(seq)

    # current shift / next jump in NLC
    k = 0
    # the nonlinear complexity / frame len
    m = 0

    #h is the current feedback function:
    h = XOR(CONST(seq[0]))

    #n is the current bit:
    for n in range(len(seq)):
        target = seq[n]

        frame = seq[n-m:n][::-1]
        predicted = h.eval(frame)

        #calculate the discrepancy:
        discrepancy = target ^ predicted

        #decrement k:
        k -= 1

        if discrepancy:
            
            #base case update
            if m == 0:
                k = n
                m = n

            # nonunique update/over half update??
            elif k < 0:
                
                # if the kmp length jumps, increase register size to accomodate
                s = max(KMP_table(seq[:n][::-1]))
                if (s > m-1):
                    k = s-(m-1)
                    m = s + 1

            # add new minterm (labels reversed so they don't have to be updated):
            variables = []
            for idx in range(m):
                if seq[n-1-idx]:
                    variables.append(VAR(idx))
                else:            
                    variables.append(NOT(VAR(idx)))
            h.add_arguments(AND(*variables))


    # relable the function inputs to be accurate
    reverse_labels = {idx: m-1-idx for idx in range(m)}
    return m, h.remap_indices(reverse_labels)


def BM_NL_iterator(seq, yield_rate = 1000, yield_corrected = True, bijective = False):
    """Stream the nonlinear fit, yielding the current result periodically.

    :param seq: A binary sequence or generator of one.
    :type seq: Iterable[int]
    :param yield_rate: Number of bits consumed between yields.
    :type yield_rate: int
    :param yield_corrected: If True, relabel the feedback's inputs to the
        state's own indexing before yielding. Ignored when bijective is True,
        whose feedback is built in that labelling already.
    :type yield_corrected: bool
    :param bijective: If True, yield the shortest fit whose update map is a
        bijection, matching `BM_NL(..., bijective=True)` on the prefix consumed
        so far.
    :type bijective: bool
    :return: Successive (length, feedback) pairs, the last covering all of seq.
    :rtype: Iterator[tuple[int, BooleanFunction]]
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
                yield _BM_NL_bijective(prefix)
        if len(prefix) % yield_rate:
            yield _BM_NL_bijective(prefix)
        return

    arr = []
    # current shift./jump in NLC
    k = 0
    # the complexity/frame len
    m = 0
    
    #seq can be a generator
    for n, bit in enumerate(seq):

        # initialize the feedback function:
        if n == 0:
            h = XOR(CONST(bit))
        
        arr.append(bit)
        target = bit


        if n % yield_rate == 0:
            if yield_corrected:
                reverse_labels = {idx: m-1-idx for idx in range(m)}
                yield m, h.remap_indices(reverse_labels)
            else:
                yield m, h

        
        #calculate the expected output (evaluate the function)
        frame = arr[n-m:n][::-1]
        predicted = h.eval(frame)

        # calculate the discrepancy:
        discrepancy = target ^ predicted
        # print(target,expected,d)

        # decrement k:
        k -= 1

        if discrepancy:

            # base case update
            if m == 0:
                k = n
                m = n

            # nonunique update/over half update??
            elif k < 0:
                
                #if the kmp length jumps, increase register size to accomodate
                s = max(KMP_table(arr[:n][::-1]))
                if (s > m-1):
                    k = s-(m-1)
                    m = s + 1

            # add a new minterm
            variables = []
            for idx in range(m):
                if arr[n-1-idx]:
                    variables.append(VAR(idx))
                else:            
                    variables.append(NOT(VAR(idx)))
            h.add_arguments(AND(*variables))

    # final yield
    reverse_labels = {idx: m-1-idx for idx in range(m)}
    yield m, h.remap_indices(reverse_labels)

