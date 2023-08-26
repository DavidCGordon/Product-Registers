from ProductRegisters.BooleanLogic import XOR, AND, NOT, CONST, VAR


# build a prefix table using KMP algorithm:
def KMP_table(seq):
    output = [0]
    pref_len = 0
    for idx in range(1,len(seq)):
        while (pref_len > 0) and (seq[pref_len] != seq[idx]):
            pref_len = output[pref_len-1]
        
        if seq[pref_len] == seq[idx]:
            pref_len += 1

        output.append(pref_len)
    return output


def BM_NL(seq):
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


def BM_NL_iterator(seq, yield_rate = 1000, yield_corrected = True):
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
    raise StopIteration
