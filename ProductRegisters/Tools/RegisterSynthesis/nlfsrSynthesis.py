def BM_NL(seq):
    #current shift./jump in NLC
    k = 0

    #the complexity/frame len
    m = 0

    #h is the current feedback:
    h = [seq[0]]

    #n is the current bit:
    for n in range(len(seq)):
        target = seq[n]

        #calculate the expected output (evaluate the function)
        predicted = 0
        for term in h:
            #terms
            if type(term) == list:
                prod = 1
                for idx,complement in term:
                    prod &= (seq[n - 1 - idx] ^ complement)
                predicted ^= prod
            #constants
            elif type(term) == int:
                predicted ^= term

        #calculate the discrepancy:
        d = target ^ predicted
        #print(target,expected,d)

        #decrement k:
        k -= 1

        if d:
            
            #base case update
            if m == 0:
                k = n
                m = n

            #nonunique update/over half update??
            elif k < 0:
                
                #if the kmp length jumps, increase register size to accomodate
                s = max(KMP_table(seq[:n][::-1]))
                if (s > m-1):
                    k = s-(m-1)
                    m = s + 1

            #construct minterm f:
            f = [(idx, 1-seq[n-1-idx]) for idx in range(m)]
            #add minterm f to h
            h.append(f)
        
        #print(k,m)

    #relable the function inputs to fit with convention.
    reverse = []
    for term in h:
        if type(term) == list:
            reverse.append([(m-1-idx,val) for idx,val in term])
        elif term == 1:
            reverse.append(1)

    return m, reverse

def ev(seq, i):
    return i - max(KMP_table(seq[:i][::-1]))

#build a prefix table using kmp algorithm:
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

def eval(seq,m,fn):
    for n in range(m,len(seq)):
        print(seq[n])

        intermeds = []

        expected = 0
        for term in fn:
            #terms
            if type(term) == list:
                prod = 1
                for idx,complement in term:
                    prod &= (seq[n - m + idx] ^ complement)
                
                intermeds.append(prod)
                expected ^= prod
            #constants
            elif type(term) == int:
                expected ^= term
        #print(intermeds)
        print(expected,"\n")
