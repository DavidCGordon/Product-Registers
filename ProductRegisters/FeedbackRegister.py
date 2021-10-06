from BitVector import BitVector

class FeedbackRegister:
    #INITIALIATION/DATA:
    def __init__(self, seed, fn):
        
        #feedback variables
        self.fn = fn
        self.size = len(fn)

        #state and seed (internal use only)
        if type(seed) == int:
            self._state = BitVector(intVal = seed, size = self.size)
        else:
            self._state = seed            
        self._seed = self._state



    #TYPE CONVERSIONS
    def __int__(self): return int(self._state)
    def __str__(self): return str(self._state)
    def __list__(self): return list(self._state)[::-1]



    #STATE MANIPULATION:
    def _bit(self, bitIdx): return (self.size-1-bitIdx % self.size)
    def __getitem__(self, key): return list(self._state)[::-1][key]
    def __setitem__(self, idx, val): self._state[self.bit(idx)] = val
    def __reverse__(self): self._state.reverse()



    #ITERATION:
    
    #iterate through bits in the register
    def __iter__(self):
        for i in range(self.size):
            yield self._state[self._bit(i)]


    #iterate register state through time:
            
    #clock the register
    def clock(self):
        nextState = BitVector(intVal = 0, size = self.size)
        
        for bitIdx in range(self.size):
            summation = 0
            for term in self.fn[bitIdx]:
                if type(term) == bool:
                    product = int(term)
                else:
                    product = 1
                    for idx in term:
                        product &= self._state[self._bit(idx)]
                summation ^=  product
            nextState[(self._bit(bitIdx))] = summation
        self._state = nextState
        
    #clock the register and XOR input:
    def input(self,inpt):
        l = len(inpt)
        if type(inpt) == list:
            inpt = BitVector(bitlist = inpt)
        inpt.pad_from_right(self.size-l)
        self.clock()
        self._state ^= inpt

    #reset to initial_state
    def reset(self):
        self._state = self._seed
        
    #generate a sequence of states, in order
    def run(self, arg = None):
        #different input cases:

        #list of inputs
        if type(arg) == list:
            yield self
            for inpt in arg:
                self.input(inpt)
                yield self

        #number of iterations to run
        elif type(arg) == int:
            for i in range(arg-1):
                yield self
                self.clock()
            yield self
                
        #no limit
        elif arg == None:
            while True:
                yield self
                self.clock()
                

    #DIAGNOSTIC AND EXTRA INFO     
    #return the period of the register:
    def period(self, lim = 2**18):
        curState = self._state
        first_state = self._state
        self.clock()
        count = 1
        while self._state != first_state:
            self.clock()
            count += 1
            if count > lim:
                return None
        self._state = curState
        return count

    #linear complexity profile of the sequence:
    def LinearComplexity(self,bit,lim):
        #generate a sequence of length lim
        seq = []
        curState = self._state
        for state in self.run(lim):
            seq.append(state[bit])
        self._state = curState

        Ls = []
        #berlekamp-massey
        N = len(seq)
        #c = current connection polynomial
        c = [1] + [0 for i in range(N-1)]
        #b = prev. connection polynomial
        b = [1] + [0 for i in range(N-1)]

        #L = len(LFSR frame) = upper bound on deg(C)
        L = 0
        #m is index of last change
        m = -1
        
        #n = bit we are correcting.
        for n in range(N):
            #calculate discrepancy from LFSR frame
            d = 0
            for i in range(L+1):
                d ^= (c[i] & seq[n-i])

            #handle discrepancy if needed
            if d != 0:
                
                #store copy of C
                temp = c[:]

                #c = c - x^(n-m)b
                shift = n-m
                for i in range (shift, N):
                    c[i] ^= b[i - shift]

                #if 2L <= n, then the polynomial is unique
                #it's safe to update
                if 2*L <= n:
                    L = n + 1 - L
                    Ls.append((n,L))
                    b = temp
                    m = n
       
        return Ls
