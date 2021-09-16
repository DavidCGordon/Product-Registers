from BitVector import BitVector

class FeedbackRegister:
    #Initialization and data:
    def __init__(self, seed, fn):
        
        #register variables:
        self.order = "MSB"
        
        #feedback variables
        self.fn = fn
        self.size = len(fn)

        #state and seed (internal use only)
        if type(seed) == int:
            self._state = BitVector(intVal = seed, size = self.size)
        else:
            self._state = seed            
        self._seed = self._state



    #Type conversions
    def __int__(self): return int(self._state)
    def __str__(self): return str(self._state)
    def __list__(self): return list(self._state)[::-1]


    
    #State manipulation methods:
    def bit(self, bitIdx): return (self.size-1-bitIdx % self.size)
    def __getitem__(self, key): return list(self._state)[::-1][key]
    def __setitem__(self, idx, val): self._state[self.bit(idx)] = val



    #ITERATION METHODS:
    #iterate through bits in the register
    def __iter__(self):
        for i in range(self.size):
            yield self._state[self.bit(i)]


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
                        product &= self._state[self.bit(idx)]
                summation ^=  product
            nextState[(self.bit(bitIdx))] = summation
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
                
            
    #return the period of the register:
    def period(self, lim = 2**18):
        first_state = self._state
        self.clock()
        count = 1
        while self._state != first_state:
            self.clock()
            count += 1
            if count > lim:
                return None
        self.reset()
        return count
