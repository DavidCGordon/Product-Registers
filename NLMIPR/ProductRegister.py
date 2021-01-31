from BitVector import BitVector

class ProductRegister:
    def __init__(self, seed, fn):
        self.fn = fn
        self.size = len(fn)
        self.state = BitVector(intVal = seed, size = self.size)

    def __iter__(self): return self
    
    def __str__(self): return self.state.__str__()

    def int_val(self): return int(self.state)

    def reverse(self, bitIdx): return (self.size-1-bitIdx % self.size)
    
    def __next__(self):
        nextState = BitVector(intVal = 0, size = self.size)
        
        for bitIdx in range(self.size):
            
            #make the first value of the function True in order to flip the values of the function

            summation = 0
            
            for term in self.fn[bitIdx]:

                if type(term) == bool:
                    product = int(term)
                else:
                    product = 1
                    for idx in term:
                        product &= self.state[self.reverse(idx)]

                summation ^=  product

            nextState[(self.reverse(bitIdx))] = summation
        
        self.state = nextState
        return self.state