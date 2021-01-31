from NLMIPR import ProductRegister
from BitVector import BitVector
from itertools import islice

class StandardMIPR():  

    def __init__(self, seed, fn):
        self.pr = ProductRegister(seed,fn)
        self.seed = self.pr.state
        self.size = fn.size

    def pad(self, block):
        block = list(block)
        block += ([0] * (self.size - len(block)))
        return block

    def chunks(self, inpt, n):
        for i in range(0, len(inpt), n):
            yield inpt[i:i + n]

    def multU(self,inpt):
        self.pr.state = inpt
        self.pr.state = next(self.pr)
        out = self.pr.state
        self.pr.state = self.seed
        return out

    def runHash(self,inpt, showInternalState = False):
        blockGen = self.chunks(list(inpt), self.size)

        for _ in range((len(inpt)-1)//self.size + 1):
            block = self.pad(next(blockGen))
            self.pr.state = next(self.pr)
            self.pr.state ^= BitVector(bitlist = block)
   
            if showInternalState: print(self.pr.state)

        outHash = self.pr.state
        
        #reset for next hash
        self.pr.state = self.seed

        return outHash
    
    def createNoise(self, inpt):
        blockGen = self.chunks(inpt, self.size)
        blocks = [self.pad(block) for block in blockGen]

        out = blocks[0]

        for i in range(1,len(blocks)):
            new = BitVector(bitlist = blocks[i]) ^ self.multU(BitVector(bitlist = blocks[i-1]))
            out += list(new)

        out += self.multU(blocks[-1])
        self.pr.state = self.seed
        return out

    def mapU(self, inpt):
        blockGen = self.chunks(inpt, self.size)
        blocks = [self.pad(block) for block in blockGen]
        out = []
        for i in range(len(blocks)):
            out += list(self.multU(blocks[i]))
        return out
