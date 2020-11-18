from BitVector import BitVector
from .BitFunction import BitFunction
from random import randint, sample

class Composition(BitFunction):
    def __init__(self, fnList):
        self.divisions = []
        shift = 0
        currFn = []
        for bitFn in fnList[::-1]:
            shift = len(currFn)
            self.divisions.append(shift)
            shiftedFn = []
            for fn in bitFn.fn:
                newfn = []
                for term in fn:
                    newTerm = []
                    for i in term:
                        #shift every term by i:
                        newTerm.append(i+shift) 
                    newfn.append(newTerm)
                shiftedFn.append(newfn)
            currFn += shiftedFn
        
        self.fn = currFn
        self.size = len(self.fn)

    #same as bitfunction str method, but adds block divisions
    def __str__(self):
        outstr = ""
        for i in range(len(self.fn)-1,-1,-1):
            outstr += str(i) + "="
            for term in self.fn[i]: 
                if len(term) == 1:
                    outstr += str(term[0]) + ","
                else:
                    outstr += "(" + ",".join(str(t) for t in term) + "),"
            outstr = outstr[:-1] + ";\n"
            if i in self.divisions:
                outstr += ("-"*30 + "\n") #adds a line of - at division markers
        return outstr

    def generateNonlinearity(self, maxTerms = 4, maxTermSize = 5):
        #iterate through the blocks:
        for d in range(len(self.divisions)-1):

            #Calculate the range of the upperblocks
            upperBlocks = range(self.divisions[d+1], self.size)

            #iterate through the bits in the current block
            for bit in range(self.divisions[d], self.divisions[d+1]):
                numTerms = randint(1,maxTerms) #randomize number of terms.
                newTermSets = []
                newTerms = []
                while len(newTerms) < numTerms:
                    termSize = randint(1,maxTermSize) #randomize number of bits in each term. 
                    newTerm = sample(upperBlocks,termSize) #select the bits in each NL term.
                    
                    if set(newTerm) not in newTermSets: #check to make sure we do not have duplicates
                        newTermSets.append(set(newTerm))
                        newTerms.append(newTerm)
                self.fn[bit] += newTerms #append new term to bit