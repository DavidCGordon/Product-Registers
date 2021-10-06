from BitVector import BitVector
from .BitFunction import BitFunction
from random import randint, sample
import numpy as np

class Composition(BitFunction):
    def __init__(self, fnList):
        self.divisions = []
        self.subregister_types = []

        shift = 0
        currFn = []
        for bitFn in fnList[::-1]:
            shift = len(currFn)
            self.divisions.append(shift)
            self.subregister_types.append(type(bitFn).__name__)
            shiftedFn = []
            for fn in bitFn.fn:
                newfn = []
                for term in fn:
                    newTerm = []
                    if type(term) == bool:
                        newTerm = term
                    else:
                        for i in term:
                            #shift every term by i:
                            newTerm.append(i+shift) 
                    newfn.append(newTerm)
                shiftedFn.append(newfn)
            currFn += shiftedFn
        
        self.fn = currFn
        self.size = len(self.fn)

    #overwrite bitfn str method, adds block divisions
    def __str__(self):
        outstr = ""
        for i in range(len(self.fn)-1,-1,-1):
            outstr += str(i) + "="
            for term in self.fn[i]:
                if type(term) == bool:
                    outstr += str(term) + ","
                elif len(term) == 1:
                    outstr += str(term[0]) + ","
                else:
                    outstr += "(" + ",".join(str(t) for t in term) + "),"
            outstr = outstr[:-1] + ";\n"
            if i in self.divisions[1:]:
                outstr += ("-"*30 + "\n") #adds a line of - at division markers
        return outstr

    def generateChaining(self, maxTerms = 4, maxTermSize = 5):
        #iterate through the blocks:
        for d in range(len(self.divisions)-1):

            #Calculate the range of the upperblocks
            upperBlocks = range(self.divisions[d+1], self.size)
            blockSize = len(upperBlocks)

            #iterate through the bits in the current block
            for bit in range(self.divisions[d], self.divisions[d+1]):
                
                newTermSets = []
                newTerms = []
                
                #randomly select whether or not there is a bool:
                if randint(0,1):
                    newTerms.append(True)

                #validate numTerms (actual number of terms)
                numTerms = None
                while numTerms is None or numTerms > blockSize:
                    numTerms = randint(1,maxTerms) #randomize number of terms.
                
                while len(newTerms) < numTerms:
                    termSize = randint(1, min(maxTermSize, blockSize)) #randomize number of bits in each term. 
                    newTerm = sample(upperBlocks,termSize) #select the bits in each NL term.
                    
                    if set(newTerm) not in newTermSets: #check to make sure we do not have duplicates
                        newTermSets.append(set(newTerm))
                        newTerms.append(newTerm)
                self.fn[bit] += newTerms #append new term to bit
    
    #assumes all subregisters are mersenne
    def fixpoint(self):
        for sub_type in self.subregister_types:
            if sub_type not in ["Mersenne","Fibonacci"]:
                raise ValueError("Feedback function must consist of Mersenne FPRs only.")
            
        #define GF(2) gaussian elimination subroutine:
        def XorSolve(eqs):
            num_row, num_col = eqs.shape
            p_row = 0
            p_col = 0

            while p_row < num_row and p_col < num_col:
                offset = np.argmax(eqs[p_row:, p_col])
                if eqs[p_row+offset, p_col] == 0:
                    p_col += 1
                else:
                    #swap rows:
                    eqs[[p_row,p_row+offset]] = eqs[[p_row+offset,p_row]]

                    need_changing = eqs[p_row+1:, p_col][:,np.newaxis]
                    change = eqs[p_row] * need_changing

                    eqs[p_row+1:] ^= change
                    
                    p_row += 1
                    p_col += 1
            return eqs
    
        #create a key of known values (bit-#: value):
        key = {}
        
        #calculate bit ranges for each block and flip
        divs = self.divisions + [self.size]
        blocks = []
        for d in range(len(self.divisions)):
            blocks.append(list(range(divs[d],divs[d+1])))
        blocks = blocks[::-1]
        
        #for each subregister (descending order):
        for block in blocks:

            #offset mapping for each variable:
            offset = {bit:(bit-block[0]) for bit in block}
            
            #number of variables in this xor system
            b_size = len(block)
            xors = np.zeros([b_size, b_size+1],dtype = np.uint8)
            for b in block:
                xors[offset[b],offset[b]] = 1


            #for each bit:
            for bit in block:
                
                cur_val = 0
                anf = self.fn[bit]

                for var in anf:
                    if type(var) == bool:
                        cur_val ^= int(var)
                    elif len(var)==1 and var[0] in block:
                        xors[offset[bit],offset[var[0]]] ^= 1
                    else:
                        prod = 1
                        for subvar in var:
                            prod &= key[subvar]
                        cur_val ^= prod
                xors[offset[bit],b_size] = cur_val
            XorSolve(xors)

            #backsubstitution:
            for bit in block[::-1]:
                eq = xors[offset[bit],offset[bit]+1:-1]
                res = xors[offset[bit],b_size]
                #print(eq, res)
                for i in range(len(eq)):
                    if eq[i] == 1:
                        res ^= key[i+bit+1]
                key[bit]=res
           
        return BitVector(bitlist = [key[i] for i in range(self.size)][::-1])
