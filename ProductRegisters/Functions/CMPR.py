from BitVector import BitVector

from ProductRegisters.Functions import FeedbackFunction
from ProductRegisters.Functions import MPR

from ProductRegisters.ANF import ANF
from ProductRegisters.Tools.RootCounting import RootExpression

from random import randint, sample
import numpy as np

from functools import cached_property

class CMPR(FeedbackFunction):

    def __init__(self, FnList):
        #list of [0..sizes.. size]
        self.num_subregisters = len(FnList)
        self.divisions = [] 

        shift = 0
        currFn = []
        for MPR in FnList[::-1]:
            
            self.divisions.append(shift)
            shiftedFn = []

            for fn in MPR.anf:
                newfn = []
                for term in fn:
                    newTerm = []
                    if type(term) == bool:
                        newTerm = term
                    else:
                        for i in term:
                            #shift every term:
                            newTerm.append(i+shift) 
                    newfn.append(newTerm)
                shiftedFn.append(newfn)
            currFn += shiftedFn

            #update the amount to shift
            shift += MPR.size
       
        self.anf = [ANF(fn) for fn in currFn]
        self.size = len(self.anf)
        self.divisions.append(self.size)

    @cached_property 
    def blocks(self):
        block_list = []
        for d in range(len(self.divisions)-1):
            bits = list(range(self.divisions[d], self.divisions[d+1]))
            block_list.append(bits)
        return block_list[::-1]

    #overwrite bitfn str method, adds block divisions
    def __str__(self):
        outstr = ""
        for i in range(len(self.anf)-1,-1,-1):
            outstr += str(i) + "="
            outstr += str(self.anf[i]) + ";\n"

            #adds a line of - at division markers
            if i in self.divisions[1:-1]:
                outstr += ("-"*30 + "\n") 
        return outstr[:-1]


    # TODO: This samples poorly for linear complexity:
    #  - we whould rewrite with better heuristics
    def generateChaining(self, maxTerms = 4, maxTermSize = 5):
        #iterate through the blocks:
        for d in range(len(self.divisions)-2):

            #Calculate the range of the upperblocks
            upperBlocks = range(self.divisions[d+1], self.size)
            blockSize = len(upperBlocks)
            #iterate through the bits in the current block
            for bit in range(self.divisions[d], self.divisions[d+1]):
                
                #an ANF of the new terms being added
                newTerms = ANF()
                
                #randomly select whether or not there is a bool:
                if randint(0,1):
                    newTerms.add(frozenset())

                # NO CLUE why I had this but I did so im not touching rn.
                # validate numTerms (actual number of terms)
                # numTerms = None
                # while numTerms is None or numTerms > blockSize:
                    
                numTerms = randint(1,maxTerms) #randomize number of terms.
                
                while len(newTerms) < numTerms:
                    #randomize number of bits in each term:
                    termSize = randint(1, min(maxTermSize, blockSize))  
                    #select the bits in each NL term:
                    newTerm = sample(upperBlocks,termSize) 
                    #union the new term into the ANF:
                    newTerms.add(newTerm)

                self.anf[bit] += newTerms #append new term to bit


    def fixpoint(self):

        #define inner GF(2) gaussian elimination subroutine:
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

        #create a key of known values (bit index: value):
        key = {}
        
        #for each subregister (descending order):
        for block in self.blocks:

            #offset mapping for each variable:
            offset = {bit:(bit-block[0]) for bit in block}
            
            #number of variables in this xor system
            b_size = len(block)
            xors = np.zeros([b_size, b_size+1],dtype = np.uint8)
            for b in block:
                xors[offset[b],offset[b]] = 1

            #for each bit, calculate chaining constant with key:
            for bit in block:
                
                cur_val = 0
                anf = self.anf[bit]

                for term in anf:
                    if term == frozenset():
                        cur_val ^= 1
                    elif len(term)==1 and tuple(term)[0] in block:
                        xors[offset[bit],offset[tuple(term)[0]]] ^= 1
                    else:
                        prod = 1
                        for var in term:
                            prod &= key[var]
                        cur_val ^= prod

                xors[offset[bit],b_size] = cur_val
            XorSolve(xors)

            #backsubstitution and filling out key:
            for bit in block[::-1]:
                eq = xors[offset[bit],offset[bit]+1:-1]
                res = xors[offset[bit],b_size]
                #print(eq, res)
                for i in range(len(eq)):
                    if eq[i] == 1:
                        res ^= key[i+bit+1]
                key[bit]=res
           
        return BitVector(bitlist = [key[i] for i in range(self.size)][::-1])
    
    @property
    def RootExpressions(self):
        block_map = {}   # map: bit -> block
        expr_table = {}  # map: block -> expression

        blocks = self.blocks
        for b in range(len(blocks)):
            
            #create inital RE for each block
            expr_table[b] = RootExpression([{len(blocks[b]) : 1}]) 

            for bit in blocks[b]:
                block_map[bit] = b # map each bit to its corresponding block            

            unified = RootExpression() #the unified expression for the entire block
            for bit in blocks[b]:
                for term in self.anf[bit]:
                    #create RE for the term:
                    if term == frozenset():
                        pass
                    else:
                        new_term = None
                        for val in term:
                            if new_term == None: #use the first value as the base for the term
                                new_term = expr_table[block_map[val]]
                            else:
                                new_term *= expr_table[block_map[val]]

                    #insert unify the RE for the new term with the rest of the block
                    unified += new_term

            expr_table[b] = unified
        out_list= [expr_table[block_map[bit]] for bit in range(self.size)]
        return out_list


    def estimate_LC(self, output_bit, locked_list = None, benchmark = False):
        # locked-list is used to cancel effects of the locked registers.
        # the locked list contains the sizes of the locked MPRs

        from time import time_ns
        t1 = time_ns()
        bitRE = self.RootExpressions[output_bit]
        t2 = time_ns()

        if benchmark:
            print(f"Root Expression Generation: {len(bitRE.anf)} terms generated in {t2-t1} ns")

        t1 = time_ns()
        #get the length of the block bit is in.
        blocks = self.blocks
        for block in blocks:
            if output_bit in block:
                blockLen = len(block)
        
        #lower the minimum if this block is locked:
        if locked_list and blockLen in locked_list:
            blockLen = 1

        #add 1 to include the "True" / 1 in GF(2) value we ignore:
        upper = bitRE.upper(locked_list) + 1

        #subtract safety_factor * size to account for some natural degeneracies:
        # TODO: need a better incorporation of the safety factor.

        lower = max(blockLen,bitRE.lower(locked_list))
        t2 = time_ns()
        if benchmark:
            print(f"Terms evaluated in {t2-t1} ns")

        return (lower,upper)