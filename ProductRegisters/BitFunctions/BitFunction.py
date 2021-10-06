from BitVector import BitVector
from random import randrange, randint, sample
from copy import deepcopy
import json

class BitFunction:
    def __init__(self, fn):
        self.fn = fn
        self.size = len(fn)

    def __getitem__(self, idx): return self.fn[idx]

    def __setitem__(self, idx, val): self.fn[idx] = val

    def __len__(self): return len(self.fn)

    def __eq__(self, other): return self.fn == other.fn

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
        return outstr[:-1]

    

    #new function that generates the same states in reverse bit order
    def flip(self):
        newFn = []
        for bitFn in self.fn[::-1]:
            newBitFn = []
            for term in bitFn:
                newTerm = []
                for var in term:
                    newTerm.append(self.size-var-1)
                newBitFn.append(newTerm)
            newFn.append(newBitFn)
        self.fn = newFn



    #return the number of gates before any optimization
    def gateSummary(self):
        xors = 0
        ands = 0
        for b in self.fn:
            for term in b:
                xors += 1
                for element in term:
                    ands += 1
                ands -= 1
        #print(f"\nMaximum # of XORS: {xors}")
        #print(f"Maximum # of ANDS: {ands}")
        return (xors, ands)

    #return true if the function is linear
    def isLinear(bitFn, allowAfine = False):
        valid = True
        for ls in bitFn.fn:
            for l in ls:
                if type(l) == bool:
                    valid &= allowAfine
                elif len(l) != 1:
                    valid = False
        return valid


    #to/from JSON methods (inherited by all bitfunction types)
    def toJSON(self, filename):
        dct = {}
        dct["Type"] = type(self).__name__
        dct["Data"] = self.__dict__
        with open(filename, "w") as f:
            f.write(json.dumps(dct))

    @classmethod
    def fromJSON(self, filename):
        with open(filename, "r") as f:
            dct = json.load(f)

            #if wrong fn type
            if dct["Type"] != self.__name__:
                expected = self.__name__
                actual = dct["Type"]
                raise ValueError(
                    f"Expected function type {expected}, but got {actual}"
                    )
            
            #otherwise remake the FN
            fn = object.__new__(self)
            fn.__dict__ = dct["Data"]
        return fn
    
