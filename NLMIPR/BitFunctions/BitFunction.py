from BitVector import BitVector
from random import randrange, randint, sample

class BitFunction:
    def __init__(self, fn):
        self.fn = fn
        self.size = len(fn)

    def __getitem__(self, idx): return self.fn[idx]

    def __setitem__(self, idx, val): self.fn[idx] = val

    def __len__(self): return len(self.fn)

    def __repr__(self):
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
    
    def gateSummary(self):
        xors = 0
        ands = 0
        for b in self.fn:
            for term in b:
                xors += 1
                for element in term:
                    ands += 1
                ands -= 1
        print(f"\nMaximum # of XORS: {xors}")
        #print(f"Maximum # of ANDS: {ands}")
        return (xors, ands)