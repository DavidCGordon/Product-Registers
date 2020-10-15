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
                if len(term) == 1:
                    outstr += str(term[0]) + ","
                else:
                    outstr += "(" + ",".join(str(t) for t in term) + "),"
            outstr = outstr[:-1] + ";\n"
        return outstr
