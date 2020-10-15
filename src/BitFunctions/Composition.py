from BitVector import BitVector
from .BitFunction import BitFunction

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
                        newTerm.append(i+shift)
                    newfn.append(newTerm)
                shiftedFn.append(newfn)
            currFn += shiftedFn
        
        self.fn = currFn
        self.size = len(self.fn)

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
                outstr += ("-"*30 + "\n")
        return outstr