from BitVector import BitVector
from .BitFunction import BitFunction

class Fibonacci(BitFunction):
    def __init__(self, size, topFunc):
        if type(topFunc) == str:
            taps = list(BitVector(intVal = int(topFunc, 16), size = size))
            topFunc = [[(size-idx)%size] for (idx, t) in enumerate(taps) if t == 1][1:]

        self.fn = [[[(i+1)%size]] for i in range(size-1)] + [([[0]]+topFunc)]
        self.size = size
    
    def __str__(self):
        outstr = ""
        for i in range(len(self.fn)-1,-1,-1):
            outstr += ("Bit " + str(i)
                + ": " + str(self.fn[i][0][0]) + " + ("
                + (" + ".join([str(term) for term in self.fn[i]][1:])) + ")\n")
        return outstr
