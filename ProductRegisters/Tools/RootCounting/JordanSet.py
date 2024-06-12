from ProductRegisters.Tools.RootCounting.JordanPartition import JP_solve

class JordanSet:
    def __init__(self,roots, mults):
        self.roots = roots # Dictionary 
        self.mults = mults # set

    # Jordan set x Jordan set = Jordan set
    def __mul__(self,other):
        new_roots = {}
        for basis,count in self.roots.items():
            if basis in new_roots:
                new_roots[basis] += count
            else:
                new_roots[basis] = count
        
        for basis,count in other.roots.items():
            if basis in new_roots:
                new_roots[basis] += count
            else:
                new_roots[basis] = count
        
        # reduce the counts in each multiplication, so that they dont get big.
        for basis in new_roots:
            new_roots[basis] = min(basis,new_roots[basis])
        
        multiplicities = set()
        for m1 in self.mults:
            for m2 in other.mults:
                multiplicities |= set([mult for mult,count in JP_solve(m1,m2,2)])

        return JordanSet(new_roots, multiplicities)

    def isFull(self):
        for b in self.roots:
            if self.roots[b] != b:
                return False
        return True
    
    def __str__(self):
        return  (
            "<" + ", ".join(f"{k}:{v}" for k,v in self.roots.items()) + 
            " (" + ",".join(str(m) for m in self.mults) +
            ")>"
        )
        
    def __copy__(self):
        return JordanSet({k:v for k,v, in self.roots.items()}, set(m for m in self.mults))
