from ProductRegisters.Tools.RootCounting.JordanPartition import JP_solve

class JordanSet:
    def __init__(self,roots,m):
        self.roots = roots # Dictionary 
        self.m = m         # Integer

    # Jordan set x Jordan set = Set of JordanSets
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
        new_roots
        
        multiplicities = JP_solve(self.m,other.m,2)
        return set([JordanSet(new_roots, m[0]) for m in multiplicities])

    def isFull(self):
        for b in self.roots:
            if self.roots[b] != b:
                return False
        return True
    
    def __str__(self):
        return  "<" +  ", ".join(f"{k}:{v}" for k,v in self.roots.items()) + f" ({self.m})>"
    
    def __copy__(self):
        return JordanSet({k:v for k,v, in self.roots.items()}, self.m)
