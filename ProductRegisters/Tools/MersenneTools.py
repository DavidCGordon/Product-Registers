from itertools import chain, combinations
from math import log

class Node:
    def __init__(self,val,left,right, count = 0):
        self.val = val
        self.left = left
        self.right = right
        self.count = count
    
    #def __str__(self): return str(self.val)
    #def __repr__(self): return str(self.val)

    @classmethod
    def leaf(self):
        return Node(None,None,None)


#first 47 mersenne exponents
MerExp = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,
          2203,2281,3217,4253,4423,9689,9941,11213,19937,
          21701,23209,44497,86243,110503,132049,216091,
          756839,859433,1257787,1398269,2976221,3021377,
          6972593,13466917,20996011,24036583,25964951,
          30402457,32582657,37156667,42643801,43112609]

#unwrap and helper fn get all solutions, given the node
def unwrap(node):
    return sorted(unwrap_helper(node), key = lambda x: x[0])

def unwrap_helper(node):
    #define empty lists to initialize vars:
    if node is None: return []
    left_sol = []
    right_sol = []

    #left = solutions not using this value:
    if node.left:
        left_sol = unwrap_helper(node.left)
    #right = solutions including this value:
    if node.right:
        if node.right.val == 0:
            #if this value is the first one create a list for it:
            right_sol = [[node.val]]
        else:
            #append this value to all solutions not using this value:
            right_sol = [x+[node.val] for x in unwrap_helper(node.right)]
    return left_sol+right_sol

def analyzeAll(lim):
    #select only necessary mersenne exponents
    i = 0
    while (i<len(MerExp) and lim >= MerExp[i]): i += 1
    MerTable = MerExp[:i]

    #creates a k+1 by n+1 table, with space for padding:
    table = [[0 for _ in range(len(MerTable) + 1)] for _ in range(lim + 1)]
    for i in range(len(table)): table[i][0] = Node(1,None,None)
    for i in range(len(table[0])): table[0][i] = Node(0,None,None,1)

    #fills in the trees with the proper relationships
    for k in range(1,len(MerTable) + 1):
        exp = MerTable[k - 1]
        for n in range(1,lim + 1):
            table[n][k] = Node(exp,None,None)

            if table[n][k-1].count:
                table[n][k].left = table[n][k-1]
                table[n][k].count += table[n][k-1].count
            if n >= exp and table[n-exp][k-1].count: 
                table[n][k].right = table[n-exp][k-1]
                table[n][k].count += table[n-exp][k-1].count

    return table #O(nk) < O(nlogn) to find table

def analyze(ns):
    out = []
    lim = max(ns)
    k_lim = 0
    while (k_lim<len(MerExp) and lim >= MerExp[k_lim]):
        k_lim += 1
    
    t = analyzeAll(lim)
    for n in ns:
        tree = t[n][k_lim]
        sols = unwrap(tree)
        sols_exps = [(s, exp(s)) for s in sols]
        out.append((n,sols_exps))
    return out

def listPossible(*args):
    if len(args) == 1:
        start = 0
        stop = args[0]
        step = 1
    elif len(args) == 2:
        start = args[0]
        stop = args[1]
        step = 1
    elif len(args) == 3:
        start = args[0]
        stop = args[1]
        step = args[2]
    else: 
        raise ValueError

    table = analyzeAll(stop) #O(nlogn) call dominates this function
    row_length = len(table[0]) 
    rows_in_range = table[start:stop+1:step]
    values = [row[row_length-1].count for row in rows_in_range]
    return [(i+start,v) for (i,v) in enumerate(values) if v][1:]


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def cyc_sizes(sizes):
	cycs = [2**a-1 for a in sizes]
	out = []
	for combo in powerset(cycs):
		prod = 1
		for c in combo:
			prod *= c
		out.append(prod)
	return out
    
def exp(sizes):
	cycs = cyc_sizes(sizes)
	tot = sum(cycs)
	exp = 0
	for c in cycs:
		exp += (c/tot)**2
	return exp

def pprint(ns):
    output = analyze(ns)
    for num in output:
        print(f"{num[0]}:\n")
        for sol in num[1]:
            power_10 = log(sol[1],10) + num[0]*log(2,10)
            coef = 10**(power_10 % 1)
            expon = int(power_10 // 1)
            print(f"\t{sol[0]}:")
            print(f"\t\tApprox. Expected Length: {coef} x 10^{expon}")
            print(f"\t\tRatio to Full Period: {sol[1]}")



#single number check:
def single(target):
    #select only necessary mersenne exponents
    i = 0
    while (i<len(MerExp) and target >= MerExp[i]): i += 1
    MerTable = MerExp[:i]

    out = []

    for sset in powerset(MerTable):
        if sum(sset) == target:
            out.append(sset)
    return out

#more performant for small numbers, but scales worse.
def analyze_(ns, sort = False):
    out = []
    for n in ns:
        sols = single(n)
        sols_and_epr = [(s, exp(s)) for s in sols]

        #if sorted by EPR:
        if sort:
             sols_and_epr = sorted(sols_and_epr, key = lambda x: x[1])

        out.append((n,sols_and_epr))
    return out
