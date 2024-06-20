from itertools import chain, combinations
from math import log

#first 47 mersenne exponents (well beyond what is ever needed) from the OEIS:
mersenne_exponents = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,
          2203,2281,3217,4253,4423,9689,9941,11213,19937,
          21701,23209,44497,86243,110503,132049,216091,
          756839,859433,1257787,1398269,2976221,3021377,
          6972593,13466917,20996011,24036583,25964951,
          30402457,32582657,37156667,42643801,43112609]

#renaming this for smalled reference
MerExp = mersenne_exponents

#small node class for DP method
class Node:
    def __init__(self,val,left,right, count = 0):
        self.val = val
        self.left = left
        self.right = right
        self.count = count

    @classmethod
    def leaf(self):
        return Node(None,None,None)


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

    return left_sol + right_sol

#O(nk) < O(nlogn)
def build_solution_table(lim):
    #select only necessary mersenne exponents:
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

    return table

# O(mlogm + n^2)
def _table_analysis(nums, sort):
    out = []
    num_exponents = 0
    while (num_exponents < len(MerExp) and max(nums) >= MerExp[num_exponents]):
        num_exponents += 1
    
    table = build_solution_table(max(nums))
    for n in nums:
        tree = table[n][num_exponents]
        solutions = unwrap(tree)
        sols_and_epr = [(s, expected_period_ratio(s)) for s in solutions]

        if sort:
            sols_and_epr = sorted(sols_and_epr, key = lambda x: x[1], reverse=True)

        out.append((n,sols_and_epr))
    return out

def list_possible(nums):
    # O(nlogn) call dominates this function
    table = build_solution_table(max(nums)) 

    # Scan through the last column to find nums with any solutions.
    row_length = len(table[0])
    return [n for n in nums if table[n][row_length-1].count > 0]



# single number brute force check:
def _single_brute_force(target):
    #select only necessary mersenne exponents
    i = 0
    while (i<len(MerExp) and target >= MerExp[i]): 
        i += 1
    MerTable = MerExp[:i]

    out = []

    for sset in powerset(MerTable):
        if sum(sset) == target:
            out.append(sset)
    return out

# more performant for small numbers, but scales worse.
# O(n^2) - idk which is better?
def _brute_force_analysis(ns, sort):
    out = []
    for n in ns:
        sols = _single_brute_force(n)
        sols_and_epr = [(s, expected_period_ratio(s)) for s in sols]

        #if sorted by EPR:
        if sort:
             sols_and_epr = sorted(sols_and_epr, key = lambda x: x[1], reverse=True)

        out.append((n,sols_and_epr))
    return out

# choose between methods:
def mersenne_combinations(targets, build_table = True, sort = True):
    if build_table:
        return _table_analysis(targets, sort)
    else: 
        return _brute_force_analysis(targets, sort)


# used to print the output of the find_solutions methods
def pretty_print_constructions(constructions):
    for num in constructions:
        print(f"{num[0]}:\n")
        for sol in num[1]:
            power_10 = log(sol[1],10) + num[0]*log(2,10)
            coef = 10**(power_10 % 1)
            exponent = int(power_10 // 1)
            print(f"\t{sol[0]}:")
            print(f"\t\tApprox. Expected Length: {coef} x 10^{exponent}")
            print(f"\t\tRatio to Full Period: {sol[1]}")








# Helper methods:
def assert_no_repeats(sizes):
    already_seen = set()
    for s in sizes:
        if s in already_seen:
            raise ValueError("Configuration-based properties cannot be determined with repeated sizes")
        already_seen.add(s)

def max_period(sizes):
    period = 1
    already_seen = set()
    for s in sizes:
        if s == 1 or s in already_seen:
            period *= 2
        else:
            already_seen.add(s)
            period *= (2**s-1)
    
    return period


#EPR using direct calculation, as a check:
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def cycle_sizes(sizes):
    assert_no_repeats(sizes)
    primes = [2**a-1 for a in sizes]
    out = []
    for combo in powerset(primes):
        prod = 1
        for c in combo:
            prod *= c
        out.append(prod)
    return out

def epr_brute_force(sizes):
    assert_no_repeats(sizes)
    cycles = cycle_sizes(sizes)
    tot = sum(cycles)
    expected_val = 0
    for c in cycles:
        expected_val += (c/tot)**2
    return expected_val

def expected_period_brute_force(sizes):
    assert_no_repeats(sizes)
    cycles = cycle_sizes(sizes)
    tot = sum(cycles)
    expected_val = 0
    for c in cycles:
        expected_val += (c**2)/tot
    return expected_val



# faster calculation proved in the paper:
def expected_period(sizes):
    assert_no_repeats(sizes)
    numerator = 1
    for s in sizes: numerator *= (((2**s)-1)**2 + 1)
    denominator = 1
    for s in sizes: denominator *= (2**s)
    return numerator/denominator

def expected_period_ratio(sizes):
    assert_no_repeats(sizes)
    numerator = 1
    for s in sizes: numerator *= (((2**s)-1)**2 + 1)
    denominator = 1
    for s in sizes: denominator *= (2**(2*s))
    return numerator/denominator