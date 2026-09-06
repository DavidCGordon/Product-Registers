import numpy as np
import numba as nb
import galois as gl
from galois import egcd
import sympy

import math

from PyPR.FeedbackFunctions import MPR
from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

from itertools import product

def modular_compose(remainders,moduli):
    if (len(remainders) != len(moduli)):
        raise ValueError()
    if len(remainders) == 0:
        raise ValueError()
 
    curr_remainder = remainders[0]
    curr_modulus = moduli[0]

    #Merge the solution with remaining equations
    for i in range(1,len(remainders)):
        merge_rem = remainders[i]
        merge_mod = moduli[i]

        common_factor = math.gcd(curr_modulus, merge_mod)
        if (curr_remainder % common_factor != merge_rem % common_factor):
            raise ValueError() # No solution exists

        # Merge the two equations
        _, bezout_current, bezout_merge = egcd(
            curr_modulus // common_factor, 
            merge_mod //  common_factor
        )

        # Potential overflow here:
        new_mod = curr_modulus // common_factor * merge_mod; # LCM of m1 and m2
        curr_remainder = (
            merge_rem * (bezout_current) * (curr_modulus//common_factor) + 
            curr_remainder * (bezout_merge) * (merge_mod//common_factor)
        ) % new_mod
        curr_modulus = new_mod

    return curr_remainder, curr_modulus

def bezout(numbers):
    """ Given Numbers, generate (g, coefficients), such that g is the gcd of numbers
    and coefficents is a bezout-like vector (coefficients dot numbers = g)
    """
    if not numbers:
        raise ValueError()
    if len(numbers) == 1:
        return (numbers[0], [1])

    # Start with the first two numbers
    a = numbers[0]
    b = numbers[1]
    g, x, y = egcd(a, b)
    coefficients = [x, y]

    # Iteratively incorporate the rest of the numbers
    for i in range(2, len(numbers)):
        c = numbers[i]
        # Calculate g' = gcd(g, c) and coefficients s, t such that g*s + c*t = g'
        g_prime, s, t = egcd(g, c)
        
        # Update coefficients for all previous numbers:
        for j in range(len(coefficients)):
            coefficients[j] *= s
        coefficients.append(t)
        g = g_prime
        
    return (g, coefficients)

def order(num, mod):
    if num == 1:
        return 1
    
    power = 1
    value = num
    while value != 1:
        value = (value * num) % mod
        power += 1
    return power

def get_field(denominator):
    return order(2,denominator)

def euler_totient(x):
    if x == 1:
        return 1
    
    output = 1
    for f,m in sympy.factorint(x).items():
        output *= (f**(m-1))*(f-1)
    return output

def carmichael_lambda(x):
    if x == 1:
        return 1
    
    facts = []
    for f,m in sympy.factorint(x).items():
        if f**m in [2,4]:
            facts.append(f**m // 2)
        elif f == 2:
            facts.append((f**(m-2)))
        else:
            facts.append((f**(m-1))*(f-1))
    
    if len(facts) == 1:
        return facts[0]
    
    g = math.gcd(*facts)
    output = 1
    for f in facts:
        output *= f

    return output // g

def num_primitive_polynomials(degree):
    return euler_totient(2**degree-1) // degree

def prime_modulus_primitivity_test(candidate,prime):
    test_numbers = [prime // f for f in sympy.factorint(prime-1).keys()]
    for n in test_numbers:
        if pow(candidate,n,prime) == 1:
            return False
    return True

def primitive_elements(prime,power,lim=1000):
    total_mod = prime**power

    base_set = set()
    for x in range(1,prime):
        if prime_modulus_primitivity_test(x,prime):
            base_set.add(x)

            if len(base_set) == lim:
                break
            
    current_set = set(base_set)
    for degree in range(1,power):
        if len(current_set) >= lim:
            break

        add = [n*(prime**degree) for n in range(1,prime)]
        curr_mod = prime**(degree+1)
        target_order = carmichael_lambda(curr_mod)

        current_set |= set([
            (x+y)%total_mod for x,y in product(current_set,add) 
            if order((x+y)%total_mod,curr_mod) == target_order
        ])
       
    return list(current_set)[:lim]

# NEW STUFF (MOVE LATER)
def simple_precompute(n):
    output = []
    for factor,mult in sympy.factorint(n).items():
        output += [factor]*mult
    return sorted(output)  

def simple_precompute2(n, field, multiples=100):
    n %= (2**field-1)
    conjugates = [(n*(2**i))%(2**field-1) for i in range(field)]

    #for p in [3,5,7,11,13,17,19,23,29][:nprimes]:
    targets = []
    for i in range(multiples):
        targets += [c + i*(2**field-1) for c in conjugates]

    targets = sorted(targets)
    best = simple_precompute(targets[0])
    best_score = sum(best)*len(best)
    for new_target in targets[1:]:
        if 2*(len(bin(new_target))-2)**2 > best_score:
            continue
        best = min(best, simple_precompute(new_target), key= lambda x: sum(x)*len(x))
    return best

def fast_decimate(starter_poly,decimation_factor,path=None):
    size = len(starter_poly)-1
    if path == None:
        path = simple_precompute2(decimation_factor,size)
    
    poly = starter_poly[:]
    F = FeedbackRegister(1,MPR(size,"1"))

    for d in path:
        M = MPR(size,poly)
        M.compile()

        F.fn = M
        F.reset()

        seq = np.array([state[0] for state in F.run(2*size*(d+10))], dtype='uint8')
        decimated_seq = seq[::d]
        lc, poly = berlekamp_massey(decimated_seq)
        poly = poly[::-1] # reverse poly to get primal

    return np.array(poly,dtype='uint8')

def generate_primitive_polynomials(base):
    size = len(base)-1

    # precompute scalars:
    components = list(sympy.factorint(2**size-1).items())
    component_moduli = [f**m for f,m in components]
    # order_2 = {(f,m): order(2,f**m) for (f,m) in components}
    # decycle = max(components, key = lambda x: (order_2[x],x))
    
    exps = []
    paths = []
    for i,(f,m) in enumerate(components):
        # if (f,m) == decycle:
        #     continue
        
        vec = [1]*(len(components))

        best_exp = None
        best_path = None
        choices = primitive_elements(f,m,lim=100)
        for choice in choices:
            vec[i] = choice
            combined = modular_compose(vec,component_moduli)[0]
            path = simple_precompute2(combined,size,multiples=100)
           
            if (best_exp == None) or sum(path) < sum(best_path): #type:ignore
                best_exp = combined
                best_path = path
        
        exps.append(best_exp)
        paths.append(best_path)
    
    # Use precomputed exponents/paths to cycle:
    seen = set()
    maxs = [f**m-2 for f,m in components if (f,m)]
    curr = [0]*len(maxs)
    idx = len(curr)-1
    poly = base
    while idx >= 0:
        if curr[idx] == maxs[idx]:
            curr[idx] = 0
            idx -= 1
            continue

        if tuple(poly) not in seen:
            seen.add(tuple(poly))
            yield poly

        curr[idx] += 1
        poly = fast_decimate(poly,exps[idx],path=paths[idx])
        idx = len(curr)-1

        print(maxs, curr, poly, tuple(poly) in seen, gl.Poly(poly[::-1]).is_primitive()) 

    # yield final polynomial
    if tuple(poly) not in seen:
        seen.add(tuple(poly))
        yield poly
 
print('start')
base = [0]*128
for i in [0,13,45,54,127]:
    base[i] = 1

# base = [0]*64
# for i in [0,13,45,54,127]:
#     base[i] = 1

#base = [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1]

print(euler_totient(2**16-1))
x = [p for p in generate_primitive_polynomials(base)]
print(len(x),num_primitive_polynomials(len(base)-1))
