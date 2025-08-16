from PyPR.BooleanLogic.Gates import NOT
from PyPR.BooleanLogic.BooleanANF import BooleanANF
from functools import cmp_to_key
import time

import os
os.system('')

# term <-> frozenset cheat sheet:
#   - (TERM)    ---   (SET)       ---   (PYTHON)
#   - times     ---   union       ---   |
#   - divides   ---   subset      ---   <=
#   - divide    ---   set minus   ---   -
#   - coprime   ---   empty intersection


# Helper Functions for readability:
def less_than_or_equal(term_1, term_2):
    if len(term_1) < len(term_2):
        return True
    if len(term_1) > len(term_2):
        return False
    return (
        sorted(term_1, reverse=True) <= 
        sorted(term_2, reverse=True)
    )

def strictly_less_than(term_1, term_2):
    return (not term_1==term_2) and less_than_or_equal(term_1,term_2)

def monomial_compare(term_1,term_2):
    """
    Compare two monomials (frozensets):
    returns:
      -   0  if term1 == term2
      -   1  if term1  > term2
      - (-1) if term1  < term2
    """
    # handle equals case first:
    if term_1 == term_2:
        return 0

    # grading based on degree
    if len(term_1) > len(term_2):
        return 1
    if len(term_1) < len(term_2):
        return -1
    
    # differentiate based on sorted:
    if (
        sorted(term_1, reverse=True) >=
        sorted(term_2, reverse=True)
    ):
        return 1
    else:
        return -1
monomial_order = cmp_to_key(monomial_compare)

def lead_term(f) -> frozenset | None:
    """
    Get the lead term of a BooleanANF
    """
    if f.terms:
        return max(f.terms, key=monomial_order)
    else:
        return None

def coprime(term_1, term_2):
    # NOT IN THEIR DEFINITION:
    # if (not term_1) or (not term_2):
    #     raise ValueError
    
    if term_1 & term_2:
        return False
    return True



# Algorithm implementation:
# (intended to match paper pseudocode as much as possible)

def bounded_basis(t,f,sigma,bounded_bases):
    # fail when lead_term(p) = t = None
    if t == None:
        return None

    # otherwise, find a bounded base which can be scaled-up
    # (This is guaranteed to exist, so it shouldn't fail)
    for sigma_prime,c,p in bounded_bases:
        if (
            lead_term(p) <= t and 
            coprime(t - lead_term(p), sigma_prime) and
            strictly_less_than((t - lead_term(p)) | sigma_prime, sigma)
        ):
            return (sigma_prime,c,p)

def gen(sigma, f, bounded_bases):
    # initialize c and p:
    sigma_0, c_0, p_0 = max(
        [(sigma_prime, c_prime, p_prime) for (sigma_prime, c_prime, p_prime) in bounded_bases if sigma_prime <= sigma],
        key = lambda x: monomial_order(x[0])
    )
    c = c_0 * BooleanANF([sigma-sigma_0])
    p = p_0 * BooleanANF([sigma-sigma_0])

    # end early on failure:
    basis = bounded_basis(lead_term(p), f, sigma, bounded_bases)
    if coprime(sigma-sigma_0, lead_term(p_0)):
        if basis == None:
            return None
        
    # reduce leading terms in p, while maintaining lead_term(c) = sigma
    while basis != None:
        sigma_prime, c_prime, p_prime = basis
        update_poly = BooleanANF([lead_term(p)-lead_term(p_prime)]) #type: ignore
        c += c_prime * update_poly
        p += p_prime * update_poly
        basis = bounded_basis(lead_term(p), f, sigma, bounded_bases)
    return sigma, c, p

def potential_critical_characters(basis, bounded_bases, length_limit):
    sigma_1, c_1, p_1 = basis
    p_1_lead_term = lead_term(p_1)
    if p_1_lead_term == None:
        raise ValueError("Zero Polynomial")

    todo = set()
    for sigma_2, c_2, p_2 in bounded_bases:
        # reuse variables to reduce lead_term() calls:
        p_2_lead_term = lead_term(p_2)
        mu_1 = p_2_lead_term - p_1_lead_term # type: ignore
        mu_2 = p_1_lead_term - p_2_lead_term # type: ignore

        if (
            # coprime check:
            not(mu_1 & sigma_1) and
            not(mu_2 & sigma_2) and 
            (mu_1 | sigma_1) != (mu_2 | sigma_2)
        ):
            candidate = max(
                (mu_1 | sigma_1), (mu_2 | sigma_2),
                key = monomial_order
            )

            # filtering here helps keep the queue small
            if len(candidate) <= length_limit:
                todo.add(candidate)
    
    for var in p_1_lead_term:
        if var not in sigma_1:
            todo.add(sigma_1 | frozenset([var]))

    return todo
            
def annihilators(f, annihilator_only=False, verbose = True):
    if verbose:
        print("Starting\n\n\n")
        start_time = time.time()

    f = BooleanANF.from_BooleanFunction(f)
    degrees = (f.degree(),0)
    annihilators = []
    bounded_bases = []

    starting_basis = (frozenset(),BooleanANF([True]),f)
    todo = potential_critical_characters(starting_basis,bounded_bases,max(degrees))
    bounded_bases.append(starting_basis)

    count = 0
    while todo:        
        # pop minimum sigma (maybe PQ here later?)
        sigma_todo = min(todo, key=monomial_order)
        todo.remove(sigma_todo)

        # update printed status:
        if verbose:
            print(
                f"\r\x1B[3A" +
                f"|    Iteration: {count+1}\n" + 
                f"|    Progress: {min(len(sigma_todo),max(degrees))}/{max(degrees)}\n" +
                f"|    Currently Found: (Degrees: {degrees} / Dimension: {len(annihilators)})\n" +
                f"|    Queue/Basis Size: {len(todo)}/{len(bounded_bases)}",
                end=''
            )
            count += 1
        
        if len(sigma_todo) > max(degrees):
            break

        # skip iteration if needed:
        new_basis = gen(sigma_todo,f,bounded_bases)
        if new_basis == None:
            continue

        # update degrees / annihilators:
        sigma,c,p = new_basis

        # Normal degree update
        if not annihilator_only:
            old_degs = sorted(degrees,reverse=True)
            new_degs = sorted((c.degree(),p.degree()), reverse=True)
            if (new_degs < old_degs):
                degrees = (c.degree(),p.degree())
                annihilators = [c]
            elif new_degs == old_degs:
                annihilators.append(c)

        # check if zero character:
        if not(p.terms): 
            # Ann-Only degree update
            if annihilator_only:
                if c.degree() < degrees[0]:
                    degrees = (c.degree(),0)
                    annihilators = [c]
                elif c.degree() == degrees[0]:
                    annihilators.append(c)
        else:
            # not a zero character => critical character:
            todo |= potential_critical_characters((sigma,c,p),bounded_bases,max(degrees))
            bounded_bases.append((sigma,c,p))

    annihilators = [a.to_BooleanFunction() for a in annihilators if a.terms]

    if verbose:
        print(f"\nFinished")
        print(f"Total time -- {time.time()-start_time}")

    return (degrees, annihilators)

def ann_iterator(
    f, annihilator_only=False, verbose = True, yield_rate = 100
):
    if verbose:
        print("Starting\n\n\n")
        start_time = time.time()

    f = BooleanANF.from_BooleanFunction(f)
    degrees = (f.degree(),0)
    annihilators = []
    bounded_bases = []

    starting_basis = (frozenset(),BooleanANF([True]),f)
    todo = potential_critical_characters(starting_basis,bounded_bases,max(degrees))
    bounded_bases.append(starting_basis)

    count = 0
    while todo:        
        # pop minimum sigma (maybe PQ here later?)
        sigma_todo = min(todo, key=monomial_order)
        todo.remove(sigma_todo)

        # update printed status:
        if verbose:
            print(
                f"\r\x1B[3A" +
                f"|    Iteration: {count+1}\n" + 
                f"|    Progress: {min(len(sigma_todo),max(degrees))}/{max(degrees)}\n" +
                f"|    Currently Found: (Degrees: {degrees} / Dimension: {len(annihilators)})\n" +
                f"|    Queue/Basis Size: {len(todo)}/{len(bounded_bases)}",
                end=''
            )
            count += 1

        if count % yield_rate == 0:
            yield (degrees, annihilators, todo, bounded_bases)
        
        if len(sigma_todo) > max(degrees):
            break

        # skip iteration if needed:
        new_basis = gen(sigma_todo,f,bounded_bases)
        if new_basis == None:
            continue

        # update degrees / annihilators:
        sigma,c,p = new_basis

        # Normal degree update
        if not annihilator_only:
            old_degs = sorted(degrees,reverse=True)
            new_degs = sorted((c.degree(),p.degree()), reverse=True)
            if (new_degs < old_degs):
                degrees = (c.degree(),p.degree())
                annihilators = [c]
            elif new_degs == old_degs:
                annihilators.append(c)

        # check if zero character:
        if not(p.terms): 
            # Ann-Only degree update
            if annihilator_only:
                if c.degree() < degrees[0]:
                    degrees = (c.degree(), 0)
                    annihilators = [c]
                elif c.degree() == degrees[0]:
                    annihilators.append(c)
        else:
            # not a zero character => critical character:
            todo |= potential_critical_characters((sigma,c,p),bounded_bases,max(degrees))
            bounded_bases.append((sigma,c,p))

    annihilators = [a.to_BooleanFunction() for a in annihilators if a.terms]

    if verbose:
        print(f"\nFinished")
        print(f"Total time -- {time.time()-start_time}")

    return (degrees, annihilators)

