from PyPR.BooleanLogic.BooleanANF import BooleanANF
from PyPR.BooleanLogic.FunctionInputs import CONST, VAR
from PyPR.BooleanLogic.Gates import XOR

from functools import cmp_to_key

import heapq
import bisect

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

class pq_node:
    def __init__(self,lead_monomial,poly):
        self.poly = poly
        self.lead_term = lead_monomial
    
    # since we dont care in the equality case, we can use leq
    # but heapq checks for __lt__, so thats what we call it
    def __lt__(self,other):
        if len(self.lead_term) < len(other.lead_term):
            return True
        if len(self.lead_term) > len(other.lead_term):
            return False
        return (
            sorted(self.lead_term, reverse=True) <= 
            sorted(other.lead_term, reverse=True)
        )

def lead_term(f) -> frozenset[int] | None:
    """
    Get the lead term of a BooleanANF
    """
    if f.terms:
        return max(f.terms, key=monomial_order)
    else:
        return None



# class GroebnerEquation():
#     def __init__(self, terms = []):
#         self.terms = terms

#     def __mul__(self, scale):
#         termset = set()
#         for term in self.terms:
#             new_term = term | scale
#             if new_term in termset:
#                 termset.remove(new_term)
#             else:
#                 termset.add(new_term)
#         return GroebnerEquation(sorted(termset), key=monomial_order)
        
#     def __add__(self, other):
#         output = []
#         merged = heapq.merge(self.terms, other.terms, key=monomial_order)
#         for i in range(len(merged)):
#             term = merged[i]

#             j = 0
#             while (i+j < len(merged)) and merged[j] == term:
#                 j += 1

#             if (j+1) % 2:
#                 output.append(term)
#         return GroebnerEquation(output)


#     def __lt__(self,other):
#         if len(self.lead_term) < len(other.lead_term):
#             return True
#         if len(self.lead_term) > len(other.lead_term):
#             return False
#         return (
#             sorted(self.lead_term, reverse=True) <= 
#             sorted(other.lead_term, reverse=True)
#         )
    



















class GroebnerEqStore:
    lead_terms: list[frozenset[int]]

    def __init__(self):
        self.num_vars = 0
        self.num_eqs = 0

        self.seen = set()
        self.queue = []
        self.equations = []
        self.lead_terms = []

        self.unknown_vars = set()
        self.solved_vars = {}

    def reduce(self, poly):
        curr_lead = lead_term(poly)

        while curr_lead:
            selected = next(
                (i for i in range(self.num_eqs) if self.lead_terms[i] <= curr_lead), 
                None # default to value to return for empty polynomials
            )

            if selected == None:
                return poly
            
            poly += (
                BooleanANF([curr_lead - self.lead_terms[selected]]) *
                self.equations[selected]
            )
            
            curr_lead = lead_term(poly)
        return poly

    def syzygy(self,i,j):
        return (
            self.equations[i] * BooleanANF([self.lead_terms[j]-self.lead_terms[i]]) +
            self.equations[j] * BooleanANF([self.lead_terms[i]-self.lead_terms[j]]) 
        )
    
    def enqueue_equation(self, equation, extra_const = 0, identifier=None, translate_ANF = True):
        equation = XOR(equation,CONST(extra_const)).compose({
            var: CONST(val) for var,val in self.solved_vars.items()
        })

        self.unknown_vars |= set(equation.idxs_used())
        equation = BooleanANF.from_BooleanFunction(equation)
        
        if equation.terms and equation not in self.seen:
            heapq.heappush(self.queue, pq_node(lead_term(equation), equation))
            self.seen.add(equation)

    def consume_queue(self, num=None, verbose = False):
        # if no num provided, consume the whole queue
        # (-1) will decrement but never reach 0
        if num == None:
            num = -1

        # main loop
        while self.queue and num:
            candidate = heapq.heappop(self.queue).poly
            num -= 1

            # reduce candidate:
            reduced = self.reduce(candidate)
            reduced_lead = lead_term(reduced)

            # don't process 0 equations
            if reduced_lead == None:
                continue

            # return false for contradictions:
            if reduced.terms == frozenset([frozenset()]):
                return False

            # add in new eq:
            insert_idx = bisect.bisect(self.lead_terms,reduced_lead)
            self.lead_terms.insert(insert_idx, reduced_lead)
            self.equations.insert(insert_idx, reduced)
            self.num_eqs += 1

            # self.lead_terms.append(reduced_lead)
            # self.equations.append(reduced)
            # self.num_eqs += 1

            # unit-propagate solved variables:
            if len(reduced_lead) <= 1:
                if verbose: print("unit propagation")
                finished = False
                while not finished:
                    new_vars = []
                    for v in self.unknown_vars:
                        reduced = self.reduce(BooleanANF([[v]]))

                        if reduced == BooleanANF([True]):
                            self.solved_vars[v] = 1
                        elif reduced == BooleanANF([]):
                            self.solved_vars[v] = 0

                    # update unknown vars:
                    for v in new_vars:
                        if verbose: print('now known: ', v, '=', self.solved_vars[v])
                        self.unknown_vars.remove(v)
                    
                    if new_vars:
                        self._simplify()
                    else:
                        finished = True

            # compute/add syzygies:
            count = 0
            for eq_idx in range(self.num_eqs):
                if eq_idx == insert_idx:
                    continue

                # coprime lead terms won't lead to good sysygies
                if not (self.lead_terms[insert_idx] & self.lead_terms[eq_idx]):
                    continue

                s_poly = self.syzygy(insert_idx, eq_idx)
                if s_poly.terms and s_poly not in self.seen:
                    heapq.heappush(self.queue, pq_node(lead_term(s_poly), s_poly))
                    self.seen.add(s_poly)
                    count += 1

            if verbose: print(f"added {count} to queue. Queue length: {len(self.queue)}")
        return True
    
    # TODO: Cleanup and justify????
    def _simplify(self):
        seen_set = set()
        new_eqs = []

        zero_set = frozenset([v for v,val in self.solved_vars.items() if val == 0])
        ones_set = frozenset([v for v,val in self.solved_vars.items() if val == 1])

        # REDUCE EQUATIONS:
        for equation in self.equations:
            reduced_terms = [
                term-ones_set for term in equation
                if not (term & zero_set)
            ]

            filtered_terms = set()
            for term in reduced_terms:
                if term not in filtered_terms:
                    filtered_terms.add(term)
                else:
                    filtered_terms.remove(term)
            
            new_equation = BooleanANF(filtered_terms,fast_init=True)
            if new_equation.terms and new_equation not in seen_set:
                seen_set.add(new_equation)
                new_eqs.append(new_equation)
        
        # May need to add in a way to filter old eqs!!
        # ADD IN LINEAR EQs:
        prefix = [BooleanANF([[v],val]) for v,val in self.solved_vars.items()]
        self.equations = prefix + new_eqs
        self.num_eqs = len(self.equations)
        self.lead_terms = []
        for eq in self.equations: 
            lt = lead_term(eq)
            if lt: 
                self.lead_terms.append(lt)
        


































# # All
# class GroebnerEqStore:
#     def __init__(self):
#         self.num_vars = 0
#         self.num_eqs = 0

#         self.queues = [[], []]
#         self.equations = []

#         self.known_vars = set()
#         self.solved_vars = {}

#     # TODO: variable naming + readability
#     # TODO: just improve this shit code

#     def syzygy(self,i,j):
#         f = self.equations[i]
#         g = self.equations[j]
#         return (
#             f.poly * BooleanANF([g.lead_term-f.lead_term]) +
#             g.poly * BooleanANF([f.lead_term-g.lead_term]) 
#         )
    
#     # full-reduce
#     # TODO: This needs attention (prob correct but NOT good)
#     def reduce(self, poly):
#         curr_term = lead_term(poly)
#         output = BooleanANF()

#         while curr_term:
#             selected = next(
#                 (eq for eq in self.equations if eq.lead_term <= curr_term), 
#                 None # default to value to return for empty polynomials
#             )

#             if selected == None:
#                 output.terms ^= frozenset([curr_term])
#                 selected = pq_node(curr_term,BooleanANF([curr_term]))
            
#             poly += (
#                 BooleanANF([curr_term - selected.lead_term]) *
#                 selected.poly
#             )
            
#             curr_term = lead_term(poly)
#         return output + poly 
    


#     def enqueue_equation(self, equation, extra_const = 0, identifier=None, translate_ANF = True):
#         poly = XOR(equation,CONST(extra_const)).compose(self.solved_vars)
#         new_equation = BooleanANF.from_BooleanFunction(poly)
#         self.known_vars |= set(poly.idxs_used())

#         if new_equation.terms:# and (new_equation not in self.seen):
#             heapq.heappush(self.queues[1], pq_node(lead_term(new_equation), new_equation))
#             #self.seen.add(new_equation)

#     def consume_queue(self, num=None, verbose = False):
#         # if no num provided, consume the whole queue
#         # (-1) will decrement but never reach 0
#         if num == None:
#             num = -1

#         # main loop
#         reorder = 0
#         while self.queues[0] or (self.queues[1] and num):
#             # pop from the highest priority queue w/ items:
#             for queue_idx in range(len(self.queues)):
#                 if self.queues[queue_idx]:
#                     candidate = heapq.heappop(self.queues[queue_idx]).poly

#                     # count only when pulling from main queue
#                     if queue_idx == 0:
#                         reorder += 1
#                     if queue_idx == 1:
#                         print("REORDER COST: ", reorder)
#                         reorder = 0
#                         num -= 1
                        
#                     break

#             # check equation validity:
#             if candidate.terms == set([frozenset()]):
#                 return False
#                 # raise ValueError("Inconsistent equations")

#             #print("C: ", candidate)
#             # reduce candidate:
#             reduced = self.reduce(candidate)
#             reduced_lead = lead_term(reduced)
#             #print("R: ", reduced)

#             if not reduced.terms:
#                 continue

#             # mark known vars:
#             if len(reduced_lead) == 1:
#                 if reduced.terms == {reduced_lead}:
#                     self.solved_vars[tuple(reduced_lead)[0]] = 0
#                 if reduced.terms == {reduced_lead,frozenset()}:
#                     self.solved_vars[tuple(reduced_lead)[0]] = 1
                
                    
#             # split and reinsert equations (queue 0)
#             node = pq_node(reduced_lead,reduced)
#             insertion_idx = bisect.bisect(self.equations, node)

#             for idx in range(insertion_idx, len(self.equations)):
#                 heapq.heappush(self.queues[0], self.equations[idx])

#             self.equations = self.equations[:insertion_idx]
#             self.equations.append(node)
        
#             # compute/add syzygies:
#             if queue_idx == 1: #VALIDATE THIS LINE
#                 count = 0
#                 for eq_idx in range(len(self.equations) - 1):
#                     # coprime lead terms won't lead to good sysygies
#                     if not (
#                         self.equations[-1].lead_term & 
#                         self.equations[eq_idx].lead_term
#                     ):
#                         continue

#                     s_poly = self.syzygy(-1, eq_idx)
#                     if s_poly.terms: # and s_poly not in self.seen:
#                         heapq.heappush(self.queues[1], pq_node(lead_term(s_poly), s_poly))
#                         #self.seen.add(s_poly)
#                         count += 1

#                 if verbose: print(f"added {count} to queue. Queue length: {len(self.queue)}")
                
#         self.num_eqs = len(self.equations)
#         return True
    
#     def solved_vars2(self):
#         new_vars = {}
#         for v in self.known_vars:
#             reduced = self.reduce(BooleanANF([[v]]))

#             if reduced == BooleanANF([True]):
#                 #self.solved_vars[v] = CONST(1)
#                 new_vars[v] = CONST(1)
#             elif reduced == BooleanANF([]):
#                 #self.solved_vars[v] = CONST(0)
#                 new_vars[v] = CONST(0)
#         return new_vars

