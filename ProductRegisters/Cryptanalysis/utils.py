from ProductRegisters.FeedbackFunctions import *
from ProductRegisters.BooleanLogic import *
import numpy as np



# some utils for making cryptanalysis / automated testing easier
# will be folded into other things as they get nicer

def get_var_map(
    feedback_fn,
    monomial_profile,
    variable_blocks,
    include_variables = True, # include all base variables
    complete_subsets = False, # ensure variable map is closed under subsets
    lexicographic = True      # sort monomials lexicographically
):
    mons_by_len = {}

    for selectors in monomial_profile.get_monomials(complete_subsets=complete_subsets):
        mon = []
        for block_idx in range(len(selectors)):
            for bit_idx in selectors[block_idx]:
                mon.append(variable_blocks[block_idx][bit_idx])
        
        mon = tuple(sorted(mon))
        if len(mon) in mons_by_len:
            mons_by_len[len(mon)].append(mon)
        else:
            mons_by_len[len(mon)] = [mon]

    # replace length 1 segment if necessary:
    if include_variables:
        mons_by_len[1] = [(i,) for i in range(feedback_fn.size)] 

    # cosmetic changes to order:
    list_segments = []
    for length, mon_list in mons_by_len.items():
        if lexicographic:
            sorted_mons = [m for m in mon_list]
            for i in range(length):
                sorted_mons = sorted(sorted_mons, key=lambda x: x[i])
            list_segments.append((length,sorted_mons))
        else:
            list_segments.append((length,mon_list))
    list_segments = sorted(list_segments, key = lambda x: x[0])

    # merging list segments into output maps
    var_idx = 0
    comb_to_idx = {}
    idx_to_comb = {}
    for segment in list_segments:
        for monomial in segment[1]:
            comb_to_idx[monomial] = var_idx
            idx_to_comb[var_idx] = monomial
            var_idx += 1

    return comb_to_idx #, idx_to_comb

from pysat.examples.rc2 import RC2
from pysat.formula import WCNF

from ProductRegisters.BooleanLogic.ANF import ANF_spec_repr

def max_sat(vars,clauses):
    rc2 = RC2(WCNF())
    for clause in clauses:
        rc2.add_clause(clause)
    
    for v in vars:
        rc2.add_clause([v], weight=1)

    model = rc2.compute()
    rc2.delete()
    return model


def ann(f):
    # anf information
    anf = ANF_spec_repr.from_BooleanFunction(f)
    potential_degrees = sorted(list(set([len(x) for x in anf])))[:-1]

    # sat solve for best total annihilator
    vs = [(v+1) for v in f.idxs_used()]
    clauses = [[-(v+1) for v in term] for term in anf]

    rc2 = RC2(WCNF())
    for clause in clauses:
        rc2.add_clause(clause)
    for v in vs:
        rc2.add_clause([v], weight=1)
    assignments = rc2.compute()
    rc2.delete()
    base_negated_vars = [abs(v)-1 for v in max_sat(vs,clauses) if v<0]

    # use total annihilator as a base for data
    best_degree = len(base_negated_vars)
    best_annihilator = AND(*(NOT(VAR(v)) for v in base_negated_vars))
    best_multiple = CONST(0)
    bits = [1]
    
    for target_degree in potential_degrees:
        # construct DNF which preserves at least 1 small term
        preservation_terms = [term for term in anf if len(term) <= target_degree]
        preservation_fns = []
        for term in preservation_terms:
            preservation_fns.append(AND(*(VAR(v) for v in term)))
        preservation_query = OR(CONST(0),*preservation_fns)

        # construct CNF which eliminates all large terms
        elimination_terms = [term for term in anf if len(term) > target_degree]
        elimination_fns = []
        for term in elimination_terms:
            elimination_fns.append(NOT(AND(*(VAR(v) for v in term))))
        elimination_query = AND(CONST(1),*elimination_fns)
        
        # construct vars and nodes for max_sat:
        total_query = AND(preservation_query,elimination_query)
        clauses, node_map, var_map = total_query.tseytin()
        clauses += [(max(node_map[total_query]),)] # assert query is true
        var_weights = [v for v in var_map.values()]

        # print(clauses)
        # for n, l in node_map.items():
        #     print(l, " - ", n.dense_str())
        # break

        #construct degree tower:
        # deg_map = {}
        # for term in nonzero_fns:
        #     d = len(term.args)
        #     if d in deg_map:
        #         deg_map[d] += [term]
        #     else:
        #         deg_map[d] = [term]

        # deg_terms = [(d, AND(*(NOT(term) for term in terms))) for d, terms in deg_map.items()]
        # deg_terms = sorted(deg_terms, key=lambda x: x[0], reverse=True)
        
        # cum_deg_fns = []
        # deg_gaps = [target_degree - deg_terms[0][0]]
        # curr_fn = CONST(1)
        # for i in range(len(deg_terms)):
        #     curr_fn = AND(deg_terms[i][1],curr_fn)
        #     cum_deg_fns.append(curr_fn)

        #     if i>0: deg_gaps.append( deg_terms[i-1][0] - deg_terms[i][0])

        # clauses, node_map, var_map = curr_fn.tseytin(clauses,node_map,var_map)
        # tower_vars = [max(node_map[fn]) for fn in cum_deg_fns]

        # solve sat instance:
        rc2 = RC2(WCNF())
        for clause in clauses:
            rc2.add_clause(clause)
        
        for v in var_weights:
            rc2.add_clause([v], weight=1)
        
        # for v, w in zip(tower_vars,deg_gaps):
        #     rc2.add_clause([v], weight=w)

        assignments = rc2.compute()
        rc2.delete()

        satisfiable = (assignments != None)
        if satisfiable:

            # calculate annihilator from sat assignment
            negated_vars = set([k for k,v in var_map.items() if assignments[v-1]<0])
            annihilator = AND(*(NOT(VAR(v)) for v in negated_vars))

            # calculate low-degree multiple and its degree
            low_multiple = XOR(*(
                AND(*(VAR(v) for v in term),*(NOT(VAR(v)) for v in negated_vars)) 
                for term in anf if not (term & negated_vars)
            ))

            low_multiple_degree = len(max(
                (term for term in low_multiple.args),
                key = lambda x: len(x.args)
            ).args)

            # update data if better degree is achieved:
            if low_multiple_degree < best_degree:
                best_degree = low_multiple_degree
                best_annihilator = annihilator
                best_multiple = low_multiple
                bits = [0,1]

    return (
        best_degree,
        best_annihilator, 
        best_multiple,
        bits
    )


def ann(f):
    # anf information
    anf = ANF_spec_repr.from_BooleanFunction(f)
    potential_degrees = sorted(list(set([len(x) for x in anf])))[:-1]

    # sat solve for best total annihilator
    vs = [(v+1) for v in f.idxs_used()]
    clauses = [[-(v+1) for v in term] for term in anf]

    rc2 = RC2(WCNF())
    for clause in clauses:
        rc2.add_clause(clause)
    for v in vs:
        rc2.add_clause([v], weight=1)
    assignments = rc2.compute()
    rc2.delete()
    base_negated_vars = [abs(v)-1 for v in max_sat(vs,clauses) if v<0]

    # use total annihilator as a base for data
    best_degree = len(base_negated_vars)
    best_annihilator = AND(*(NOT(VAR(v)) for v in base_negated_vars))
    best_multiple = CONST(0)
    bits = [1]
    
    for target_degree in potential_degrees:
        # construct DNF which preserves at least 1 small term
        preservation_terms = [term for term in anf if len(term) <= target_degree]
        preservation_fns = []
        for term in preservation_terms:
            preservation_fns.append(AND(*(VAR(v) for v in term)))
        preservation_query = OR(CONST(0),*preservation_fns)

        # construct CNF which eliminates all large terms
        elimination_terms = [term for term in anf if len(term) > target_degree]
        elimination_fns = []
        for term in elimination_terms:
            elimination_fns.append(NOT(AND(*(VAR(v) for v in term))))
        elimination_query = AND(CONST(1),*elimination_fns)
        
        # construct vars and nodes for max_sat:
        total_query = AND(preservation_query,elimination_query)
        clauses, node_map, var_map = total_query.tseytin()
        clauses += [(max(node_map[total_query]),)] # assert query is true
        var_weights = [v for v in var_map.values()]

        # solve sat instance:
        rc2 = RC2(WCNF())
        for clause in clauses:
            rc2.add_clause(clause)
        
        for v in var_weights:
            rc2.add_clause([v], weight=1)

        assignments = rc2.compute()
        rc2.delete()

        satisfiable = (assignments != None)
        if satisfiable:

            # calculate annihilator from sat assignment
            negated_vars = set([k for k,v in var_map.items() if assignments[v-1]<0])
            annihilator = AND(*(NOT(VAR(v)) for v in negated_vars))

            # calculate low-degree multiple and its degree
            low_multiple = XOR(*(
                AND(*(VAR(v) for v in term),*(NOT(VAR(v)) for v in negated_vars)) 
                for term in anf if not (term & negated_vars)
            ))

            low_multiple_degree = len(max(
                (term for term in low_multiple.args),
                key = lambda x: len(x.args)
            ).args)

            # update data if better degree is achieved:
            if low_multiple_degree < best_degree:
                best_degree = low_multiple_degree
                best_annihilator = annihilator
                best_multiple = low_multiple
                bits = [0,1]

    return (
        best_degree,
        best_annihilator, 
        best_multiple,
        bits
    )


# def ann_equivs(f, target_degree, annihilator_cost):
#     anf = ANF_spec_repr.from_BooleanFunction(f)

#     # construct DNF which preserves at least 1 small term
#     preservation_terms = [term for term in anf if len(term) <= target_degree]
#     preservation_fns = []
#     for term in preservation_terms:
#         preservation_fns.append(AND(*(VAR(v) for v in term)))
#     preservation_query = OR(CONST(0),*preservation_fns)

#     # construct CNF which eliminates all large terms
#     elimination_terms = [term for term in anf if len(term) > target_degree]
#     elimination_fns = []
#     for term in elimination_terms:
#         elimination_fns.append(NOT(AND(*(VAR(v) for v in term))))
#     elimination_query = AND(CONST(1),*elimination_fns)
    
#     # construct vars and nodes for max_sat:
#     total_query = AND(preservation_query,elimination_query)
#     clauses, node_map, var_map = total_query.tseytin()
#     clauses += [(max(node_map[total_query]),)] # assert query is true
#     var_weights = [v for v in var_map.values()]

#     # initialize solver:
#     rc2 = RC2(WCNF())
#     for clause in clauses:
#         rc2.add_clause(clause)
    
#     for v in var_weights:
#         rc2.add_clause([v], weight=1)

    
#     for assignments in rc2.enumerate():
#         cost = rc2.cost
#         if cost != annihilator_cost:
#             break

#         # calculate annihilator from sat assignment
#         negated_vars = set([k for k,v in var_map.items() if assignments[v-1]<0])
#         annihilator = AND(*(NOT(VAR(v)) for v in negated_vars))

#         # calculate low-degree multiple and its degree
#         low_multiple = XOR(*(
#             AND(*(VAR(v) for v in term),*(NOT(VAR(v)) for v in negated_vars)) 
#             for term in anf if not (term & negated_vars)
#         ))

#         low_multiple_degree = len(max(
#             (term for term in low_multiple.args),
#             key = lambda x: len(x.args)
#         ).args)

#         yield annihilator, low_multiple
#     rc2.delete()


def ann2(f):
    # anf information
    anf = ANF_spec_repr.from_BooleanFunction(f)

    # sat solve for best total annihilator
    vs = [(v+1) for v in f.idxs_used()]
    clauses = [[-(v+1) for v in term] for term in anf]

    rc2 = RC2(WCNF())
    for clause in clauses:
        rc2.add_clause(clause)
    for v in vs:
        rc2.add_clause([v], weight=1)
    assignments = rc2.compute()
    rc2.delete()
    base_negated_vars = [abs(v)-1 for v in max_sat(vs,clauses) if v<0]

    # use total annihilator as a base for data
    best_degree = len(base_negated_vars)
    best_annihilator = AND(*(NOT(VAR(v)) for v in base_negated_vars))
    best_multiple = CONST(0)
    bits = [1]
    

    # construct DNF which preserves at least 1 small term
    preservation_fns = []
    for term in anf:
        preservation_fns.append(AND(*(VAR(v) for v in term)))
    preservation_query = OR(CONST(0),*preservation_fns)

    #construct degree tower to cancel out high degree terms:
    deg_map = {}
    for term in preservation_fns:
        d = len(term.args)
        if d in deg_map:
            deg_map[d] += [term]
        else:
            deg_map[d] = [term]

    deg_terms = [(d, AND(*(NOT(term) for term in terms))) for d, terms in deg_map.items()]
    deg_terms = sorted(deg_terms, key=lambda x: x[0], reverse=True)
    
    cum_deg_fns = []
    deg_gaps = [0]
    curr_fn = CONST(1)
    for i in range(len(deg_terms)):
        curr_fn = AND(deg_terms[i][1],curr_fn)
        cum_deg_fns.append(curr_fn)

        if i>0: deg_gaps.append( deg_terms[i-1][0] - deg_terms[i][0])

    # construct sat clauses
    clauses, node_map, var_map = preservation_query.tseytin()
    clauses, node_map, var_map = curr_fn.tseytin(clauses,node_map,var_map)
    clauses += [(max(node_map[preservation_query]),)]
    normal_vars = [v for v in var_map.values()]
    tower_vars = [max(node_map[fn]) for fn in cum_deg_fns]

    # solve sat instance:
    rc2 = RC2(WCNF())
    for clause in clauses:
        rc2.add_clause(clause)
    
    for v in normal_vars:
        rc2.add_clause([v], weight=1)
    
    for v, w in zip(tower_vars,deg_gaps):
        rc2.add_clause([v], weight=w)

    assignments = rc2.compute()
    rc2.delete()

    satisfiable = (assignments != None)
    if satisfiable:

        # calculate annihilator from sat assignment
        negated_vars = set([k for k,v in var_map.items() if assignments[v-1]<0])
        annihilator = AND(*(NOT(VAR(v)) for v in negated_vars))

        # calculate low-degree multiple and its degree
        low_multiple = XOR(*(
            AND(*(VAR(v) for v in term),*(NOT(VAR(v)) for v in negated_vars)) 
            for term in anf if not (term & negated_vars)
        ))

        low_multiple_degree = len(max(
            (term for term in low_multiple.args),
            key = lambda x: len(x.args)
        ).args)

        # update data if better degree is achieved:
        if low_multiple_degree < best_degree:
            best_degree = low_multiple_degree
            best_annihilator = annihilator
            best_multiple = low_multiple
            bits = [0,1]

    return (
        best_degree,
        best_annihilator, 
        best_multiple,
        bits
    )

