from ProductRegisters.BooleanLogic import BooleanFunction,XOR,AND,OR,NOT,CONST,VAR
from ProductRegisters.BooleanLogic.ANF import ANF_spec_repr
from pysat.formula import CNF,WCNF
from pysat.solvers import Solver
from pysat.examples.rc2 import RC2



def tseytin(self, prev_clauses=None, prev_node_labels=None, prev_variable_labels=None):
    if not prev_clauses: clauses = dict.fromkeys([(1,)])
    else: clauses = dict.fromkeys([tuple(x) for x in prev_clauses])

    node_labels, variable_labels, = self.tseytin_labels(prev_node_labels, prev_variable_labels)
    clauses.update(dict.fromkeys(self.tseytin_clauses(node_labels)))
    return list(clauses.keys()), node_labels, variable_labels

def tseytin_clauses(self, label_map):
    visited = set()
    stack = [self]

    clauses = {}
    while stack:
        curr_node = stack[-1]

        # if current node has no children, or one child has been visited: 
        # you are moving back up the tree
        if curr_node in visited:
            stack.pop()

        elif curr_node.is_leaf():
            visited.add(curr_node)
            stack.pop()
            
        elif all([arg in visited for arg in curr_node.args]):
            curr_labels = label_map[curr_node]
            arg_labels = [label_map[arg][-1] for arg in curr_node.args]
            clauses.update(dict.fromkeys(
                type(curr_node).tseytin_unroll(curr_labels,arg_labels)
            ))

            visited.add(curr_node)
            stack.pop()

        else:
            for child in reversed(curr_node.args):
                stack.append(child)
    return list(clauses.keys())
    

# iterative implementation allows use of a counter in the recursion, 
# rather than passing the next available index through the parameters 
# and return. This is cleaner and easier to reason about in this case
def tseytin_labels(self,node_labels=None,variable_labels=None):
    stack = [self]

    # initialize index maps if needed
    if node_labels == None and variable_labels == None:
        next_available_index = 2
        variable_labels = {}
        node_labels = {}
    
    # if only 1 is passed in, raise an error
    elif node_labels == None:
        raise ValueError("Missing nodelabels")
    elif variable_labels == None:
        raise ValueError("Missing variable labels")

    # if both passed in, justn set the next index
    else:
        next_available_index = max([max(ls) for ls in node_labels.values()]) + 1

    while stack:
        curr_node = stack[-1]

        #don't visit nodes twice:
        if curr_node in node_labels:
            stack.pop()

        # handle VAR and CONST Nodes
        # each has it's own implementation in _tseytin_labels
        elif curr_node.is_leaf():
            next_available_index = curr_node._tseytin_labels(
                node_labels,
                variable_labels,
                next_available_index
            )
            stack.pop()

        # handle gate nodes
        elif all([arg in node_labels for arg in curr_node.args]):
            num_gate_labels = max(1,len(curr_node.args)-1)
            node_labels[curr_node] = [next_available_index + i for i in range(num_gate_labels)]
            next_available_index += num_gate_labels
            stack.pop()

        # place children in the stack to handle later
        else:
            for child in reversed(curr_node.args):
                stack.append(child)

    return node_labels,variable_labels





def satisfiable(self, verbose = False, solver_name = "cadical195", time_limit = None):
    clauses, node_map, var_map = self.tseytin()
    clauses += [(node_map[self][-1],)]
    num_variables = node_map[self][-1] + 1
    num_clauses = len(clauses)
    
    if verbose:
        print("Tseytin finished")
        print(f'Number of variables: {num_variables}')
        print(f'Number of clauses: {num_clauses}')

    cnf = CNF(from_clauses=clauses)
    with Solver(name = solver_name, bootstrap_with=cnf, use_timer=True) as solver:
        satisfiable = solver.solve()
        assignments = solver.get_model()

    if verbose:
        print(solver.time())

    if satisfiable:
        return {k: (assignments[v-1]>0) for k,v in var_map.items()}
    else:
        return None

    
def enumerate_models(self, solver_name = 'cadical195', verbose = False):
    clauses, node_map, var_map = self.tseytin()
    clauses += [(node_map[self][-1],)]
    num_variables = len(node_map)
    num_clauses = len(clauses)
    cnf = CNF(from_clauses=clauses)

    if verbose:
        print(cnf.nv, len(cnf.clauses))
        print("Tseytin finished")
        print(f'Number of variables: {num_variables}')
        print(f'Number of clauses: {num_clauses}')

    with Solver(name = solver_name, bootstrap_with=cnf, use_timer=True) as solver:
        for assignment in solver.enum_models():
            yield {k: (assignment[v-1]>0) for k,v in var_map.items()}













def monomial_annihilator(self):
    # construct max-sat instance
    anf = ANF_spec_repr.from_BooleanFunction(self)
    vs = [(v+1) for v in self.idxs_used()]
    clauses = [[-(v+1) for v in term] for term in anf]

    rc2 = RC2(WCNF())
    for clause in clauses:
        rc2.add_clause(clause)
    for v in vs:
        rc2.add_clause([v], weight=1)
    assignments = rc2.compute()
    rc2.delete()
    base_negated_vars = [abs(v)-1 for v in assignments if v<0]

    # derive data from it
    degree = len(base_negated_vars)
    annihilator = AND(*(NOT(VAR(v)) for v in base_negated_vars))
    multiple = CONST(0)

    return (
        degree,
        annihilator,
        multiple
    )



def low_degree_multiple(self):
    # anf information
    anf = ANF_spec_repr.from_BooleanFunction(self)
    potential_degrees = sorted(list(set([len(x) for x in anf])))[:-1]

    # use best total annihilator as a base
    best_degree,best_annihilator,best_multiple = self.monomial_annihilator
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
    )


def low_degree_multiple2(self):
    # anf information
    anf = ANF_spec_repr.from_BooleanFunction(self)

    # use best total annihilator as a base
    best_degree,best_annihilator,best_multiple = self.monomial_annihilator
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





BooleanFunction.tseytin = tseytin
BooleanFunction.tseytin_labels = tseytin_labels
BooleanFunction.tseytin_clauses = tseytin_clauses
BooleanFunction.sat = satisfiable
BooleanFunction.enum_models = enumerate_models
BooleanFunction.low_degree_multiple = low_degree_multiple

def functionally_equivalent(self,other):
    if (XOR(self,other).sat()) == None: 
        return True
    else:
        return False

BooleanFunction.functionally_equivalent = functionally_equivalent