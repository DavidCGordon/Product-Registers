from ProductRegisters.BooleanLogic.BooleanFunction import BooleanFunction
from pysat.formula import CNF
from pysat.solvers import Solver

def tseytin(self):
    binary_tree = self.binarize()
    node_labels, variable_labels, = binary_tree.tseytin_labels()

    clauses = [[1], node_labels[binary_tree]]
    clauses += binary_tree.tseytin_clauses(node_labels)
    return clauses, node_labels, variable_labels

def tseytin_clauses(self, label_map):
    visited = set()
    stack = [self]
    last = None

    clauses = {}
    while stack:
        curr_node = stack[-1]

        # if current node has no children, or one child has been visited: 
        # you are moving back up the tree
        if curr_node in visited:
            last=stack.pop()

        elif curr_node.is_leaf():
            visited.add(curr_node)
            last = stack.pop()
            
        elif last in curr_node.args:
            try:
                curr_label = label_map[curr_node][0]
                arg_labels = [label_map[arg][0] for arg in curr_node.args]
            except:
                for arg in curr_node.args:
                    print(arg, arg in visited)
                    raise KeyError

            clauses.update(dict.fromkeys(
                type(curr_node).tseytin_formula(*arg_labels, curr_label)
            ))
            #print(clauses)

            visited.add(curr_node)
            last = stack.pop()

        else:
            for child in reversed(curr_node.args):
                stack.append(child)
    #print(clauses)
    return list(clauses.keys())
    

# iterative implementation allows use of a counter in the recursion, 
# rather than passing the next available index through the parameters 
# and return. This is cleaner and easier to reason about in this case
def tseytin_labels(self):
    stack = [self]
    last = None

    next_available_index = 2
    variable_labels = {}
    node_labels = {}
    

    while stack:
        curr_node = stack[-1]

        #don't visit nodes twice:
        if curr_node in node_labels:
            last=stack.pop()

        # handle VAR and CONST Nodes
        # each has it's own implementations
        elif curr_node.is_leaf():
            next_available_index = curr_node._tseytin_labels(
                node_labels,
                variable_labels,
                next_available_index
            )

        # handle gate nodes
        elif last in curr_node.args:
            num_self_labels = max(1,len(self.args)-1)
            node_labels[curr_node] = [next_available_index + i for i in range(num_self_labels)]
            next_available_index += 1
            last = stack.pop()

        # place children in the stack to handle later
        else:
            for child in reversed(curr_node.args):
                stack.append(child)

    return node_labels,variable_labels

def satisfiable(self):
    clauses, node_map, var_map = self.tseytin()
    num_variables = len(node_map)
    num_clauses = len(clauses)
    cnf = CNF(from_clauses=clauses)

    print(cnf.nv, len(cnf.clauses))
    print("Tseytin finished")
    print(f'Number of variables: {num_variables}')
    print(f'Number of clauses: {num_clauses}')

    solver = Solver(
        bootstrap_with=cnf,
        use_timer=True,
    )

    solvable = solver.solve()
    assignments = solver.get_model()

    print(solver.time())
    if solvable:
        return {k: (assignments[v-1]>0) for k,v in var_map.items()}
    else:
        return None

BooleanFunction.tseytin = tseytin
BooleanFunction.tseytin_labels = tseytin_labels
BooleanFunction.tseytin_clauses = tseytin_clauses
BooleanFunction.sat = satisfiable