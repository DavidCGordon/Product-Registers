from collections.abc import Iterator
from typing import Any

from PyPR.BooleanLogic import BooleanFunction,XOR
from pysat.formula import CNF
from pysat.solvers import Solver


def tseytin(self, 
    prev_clauses: list[tuple[int]] | None = None, 
    prev_node_labels: dict[BooleanFunction,list[int]] | None = None, 
    prev_variable_labels: dict[int,int] | None = None,
) -> tuple[
    list[tuple[int]],
    dict[BooleanFunction,list[int]],
    dict[int,int]
]:
    """Generate constraints corresponding to the tseytin transformation of the function.

    This function primarily deals with three objects:
    - **Clauses**: A list of constraints, suitable for a sat solver. Each constraint is a tuple
        of integers, which represents an OR clause which must be satisfied (CNF). Negative integers
        represent the negation of a variable.
    - **Node Labels**: The Tseytin transform introduces new variables for each gate in the function.
        This dict maps every node to the list of variables which are used to encode it, and can be used 
        to convert the sat solution back to values in the wires of the circuit, or to impose additional
        conditions based on extra information.
    - **Variable Labels**: Because the input variables are always supposed to be the same (even
        in different VAR `nodes`), we use this dict to map input variables to the corresponding
        variables in the sat solver. This can be used to convert the sat solution back to a satisfying
        assignment of input variables, or to impose additional conditions based on extra information.

    When called with no inputs, this function will create these three from scratch. However,
    the function can also recieve the output of previous calls as input, and will continue to build
    out the constraints, reusing variables where possible. This means the same function can also be
    used to add onto existing output. For example:
    ```python
    tseytin_data = function_1.tseytin()
    tseytin_data = function_2.tseytin(*tseytin_data)
    tseytin_data = function_3.tseytin(*tseytin_data)
    ...
    ```
    Or equivalently called in a for loop will create a set of constraints which describes all functions
    and which reuses variables where appropriate. This updating lets you create complicated queries,
    if you need to interact with the sat solver more deeply than just calling `fn.sat()`

    
    :param prev_clauses: A list of constraints, suitable for a sat solver. Each constraint is a tuple
        of integers, which represents an OR clause which must be satisfied (CNF). Negative integers
        represent the negation of a variable.
    :type prev_clauses: list[tuple[int]]
    :param prev_node_labels: The Tseytin transform introduces new variables for each gate in the function.
        This dict maps every node to the list of variables which are used to encode it, and can be used to
        convert the sat solution back to values in the wires of the circuit, or to impose additional
        conditions based on extra information.
    :type prev_node_labels: dict[BooleanFunction,list[int]]
    :param prev_variable_labels: Because the input variables are always supposed to be the same (even
        in different VAR `nodes`), we use this dict to map input variables to the corresponding
        variables in the sat solver. This can be used to convert the sat solution back to a satisfying
        assignment of input variables, or to impose additional conditions based on extra information.
    :type prev_variable_labels: dict[int,int]
    :return: (Clauses, Node Labels, Variable Labels)
    :rtype: tuple[
        list[tuple[int]],
        dict[BooleanFunction,list[int]],
        dict[int,int]
    ]
    """ 
    clauses: dict[tuple,None]
    if not prev_clauses: clauses = dict.fromkeys([(1,)])
    else: clauses = dict.fromkeys([tuple(x) for x in prev_clauses])

    node_labels, variable_labels, = self.tseytin_labels(prev_node_labels, prev_variable_labels)
    clauses.update(dict.fromkeys(self.tseytin_clauses(node_labels)))
    return list(clauses.keys()), node_labels, variable_labels
    
def tseytin_labels(self,
    node_labels: dict[BooleanFunction,list[int]] | None = None,
    variable_labels: dict[int,int] | None = None
) -> tuple[
    dict[BooleanFunction,list[int]],
    dict[int,int]
]:
    """Generate the labels for the tseytin transformation of the function.

    Tseytin primarily deals with three objects: Clauses, Node Labels and variable labels. This
    function is responsible for generating the latter two, which are then used to create clauses.
    In more detail, this  method produces:
    - **Node Labels**: The Tseytin transform introduces new variables for each gate in the function.
        This dict maps every node to the list of variables which are used to encode it, and can be used 
        to convert the sat solution back to values in the wires of the circuit, or to impose additional
        conditions based on extra information.
    - **Variable Labels**: Because the input variables are always supposed to be the same (even
        in different VAR `nodes`), we use this dict to map input variables to the corresponding
        variables in the sat solver. This can be used to convert the sat solution back to a satisfying
        assignment of input variables, or to impose additional conditions based on extra information.

    When called with no inputs, this function will create these two from scratch. However,
    the function can also recieve the output of previous calls as input, and will continue to build
    out the constraints, reusing variables where possible. This means the same function can also be
    used to add onto existing output. For example:
    ```python
    labels = function_1.tseytin_labels()
    labels = function_2.tseytin_labels(*labels)
    labels = function_3.tseytin_labels(*labels)
    ...
    ```
    Or equivalently called in a for loop will create a set of labels which covers all functions
    and which reuses variables where appropriate. This is an important part of the `tseytin` function.
    
    :param prev_node_labels: The Tseytin transform introduces new variables for each gate in the function.
        This dict maps every node to the list of variables which are used to encode it, and can be used to
        convert the sat solution back to values in the wires of the circuit, or to impose additional
        conditions based on extra information.
    :type prev_node_labels: dict[BooleanFunction,list[int]]
    :param prev_variable_labels: Because the input variables are always supposed to be the same (even
        in different VAR `nodes`), we use this dict to map input variables to the corresponding
        variables in the sat solver. This can be used to convert the sat solution back to a satisfying
        assignment of input variables, or to impose additional conditions based on extra information.
    :type prev_variable_labels: dict[int,int]
    :return: (Node Labels, Variable Labels)
    :rtype: tuple[
        dict[BooleanFunction,list[int]],
        dict[int,int]
    ]
    """
    stack = [self]

    # initialize index maps if needed
    if node_labels == None and variable_labels == None:
        next_available_index = 2
        variable_labels = {}
        node_labels = {}
    
    # if only 1 is passed in, raise an error
    elif node_labels == None:
        raise ValueError("Missing node labels")
    elif variable_labels == None:
        raise ValueError("Missing variable labels")
    else:
        # if both passed in, just set the next index
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

def tseytin_clauses(self, 
    label_map: dict[BooleanFunction, list[int]]
) -> list[tuple[int]]:
    """Generate constraints corresponding to the tseytin transformation of the function.

    `Tseytin` primarily deals with three objects: Clauses, Node Labels and variable labels. This
    function is responsible for generating the clauses, and requires the node labels. The Clauses 
    returned are list of constraints, suitable for a sat solver. Each constraint is a tuple
    of integers, which represents an OR clause which must be satisfied (CNF). Negative integers
    represent the negation of a variable.

    :param node_labels: The Tseytin transform introduces new variables for each gate in the function.
        This dict maps every node to the list of variables which are used to encode it, and can be used to
        convert the sat solution back to values in the wires of the circuit, or to impose additional
        conditions based on extra information.
    :type node_labels: dict[BooleanFunction,list[int]]
    :return: clauses which encode the given function for a sat solver.
    :rtype: list[tuple[int]],
    """
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

            # Why did I need to use dict?
            clauses.update(dict.fromkeys(
                type(curr_node).tseytin_unroll(curr_labels,arg_labels)
            ))

            visited.add(curr_node)
            stack.pop()

        else:
            for child in reversed(curr_node.args):
                stack.append(child)
    return list(clauses.keys())

def satisfiable(self,
    solver_name: str = "cadical195", 
    verbose: bool = False
) -> dict[int,bool] | None:
    """Solve the SAT problem for a given BooleanFunction

    First, the `tseytin` method is used to build a SAT encoding of the function.
    Then, the output is manually asserted to be true, and the resulting clauses are
    used as input for a SAT solver provided by the PySAT package. If more complicated
    sat-based procedures are needed, `tseytin` is exposed, and more complicated instances
    can be created manually. However, because direct SAT solving is the most common
    application, this function provides a much simpler interface, abstracting away
    details of the encoding. 

    :param solver_name: A string giving the name of a sat solver provided by PySAT
        which will be used as the solver, defaults to "cadical195", which we observed
        to work well experimentally.
    :type solver_name: str, optional
    :param verbose: if True, print statistics and timings for debugging, defaults to False
    :type verbose: bool, optional

    :return: if the function is unsatisfiable, return None. otherwise, returns a dictionary
        which maps variables to their boolean values in a satisfying assignment. Any variables 
        which don't appear in the dict are "don't care" .
    :rtype: dict[int,bool] | None
    """
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
        assignments: Any = solver.get_model()

    if verbose:
        print(solver.time())

    if satisfiable:
        return {k: (assignments[v-1]>0) for k,v in var_map.items()}
    else:
        return None
  
def enumerate_models(self,
    solver_name: str = 'cadical195',
    verbose: bool = False
) -> Iterator[dict[int,bool]]:
    """Enumerate solutions to the SAT problem for a given BooleanFunction

    First, the `tseytin` method is used to build a SAT encoding of the function.
    Then, the output is manually asserted to be true, and the resulting clauses are
    used as input for a SAT solver provided by the PySAT package. These solvers include
    the ability to enumerate solutions, and this method lifts that to match the PyPR
    interface for sat. 
    
    If more complicated sat-based procedures are needed, `tseytin` is exposed, and 
    more complicated instances can be created manually. However, because direct SAT
    solving and model enumeration are the most common applications, this function provides
    a much simpler interface, abstracting away details of the encoding. 

    :param solver_name: A string giving the name of a sat solver provided by PySAT
        which will be used as the solver, defaults to "cadical195", which we observed
        to work well experimentally.
    :type solver_name: str, optional
    :param verbose: if `True` print statistics and timings for debugging, defaults to False
    :type verbose: bool, optional

    :return: if the function is unsatisfiable, return None. otherwise, on each iteration,
        return a dictionary which maps variables to their boolean values in a satisfying 
        assignment. Any variables which don't appear in the dict are "don't care".
    :rtype: dict[int,bool] | None
    """
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
        for assignment in solver.enum_models(): # type: ignore (this is from bad typing in pysat)
            yield {k: (assignment[v-1]>0) for k,v in var_map.items()}

def functionally_equivalent(self, 
    other: BooleanFunction 
) -> bool:
    """Enumerate solutions to the SAT problem for a given BooleanFunction

    First, the `tseytin` method is used to build a SAT encoding of the function.
    Then, the output is manually asserted to be true, and the resulting clauses are
    used as input for a SAT solver provided by the PySAT package. These solvers include
    the ability to enumerate solutions, and this method lifts that to match the PyPR
    interface for sat. 
    
    If more complicated sat-based procedures are needed, `tseytin` is exposed, and 
    more complicated instances can be created manually. However, because direct SAT
    solving and model enumeration are the most common applications, this function provides
    a much simpler interface, abstracting away details of the encoding. 

    :param solver_name: A string giving the name of a sat solver provided by PySAT
        which will be used as the solver, defaults to "cadical195", which we observed
        to work well experimentally.
    :type solver_name: str, optional
    :param verbose: if `True` print statistics and timings for debugging, defaults to False
    :type verbose: bool, optional

    :return: if the function is unsatisfiable, return None. otherwise, on each iteration,
        return a dictionary which maps variables to their boolean values in a satisfying 
        assignment. Any variables which don't appear in the dict are "don't care".
    :rtype: dict[int,bool] | None
    """
    return ((satisfiable(XOR(self,other))) == None)

# add functions to BooleanFunction class
BooleanFunction.tseytin = tseytin
BooleanFunction.tseytin_labels = tseytin_labels
BooleanFunction.tseytin_clauses = tseytin_clauses
BooleanFunction.sat = satisfiable
BooleanFunction.enum_models = enumerate_models
BooleanFunction.functionally_equivalent = functionally_equivalent