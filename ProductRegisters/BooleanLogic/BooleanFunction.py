from numba import njit
import types
import json


class BooleanFunction:
    def __init__(self):
        self.args = None
        self.arg_limit = None        
    
    @classmethod
    def _node_copy(self, fn, child_copies):
        """Helper function which specifies how to copy a specific node in a 
        boolean function.

        Given a dictionary of child copies, form a new node and return it.
        This is similar to recursion but is used with iteration and memoization
        in `.copy()` and `.__copy__()`. This is also the function to override
        if you want custom behavior (such as copying custom attributes) for
        a child class which extends boolean function.

        Args:
            fn (BooleanFunction): the node you wish to return a copy of.
            child_copies(Dict): A dict containing copies of necessary nodes,
            which are needed to make the new copy.                
        
        Returns:
            output (BooleanFunction): A copy of the node.
        """
        return type(fn)(
            *(child_copies[arg] for arg in fn.args),
            arg_limit = fn.arg_limit
        )

    def __copy__(self):
        """Creates a copy of a BooleanFunction.

        Creates a copy of the boolean function on which it is called.
        the output DAG structure should be identical to the DAG structure
        of the function on which it is called.
        
        Returns:
            output (BooleanFunction): A copy of the boolean function.
        """
        copies = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in copies:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                copies[curr_node] = type(curr_node)._node_copy(curr_node,copies)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                copies[curr_node] = curr_node.__copy__()
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return copies[self]
    
    def copy(self):
        """An alias of `__copy__()` which creates a copy of a BooleanFunction.

        Creates a copy of the boolean function on which it is called.
        the output DAG structure should be identical to the DAG structure
        of the function on which it is called.
        
        Returns:
            output (BooleanFunction): A copy of the boolean function.
        """
        return self.__copy__()

    def functionally_equivalent(self,other):
        """Determines if two functions have the same truth table.

        This function determines if two functions are equivalent,
        not just structurally, but semantically. This is accomplished
        through the tseytin transform and a SAT solver, as the
        problem is NP-complete.

        Args:
            Other (BooleanFunction): the boolean function to compare to

        Returns:
            output (Bool): whether or not the two functions have the 
            same truth table.
        
        """
        raise NotImplementedError # defined in CNF.py


    def add_arguments(self, *remove_args):
        """Adds one or more arguments to a BooleanFunction.

        Args:
            new_args (BooleanFunction): A list of BooleanFunction objects to
            be added as children.

        Raises:
            ValueError: if the number of arguments would make len(self.args)
            greater than self.arg_limit
        """
        if (not self.arg_limit) or (len(self.args) + len(remove_args) <= self.arg_limit):
            self.args = tuple(list(self.args) + list(remove_args))
        else:
            raise ValueError(f"{type(self)} object supports at most {self.arg_limit} arguments")
        
    def remove_arguments(self, *remove_args):
        """Removes one or more arguments from a BooleanFunction.

        Args:
            new_args (BooleanFunction): A list of BooleanFunction objects to
            be added as children.
        """
        if not remove_args:
            self.args = []
        else:
            self.args = tuple(x for x in self.args if x not in remove_args)


    def subfunctions(self):
        """Returns a list of subfunctions

        A subfunction is a BooleanFunction object which is referenced 
        more than one time by BooleanFunctions above it in the DAG. The
        subfunctions are returned in the order they are first encountered
        by a DFS (postorder) traveral. This means the returned array is
        topologically sorted.

        Returns:
            subfunctions (List[BooleanFunction]): A topologically sorted
            list of BooleanFunction which are referenced by multiple parents
        """
        subfuncs = []
        visited = set()
        stack = [self]
        last = None

        # used to sort subfuncs by when they appear:
        order = {}
        idx = 0
        
        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting an already visited gate:
            elif curr_node in visited and not curr_node.is_leaf():
                if curr_node not in subfuncs:
                    subfuncs.append(curr_node)
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
                visited.add(curr_node)
                order[curr_node] = idx
                idx += 1

                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                visited.add(curr_node)
                order[curr_node] = idx
                idx += 1
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return sorted(subfuncs, key = lambda x: order[x])

    def inputs(self):
        """Returns a list of the input nodes (e.g. VAR or CONST)

        The list returns all leaf BooleanFunction objects, and does not
        merge or exclude any objects (even if they are semantically equivalent)
        In other words, two VAR objects might both be included if they are
        distinct objects, even if they reference the same variable.

        Returns:
            inputs (List[BooleanFunction]): A list of all the leaf nodes of
            the DAG, which serve as inputs to the
        """
        leaves = []
        visited = set()
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in visited:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
                visited.add(curr_node)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                leaves.append(curr_node)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return leaves


    def pretty_str(self):
        subfuncs = self.subfunctions() + [self]
        fn_strings = {root:f"(subfunction {i+1})" for i,root in enumerate(subfuncs)}
        fn_strings[self] = f"(Main Function)"
        pretty_strings = []

        for root in subfuncs:
            # no need for visited set (subfuncs is the same info):
            stack = [root]
            last = None

            pstr = f"{fn_strings[root]} = (\n"
            indent_lvl = 0

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                if curr_node in subfuncs and curr_node != root:
                    pstr += (
                        "   |" * indent_lvl + "   " + 
                        fn_strings[curr_node] + '\n'
                    )

                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    pstr += (
                        "   |" * indent_lvl + "   " +
                        curr_node.dense_str() + '\n'
                    )

                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    pstr += ("   |" * (indent_lvl-1) + "   )\n")
                    indent_lvl -= 1
                    last = stack.pop()
                    continue
                
                # before moving down to children:
                else:
                    # print prefix for functions with args:
                    pstr += (
                        "   |" * (indent_lvl) + "   "
                        f"{type(curr_node).__name__} (\n"
                    )

                    indent_lvl += 1
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            # finish string and add it in:
            pstr += ")\n"
            pretty_strings.append(pstr)

        return "\n\n".join(pretty_strings)
    
    def dense_str(self):
        subfuncs = self.subfunctions()
        fn_strings = {
            root:f"(subfunction {i+1})" 
            for i,root in enumerate(subfuncs)
        }
        
        # no need for visited set (subfuncs is the same info)
        stack = [self]
        last = None
        out_str = ""

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a subfunc:
            if curr_node in subfuncs:
                out_str += (fn_strings[curr_node] + ",")
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                out_str += (curr_node.dense_str() + ",")
                last = stack.pop()
                continue
                
            # moving up the tree after finishing children:
            elif last == False:
                # strip trailing comma
                if out_str[-1] == ',':
                    out_str = out_str[:-1]
                
                out_str += "),"
                last = stack.pop()
                continue
            
            # before moving down to children:
            else:
                # print prefix for functions with args:
                out_str += f"{type(curr_node).__name__}("

                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                continue

        # strip trailing comma
        if out_str[-1] == ',':
            out_str = out_str[:-1]
        return out_str


    def generate_c(self):
        pass

    def generate_VHDL(self):
        pass 

    def _generate_VHDL(self):
        pass

    def generate_VHDL(self,
        output_name = 'output',
        subfunction_prefix = 'fn', 
        array_name = 'array',
        overrides = {}
    ):
        subfuncs = self.subfunctions() + [self]
        fn_strings = {
            root:f"{subfunction_prefix}_{i+1}" 
            for i,root in enumerate(subfuncs)
        }

        fn_strings[self] = f"{output_name}"
        subfunction_lines = []
        vhdl_strings = {}

        for root in subfuncs:
            # dont rederive an expression we already have:
            if root in overrides:
                continue

            # no need for visited set (subfuncs is the same info):
            stack = [root]
            last = None

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # override node already has a string
                elif curr_node in overrides:
                    vhdl_strings[curr_node] = overrides[curr_node]
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                elif curr_node in subfuncs and curr_node != root:
                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    vhdl_strings[curr_node] = curr_node._generate_VHDL(vhdl_strings,array_name)
                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    vhdl_strings[curr_node] = curr_node._generate_VHDL(vhdl_strings,array_name)
                    last = stack.pop()
                    continue
            
                else:
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            subfunction_lines.append(f"{fn_strings[root]} <= {vhdl_strings[root]}" + ";")
            vhdl_strings[root] = fn_strings[root]
        return subfunction_lines

    def _generate_python(self):
        pass

    def generate_python(self,
        output_name = 'output',
        subfunction_prefix = 'fn', 
        array_name = 'array',
        overrides = {}
    ):
        subfuncs = self.subfunctions() + [self]
        fn_strings = {
            root:f"{subfunction_prefix}_{i+1}" 
            for i,root in enumerate(subfuncs)
        }

        fn_strings[self] = f"{output_name}"
        subfunction_lines = []
        py_strings = {}

        for root in subfuncs:
            # dont rederive an expression we already have:
            if root in overrides:
                continue


            # no need for visited set (subfuncs is the same info):
            stack = [root]
            last = None

            while stack:
                curr_node = stack[-1]

                # dont interact with sentinel values
                if curr_node == False:
                    last = stack.pop()
                    continue

                # override node already has a string
                elif curr_node in overrides:
                    py_strings[curr_node] = overrides[curr_node]
                    last = stack.pop()
                    continue

                # hitting a subfunc contained in the current one:
                elif curr_node in subfuncs and curr_node != root:
                    last = stack.pop()
                    continue

                # hitting a leaf:
                elif curr_node.is_leaf():
                    py_strings[curr_node] = curr_node._generate_python(py_strings,array_name)
                    last = stack.pop()
                    continue
                    
                # moving up the tree after finishing children:
                elif last == False:
                    py_strings[curr_node] = curr_node._generate_python(py_strings,array_name)
                    last = stack.pop()
                    continue
            
                else:
                    # set up children to process:
                    stack.append(False) # sentinel value
                    for child in reversed(curr_node.args):
                        stack.append(child)
                    continue

            subfunction_lines.append(f"{fn_strings[root]} = {py_strings[root]}")
            py_strings[root] = fn_strings[root]
        return subfunction_lines

    def generate_tex(self):
        pass


    def remap_indices(self, index_map, in_place = False):
        if in_place: fn = self
        else: fn = self.copy()

        for leaf in fn.inputs():
            leaf._remap_indices(index_map)
        return fn
    
    def remap_constants(self, constant_map, in_place = False):
        if in_place: fn = self
        else: fn = self.copy()

        for leaf in fn.inputs():
            leaf._remap_constants(constant_map)
        return fn

    def shift_indices(self, shift_amount, in_place = False):
        if in_place: fn = self
        else: fn = self.copy()

        for leaf in fn.inputs():
            leaf._shift_indices(shift_amount)
        return fn

    def compose(self, input_map, in_place = False):
        new_nodes = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            if curr_node in new_nodes:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
            
                # new node:
                if in_place:
                    curr_node.args = [new_nodes[arg] for arg in curr_node.args]
                    new_nodes[curr_node] = curr_node
                else:
                    new_nodes[curr_node] = type(curr_node)._node_copy(curr_node,new_nodes)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                new_nodes[curr_node] = curr_node._compose(input_map,in_place)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return new_nodes[self]
   
    def _merge_redundant(self, cache, subfunctions, in_place = False):
        if len(self.args) == 1:
            return cache[self.args[0]]
        elif in_place:
            self.args = [cache[arg] for arg in self.args]
            return self
        else:
            return type(self)._node_copy(self, cache)
    
    def merge_redundant(self, in_place = False):
        """Performs some basic simplifications on a BooleanFunction

        Performs several basic/heuristic simplifications. By default, first removes 
        non-unary functions with only one input, as these have no effect 
        on the output. Secondly, merge any associative gates (e.g. XOR, AND, OR),
        unless the child is a subfunction (does not include leaves).
        Even these basic simplifications can result in a dramatically
        simplified function when applied recursively, and help clean up 
        structures created during function composition or other modifications.
        custom simplifications can be added by overriding the `_merge_redundant()`
        helper function.

        Args:
            in_place (Bool): If set to false (by default), the method will return a
            new function, leaving the orignial unmodified (highly recommended usage). 
            If set to true, it will modify the function in place as much as is possible.
            However, because some simplifications change the root node (and thus 
            cannot be done in place), these changes are ommitted. The fully simplified
            function will still be returned, but it is possible for the two to be different
            (i.e. it is possible that `fn = fn.merge_redundant(in_place=True)` is more
            simplified than, and thus not the same as, `fn.merge_redundant(in_place=True)`).

        Returns:
            output (BooleanFunction): A simplified boolean function.
        """
        
        subfunctions = self.subfunctions()

        new_nodes = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            if curr_node in new_nodes:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:
            
                # new node:
                new_nodes[curr_node] = curr_node._merge_redundant(
                    new_nodes, subfunctions, in_place = in_place,
                )

            # hitting a leaf:
            elif curr_node.is_leaf():
                new_nodes[curr_node] = curr_node
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return new_nodes[self]



    def _binarize(self):
        raise NotImplementedError
    
    def binarize(self):
        new_nodes = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in new_nodes:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                new_nodes[curr_node] = curr_node._binarize(new_nodes)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                new_nodes[curr_node] = curr_node.__copy__()
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return new_nodes[self]


    def eval(self, array):
        values = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in values:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                values[curr_node] = curr_node._eval(array, values)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                values[curr_node] = curr_node._eval(array, values)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return values[self]
    
    def eval_ANF(self, array):
        values = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in values:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False:

                # create a deep copy:
                values[curr_node] = curr_node._eval_ANF(array, values)
                last = stack.pop()
                continue

            # hitting a leaf:
            elif curr_node.is_leaf():
                # Overwritten in Inputs.py
                values[curr_node] = curr_node._eval_ANF(array, values)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return values[self]     

    def compile(self):
        self._compiled = None
        python_body = "\n    ".join(self.generate_python())

        exec(f"""
@njit(parallel=True)
def _compiled(array):
    {python_body}
    return output
self._compiled = _compiled
""")
        return self._compiled


    @classmethod
    def construct_ANF(self,nested_iterable):
        raise NotImplementedError  # defined in ANF.py

    def translate_ANF(self):
        raise NotImplementedError  # defined in ANF.py

    def anf_str(self):
        raise NotImplementedError # defined in ANF.py
    
    def degree(self):
        raise NotImplementedError # defined in ANF.py

    def monomial_count(self):
        raise NotImplementedError # defined in ANF.py

   
    def generate_ids(self):
        # generate a unique ID for every node in the tree: 
        next_available_index = 0
        node_labels = {}
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in node_labels:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            # or hitting a leaf
            elif last == False or curr_node.is_leaf():

                # create a deep copy:
                node_labels[curr_node] = next_available_index
                next_available_index += 1
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return node_labels
    
    def _JSON_entry(self,node_ids):
        # copy class name and non-nested data
        JSON_object = {
            'class': type(self).__name__,
            'data': self.__dict__.copy()
        }

        # recurse on any children/nested data:
        if 'args' in JSON_object['data']:
            JSON_object['data']['args'] = [node_ids[arg] for arg in self.args]

        # ignore the compiled version (not serializable)
        if '_compiled' in JSON_object['data']:
            del JSON_object['data']['_compiled']

        return JSON_object
    
    def to_JSON(self):
        node_ids = self.generate_ids()
        num_nodes = max(node_ids.values())+1
        json_node_list = [None for i in range(num_nodes)]

        for node,id in node_ids.items():
            json_node_list[id] = node._JSON_entry(node_ids)
        
        return json_node_list
    
    
    @classmethod
    def from_JSON(self, json_node_list):
        # parse object class and data
        num_nodes = len(json_node_list)
        parsed_functions = [None for i in range(num_nodes)]
        for node_id in range(num_nodes):
            node_data = json_node_list[node_id]
            
            # create information for the python object for this node
            object_class = None
            object_data = node_data['data']

            # self is the BooleanFunction type
            # this finds the appropriate subclass of BooleanFunction for the node
            for subcls in self.__subclasses__():
                if subcls.__name__ == node_data['class']:
                    object_class = subcls

            # throw a better error if no class found
            if object_class == None:
                raise TypeError(f"Type \'{node_data['class']}\' is not a valid BooleanFunction")

            # put data into new object and add it to the parsed functions
            new_node = object.__new__(object_class)
            for key,value in object_data.items():
                if key == 'args':
                    new_node.args = [parsed_functions[child_id] for child_id in value]
                else:
                    setattr(new_node,key,value)
            parsed_functions[node_id] = new_node
        
        # the root node is the last one in the list:
        return parsed_functions[-1]
    
    def to_file(self, filename):
        # .json files only:
        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    @classmethod
    def from_file(self, filename):
        # .json files only:
        with open(filename, 'r') as f:
            return BooleanFunction.from_JSON(json.loads(f.read()))
        

    def is_leaf(self):
        return False
    
    def max_idx(self):
        return max((arg.max_idx() for arg in self.args), default=-1)
    
    def idxs_used(self):
        return set().union(*(arg.idxs_used() for arg in self.args))

    def condense_idxs(self):
        return self.remap_indices(
            {v:i for i,v in enumerate(sorted(self.idxs_used()))}
        )

    def num_nodes(self):
        visited = set()
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in visited:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False or curr_node.is_leaf():
                visited.add(curr_node)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return len(visited)
   
    def component_count(self):
        components = {}
        visited = set()
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # dont interact with sentinel values
            if curr_node == False:
                last = stack.pop()
                continue

            # hitting a visited node while travelling down:
            elif curr_node in visited:
                last=stack.pop()
                continue

            # moving up the tree after finishing children:
            elif last == False or curr_node.is_leaf():
                name = type(curr_node).__name__
                if name in components:
                    components[name] += 1
                else:
                    components[name] = 1
                    
                visited.add(curr_node)
                last = stack.pop()
                continue

            # before moving down to children:
            else:
                # set up children to process:
                stack.append(False) # sentinel value
                for child in reversed(curr_node.args):
                    stack.append(child)
                    continue

        return components
