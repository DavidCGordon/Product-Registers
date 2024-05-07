from numba import njit
from memoization import cached
import types
import json


import inspect
from functools import cache, wraps


# TODO: recursion on a dense DAG causes a lot of blow up.
# this is fixed in some cases (JSON storage/reconstruction, CNF construction)
# but still causes issues when calling some functions like eval/print
# might need to change the API to restrict to trees / make the DAG structure
# more explicit. 

class BooleanFunction:
    def __init__(self):
        self.args = None
        self.arg_limit = None        

    def __copy__(self):
        return type(self)(
            *(arg.__copy__() for arg in self.args),
            arg_limit = self.arg_limit
        )
        
    def functionally_equivalent(self,other):
        raise NotImplementedError # defined in CNF.py

    def is_leaf(self):
        return False
    
    def max_idx(self):
        return max((arg.max_idx() for arg in self.args), default=-1)
    
    def idxs_used(self):
        return set().union(*(arg.idxs_used() for arg in self.args))

    # add additional arguments to a function, in place.
    def add_arguments(self, *new_args):
        if (not self.arg_limit) or (len(self.args) + len(new_args) <= self.arg_limit):
            self.args = tuple(list(self.args) + list(new_args))
        else:
            raise ValueError(f"{type(self)} object supports at most {self.arg_limit} arguments")
        
    def remove_arguments(self, *rem_args):
        self.args = tuple(x for x in self.args if x not in rem_args)


    def pretty_lines(self):
        start_line = [f"{type(self).__name__} ("]
        arg_lines = []
        end_line = [")"]

        # Recurse on children, and add depth and add to arg_lines
        arg_groups = [arg.pretty_lines() for arg in self.args]
        for group in arg_groups:
            if len(group) > 1:
                arg_lines += (
                    ["   "  + group[0]] + 
                    ["   |" + line for line in group[1:-1]] +
                    ["   "  + group[-1]])
            else:
                arg_lines += ["   "  + group[0]] 

        return start_line + arg_lines + end_line
    

    def pretty_str(self):
        return "\n".join(self.pretty_lines())

    def dense_str(self):
        return f"{type(self).__name__}({','.join(arg.dense_str() for arg in self.args)})"



    def generate_c(self):
        pass

    def generate_VHDL(self):
        pass 

    def generate_python(self):
        pass


    @cached
    def remap_constants(self, const_map):
        return type(self)(*(arg.remap_constants(const_map) for arg in self.args))

    @cached
    def remap_indices(self, index_map):
        return type(self)(*(arg.remap_indices(index_map) for arg in self.args))

    @cached
    def shift_indices(self, shift_amount):
        return type(self)(*(arg.shift_indices(shift_amount) for arg in self.args))

    @cached
    def compose(self, input_map):
        return type(self)(*(arg.compose(input_map) for arg in self.args))

    @cached
    def inputs(self):
        return set().union(*(arg.inputs() for arg in self.args))


    def _binarize(self):
        raise NotImplementedError
    


    def eval(self, array):
        raise NotImplementedError # overwritten per function
    
    def eval_ANF(self, array):
        raise NotImplementedError # overwritten per function        


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


    def component_count(self):
        # get and merge counts from children
        dicts = [arg.component_count() for arg in self.args]
        unified_keys = set.union(*(set(d.keys()) for d in dicts))
        output = {}
        for key in unified_keys:
            output[key] = 0
            for d in dicts:
                # 0 as a default value (if key not in d)
                output[key] += d.get(key, 0)

        # add in this one:
        this_component = type(self).__name__

        if this_component in output:
            output[this_component] += 1
        else:
            output[this_component] = 1
        return output


    def compile(self):
        self._compiled = None

        exec(f"""
@njit(parallel=True)
def _compiled(currstate):
    return {self.generate_python()}
self._compiled = _compiled
""")
        return self._compiled


    def to_JSON(self):
        node_ids = self.generate_ids()
        num_nodes = max(node_ids.values())+1
        json_node_list = [None for i in range(num_nodes)]

        for node,id in node_ids.items():
            json_node_list[id] = node._JSON_entry(node_ids)
        
        return json_node_list

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
        
        

        # throw a better error if no class found
        if object_class == None:
            raise TypeError(f"Type \'{JSON_object['class']}\' is not a valid BooleanFunction")

        
    # json files only:
    def to_file(self, filename):
        with open(filename, 'w') as f:
            f.write(json.dumps(self.to_JSON(), indent = 2))

    # json files only:
    @classmethod
    def from_file(self, filename):
        with open(filename, 'r') as f:
            return BooleanFunction.from_JSON(json.loads(f.read()))

    @cached
    def num_nodes(self):
        return sum(arg.num_nodes() for arg in self.args)
    
    def sample_template(self):
        collected_outputs = []

        for source in self.args:
            for output in source.sample_template():
                if output != None:
                    collected_outputs.append(output)
        
        yield type(self)(*collected_outputs)

    # generate a unique ID for every node in the tree:
    def generate_ids(self):
        stack = [self]
        next_available_index = 0
        node_labels = {}

        while stack:
            curr_node = stack[-1]

            # don't visit nodes twice:
            if curr_node in node_labels:
                stack.pop()

            # add this node when all children have been added (moving up tree)
            elif curr_node.is_leaf() or all([arg in node_labels for arg in curr_node.args]):
                node_labels[curr_node] = next_available_index
                next_available_index += 1
                stack.pop()

            # otherwise place children in the stack to handle later
            else:
                for child in reversed(curr_node.args):
                    stack.append(child)

        return node_labels



# TODO: make this fit in better with the rest of the library, rather than being a 1-off
def iterative_recursion(fn):
    def traversal(self, *args, **kwargs):
        visited = set()
        stack = [self]
        last = None

        while stack:
            curr_node = stack[-1]

            # if current node has no children, or one child has been visited: 
            # you are moving back up the tree
            if curr_node in visited:
                last=stack.pop()

            elif curr_node.is_leaf() or (last in curr_node.args):
                getattr(curr_node, f"{fn.__name__}")(*args, **kwargs)
                visited.add(curr_node)
                last = stack.pop()

            else:
                for child in reversed(curr_node.args):
                    stack.append(child)

        return getattr(last, f"{fn.__name__}")(*args, **kwargs)
    return traversal


BooleanFunction.binarize = iterative_recursion(BooleanFunction._binarize)
