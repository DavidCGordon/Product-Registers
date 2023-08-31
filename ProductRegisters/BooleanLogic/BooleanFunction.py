from numba import njit
from memoization import cached
import types
import json


import inspect
from functools import cache, wraps


class BooleanFunction:
    def __init__(self):
        self.args = None
        self.arg_limit = None        

    @cached
    def __copy__(self):
        return type(self)(
            *(arg.__copy__() for arg in self.args),
            arg_limit = self.arg_limit
        )
        
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



    #TODO: JSON methods break DAG into trees
    def to_JSON(self):
        # copy class name and non-nested data
        JSON_object = {
            'class': type(self).__name__,
            'data': self.__dict__.copy()
        }

        # recurse on any children/nested data:
        if 'args' in JSON_object['data']:
            JSON_object['data']['args'] = [arg.to_JSON() for arg in self.args]

        # ignore the compiled version (not serializable)
        if '_compiled' in JSON_object['data']:
            del JSON_object['data']['_compiled']

        return JSON_object
    
    @classmethod
    def from_JSON(self, JSON_object):
        # parse object class and data
        object_data = JSON_object['data']
        object_class = None
        for subcls in self.__subclasses__():
            if subcls.__name__ == JSON_object['class']:
                object_class = subcls

        # throw a better error if no class found
        if object_class == None:
            raise TypeError(f"Type \'{JSON_object['class']}\' is not a valid BooleanFunction")

        # put data into new object
        output = object.__new__(object_class)
        for key,value in object_data.items():
            if key == 'args':
                output.args = [BooleanFunction.from_JSON(child) for child in value]
            else:
                setattr(output,key,value)
            
        return output

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
