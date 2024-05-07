from ProductRegisters.BooleanLogic.BooleanFunction import BooleanFunction

class CONST(BooleanFunction):
    def __init__(self, value):
        self.value = value

    def is_leaf(self): 
        return True
    def max_idx(self): 
        return -1
    def idxs_used(self):
        return set()


    def eval(self, array):
        return self.value
    def eval_ANF(self, array):
        return self.value


    def generate_c(self):
        return f"{self.value}"
    def generate_VHDL(self):
        return f" '{self.value}' "
    def generate_python(self):
        return f"{self.value}"
    def generate_JSON(self):
        return {
            'class': 'CONST',
            'data': {
                'value': self.value,
            }
        }
    


    # overwriting BooleanFunction
    def pretty_lines(self,depth = 0):
        return [f"CONST({self.value})"]
    def dense_str(self):
        return f"CONST({self.value})"


    # overwriting BooleanFunction
    def remap_constants(self, const_map):
        return CONST(const_map[self.value])
    def remap_indices(self, index_map):
        return CONST(self.value)
    def shift_indices(self, shift_amount):
        return CONST(self.value)


    # overwriting BooleanFunction
    def __copy__(self):
        return type(self)(self.value)

    # overwriting BooleanFunction
    def compose(self, input_map):
        return CONST(self.value)


    # overwriting BooleanFunction
    def component_count(self):
        return {"CONST":1}

    #overwriting BooleanFunction
    def inputs(self):
        return {self}
    
    def _binarize(self):
        return self
    
    def _tseytin_labels(self,node_labels, variable_labels, next_idx):
        if self.value == 0:
            node_labels[self] = [-1]
        elif self.value == 1:
            node_labels[self] = [1]
        else:
            raise ValueError("Bad Constant")
        return next_idx

    def _tseytin_clauses(self, label_map):
        return []

    def num_nodes(self):
        return 1




class VAR(BooleanFunction):
    def __init__(self, index):
        self.index = index

    def is_leaf(self): 
        return True
    def max_idx(self): 
        return self.index
    def idxs_used(self): 
        return set(self.index)


    def eval(self, array):
        return array[self.index]
    def eval_ANF(self, array):
        return array[self.index]
    
    def generate_c(self):
        return f"(*currstate)[{self.index}]"
    def generate_VHDL(self):
        return f"currstate({self.index})"
    def generate_python(self):
        return f"currstate[{self.index}]"
    def generate_JSON(self):
        return {
            'class': 'VAR',
            'data': {
                'index': self.index,
            }
        }


    # overwriting BooleanFunction str methods
    def pretty_lines(self,depth = 0):
        return [f"VAR({self.index})"]
    def dense_str(self):
        return f"VAR({self.index})"
    

    # overwriting BooleanFunction
    def remap_constants(self, constant_map):
        return VAR(self.index)
    def remap_indices(self, index_map):
        return VAR(index_map[self.index])
    def shift_indices(self, shift_amount):
        return VAR(self.index + shift_amount)


     # overwriting BooleanFunction
    def __copy__(self):
        return type(self)(self.index)


    # overwriting BooleanFunction
    def compose(self, input_map):
        return input_map[self.index]


    # overwriting BooleanFunction
    def component_count(self):
        return {"VAR":1}

    # overwriting BooleanFunction
    def inputs(self):
        return {self}



    def _binarize(self):
        return self
    
    def _tseytin_labels(self, node_labels, variable_labels, next_idx):
        if self.index in variable_labels:
            node_labels[self] = [variable_labels[self.index]]
            return next_idx
        
        variable_labels[self.index] = next_idx
        node_labels[self] = [next_idx]
        return next_idx + 1
    
    def _tseytin_clauses(self, label_map):
        return []

    def num_nodes(self):
        return 1