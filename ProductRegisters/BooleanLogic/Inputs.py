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


    def _eval(self, array, cache):
        return self.value
    def _eval_ANF(self, array, cache):
        return self.value


    def generate_c(self):
        return f"{self.value}"
    def _generate_VHDL(self, cache, array_name):
        return f" '{self.value}' "
    def _generate_python(self, cache, array_name):
        return f"{self.value}"
    def generate_tex(self):
        return f"{self.value}"
    def generate_JSON(self):
        return {
            'class': 'CONST',
            'data': {
                'value': self.value,
            }
        }
    


    # overwriting BooleanFunction
    def dense_str(self):
        return f"CONST({self.value})"


    # overwriting BooleanFunction
    def _remap_constants(self, const_map):
        for key,new_value in const_map:
            try:
                if self.value == key:
                    self.value = new_value
                    break
            except:
                pass
        
    def _remap_indices(self, index_map):
        pass
    def _shift_indices(self, shift_amount):
        pass


    # overwriting BooleanFunction
    def __copy__(self):
        return type(self)(self.value)

    # overwriting BooleanFunction
    def _compose(self, input_map, in_place = False):
        if in_place:
            return self
        else:
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
        return set([self.index])


    def _eval(self, array, cache):
        return array[self.index]
    def _eval_ANF(self, array, cache):
        return array[self.index]
    
    def generate_c(self):
        return f"(*currstate)[{self.index}]"
    def _generate_VHDL(self, cache, array_name):
        return f"{array_name}({self.index})"
    def _generate_python(self, cache, array_name):
        return f"{array_name}[{self.index}]"
    def generate_tex(self):
        return f"c_{{{self.index}}}[t]"
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
    def _remap_constants(self, constant_map):
        pass
    def _remap_indices(self, index_map):
        if self.index in index_map:
            self.index = index_map[self.index]
    def _shift_indices(self, shift_amount):
        self.index = self.index + shift_amount


     # overwriting BooleanFunction
    def __copy__(self):
        return VAR(self.index)


    # overwriting BooleanFunction
    def _compose(self, input_map, in_place = False):
        try:
            return input_map[self.index]
        except:
            if in_place:
                return self
            else:
                return VAR(self.index)


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