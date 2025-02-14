from ProductRegisters.BooleanLogic import BooleanFunction
from ProductRegisters.BooleanLogic.Inputs import VAR
from functools import reduce, cache, wraps
from itertools import product

# unified interface for inverting both bool-like objects and custom objects like
# RootExpressions, functions, monomials, etc. 
def invert(bool_like):
    # if not is already implemented for this type
    if hasattr(bool_like,'__bool__'):
        return not bool_like
    
    # for custom objects like functions or root expressions
    # we can overwrite __invert__() in the desired way
    else:
        return bool_like.__invert__()

class XOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, array, cache):
        return reduce(
            lambda a, b: a ^ b,
            (cache[arg] for arg in self.args)
        )
    def _eval_ANF(self, array, cache):
        return reduce(
            lambda a, b: a ^ b,
            (cache[arg] for arg in self.args)
        )
        
    def generate_c(self):
        return "(" + " ^ ".join(arg.generate_c() for arg in self.args) + ")"
    def _generate_VHDL(self, cache, array_name):
        return "(" + " XOR ".join(cache[arg] for arg in self.args) + ")"
    def _generate_python(self, cache, array_name):
        return "(" + " ^ ".join(cache[arg] for arg in self.args) + ")"
    def generate_tex(self):
        return " \\oplus \\,".join(arg.generate_tex() for arg in self.args)


    def _merge_redundant(self, cache, subfunctions, in_place = False, p = False):
        if len(self.args) == 1:
            return cache[self.args[0]]
        
        # merge nested args:
        new_args = {} # Dicts maintain order
        for arg in self.args:
            if arg in subfunctions:
                if cache[arg] in new_args:
                    new_args[cache[arg]] += 1
                else:
                    new_args[cache[arg]] = 1
                continue

            if type(cache[arg]) == XOR:
                for nested_arg in cache[arg].args:
                    if nested_arg in new_args:
                        new_args[nested_arg] += 1
                    else:
                        new_args[nested_arg] = 1

            elif cache[arg] in new_args:
                new_args[cache[arg]] += 1
            else:
                new_args[cache[arg]] = 1
        new_args = [arg for arg, count in new_args.items() if (count%2)]

        if in_place:
            self.args = list(new_args)
            return self
        else:
            return XOR(*new_args, arg_limit = self.arg_limit)


    def _binarize(self,new_nodes):
        return reduce(
            lambda a, b: XOR(a,b),
            (new_nodes[arg] for arg in self.args)
        )
 
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (-a,-b,-c),
            (a,b,-c),
            (a,-b,c),
            (-a,b,c)
        ]
    
    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 0: # empty gate (e.g. XOR() ) => assert unsat using this gate
            return [(gate_labels[0]),(-gate_labels[0])]
        
        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-gate_labels[0],arg_labels[0]),(gate_labels[0],-arg_labels[0])]
        
        else:
            # initial node uses 2 args
            clauses = self.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])

            # afterward each node uses the previous + the next arg
            for i in range(1,len(gate_labels)):
                clauses += self.tseytin_formula(
                    gate_labels[i-1], arg_labels[i+1], gate_labels[i]
                )
            return clauses






class AND(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, array, cache):
        return reduce(
            lambda a, b: a & b,
            (cache[arg] for arg in self.args)
        )
    def _eval_ANF(self, array, cache):
        return reduce(
            lambda a, b: a & b,
            (cache[arg] for arg in self.args)
        )
        
    
    def generate_c(self):
        return "(" + " & ".join(arg.generate_c() for arg in self.args) + ")"
    def _generate_VHDL(self, cache, array_name):
        return "(" + " AND ".join(cache[arg] for arg in self.args) + ")"
    def generate_tex(self):
        return "".join(arg.generate_tex() for arg in self.args)
    def _generate_python(self, cache, array_name):
        return "(" + " & ".join(cache[arg] for arg in self.args) + ")"


    def _merge_redundant(self, cache, subfunctions, in_place = False, p = False):
        if len(self.args) == 1:
            return cache[self.args[0]]
        
        seen = set()
        new_args = []
        for arg in self.args:
            if arg in subfunctions:
                new_args.append(cache[arg])
                continue

            if type(cache[arg]) == AND:
                for nested_arg in cache[arg].args:
                    if nested_arg not in seen:
                        seen.add(nested_arg)
                        new_args.append(nested_arg)
                
            elif cache[arg] not in seen:
                seen.add(cache[arg])
                new_args.append(cache[arg])

        if in_place:
            self.args = list(new_args)
            return self
        else:
            return AND(*new_args, arg_limit = self.arg_limit)


    def _binarize(self, new_nodes):
        return reduce(
            lambda a, b: AND(a,b),
            (new_nodes[arg] for arg in self.args)
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (-a,-b,c),
            (a,-c),
            (b,-c)
        ]
   
    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(gate_labels[0]),(-gate_labels[0])]
        
        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-gate_labels[0],arg_labels[0]),(gate_labels[0],-arg_labels[0])]
        
        else:
            # initial node uses 2 args
            clauses = self.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])

            # afterward each node uses the previous + the next arg
            for i in range(1,len(gate_labels)):
                clauses += self.tseytin_formula(
                    gate_labels[i-1], arg_labels[i+1], gate_labels[i]
                )
            return clauses
    
    

    


class OR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args
        
    def _eval(self, array, cache):
        return reduce(
            lambda a, b: a | b,
            (cache[arg] for arg in self.args)
        )
    def _eval_ANF(self, array, cache):
        return invert(reduce(
            lambda a, b: a & b,
            (invert(cache[arg]) for arg in self.args)
        ))
    
    def generate_c(self):
        return "(" + " | ".join(arg.generate_c() for arg in self.args) + ")"
    def _generate_VHDL(self, cache, array_name):
        return "(" + " OR ".join(cache[arg] for arg in self.args) + ")"
    def generate_tex(self):
        return " \\vee ".join(arg.generate_tex() for arg in self.args)
    def _generate_python(self, cache, array_name):
        return "(" + " | ".join(cache[arg] for arg in self.args) + ")"


    def _merge_redundant(self, cache, subfunctions, in_place = False, p = False):
        if len(self.args) == 1:
            return cache[self.args[0]]
        
        seen = set()
        new_args = []
        for arg in self.args:
            if arg in subfunctions:
                new_args.append(cache[arg])
                continue
            
            if type(cache[arg]) == OR:
                for nested_arg in cache[arg].args:
                    if nested_arg not in seen:
                        seen.add(nested_arg)
                        new_args.append(nested_arg)
                
            elif cache[arg] not in seen:
                seen.add(cache[arg])
                new_args.append(cache[arg])

        if in_place:
            self.args = list(new_args)
            return self
        else:
            return OR(*new_args, arg_limit = self.arg_limit)


    def _binarize(self, new_nodes):
        return reduce(
            lambda a, b: OR(a,b),
            (new_nodes[arg] for arg in self.args)
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (a,b,-c),
            (-a,c),
            (-b,c)
        ]

    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(gate_labels[0]),(-gate_labels[0])]
        
        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-gate_labels[0],arg_labels[0]),(gate_labels[0],-arg_labels[0])]
        
        else:
            # initial node uses 2 args
            clauses = self.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])

            # afterward each node uses the previous + the next arg
            for i in range(1,len(gate_labels)):
                clauses += self.tseytin_formula(
                    gate_labels[i-1], arg_labels[i+1], gate_labels[i]
                )
            return clauses




class XNOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, array, cache):
        return invert(reduce(
            lambda a, b: a ^ b,
            (cache[arg] for arg in self.args)
        ))
    def _eval_ANF(self, array, cache):
        return invert(reduce(
            lambda a, b: a ^ b,
            (cache[arg] for arg in self.args)
        ))
        
    def generate_c(self):
        return "(!(" + " ^ ".join(arg.generate_c() for arg in self.args) + "))"
    def _generate_VHDL(self, cache, array_name):
        return "(" + " XNOR ".join(cache[arg] for arg in self.args) + ")"
    def _generate_python(self, cache, array_name):
        return "(1-(" + " ^ ".join(cache[arg] for arg in self.args) + "))"


    def _binarize(self, new_nodes):
        return XNOR(
            reduce(
                lambda a, b: XOR(a,b),
                (new_nodes[arg] for arg in self.args[:-1])
            ), 
            new_nodes[self.args[-1]]
        )
 
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (a,b,c),
            (-a,-b,c),
            (-a,b,-c),
            (a,-b,-c)
        ]
    
    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 0: # empty gate (e.g. XOR() ) => assert unsat using this gate
            return [(gate_labels[0]),(-gate_labels[0])]
        
        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-gate_labels[0],arg_labels[0]),(gate_labels[0],-arg_labels[0])]
        
        elif len(arg_labels) == 2: # 2-arg => just use formula
            return self.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])
        
        else:
            # initial node uses 2 args
            clauses = XOR.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])

            # afterward each node uses the previous + the next arg
            # using the associative operation and negating
            for i in range(1,len(gate_labels)-1):
                clauses += XOR.tseytin_formula(
                    gate_labels[i-1], arg_labels[i+1], gate_labels[i]
                )

            # use negation for the final output
            idx = len(gate_labels)-1
            clauses += self.tseytin_formula(
                gate_labels[idx-1], arg_labels[idx+1], gate_labels[idx]
            )
            
            return clauses



class NAND(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, array, cache):
        return invert(reduce(
            lambda a, b: a & b,
            (cache[arg] for arg in self.args)
        ))
    def _eval_ANF(self, array, cache):
        return invert(reduce(
            lambda a, b: a & b,
            (cache[arg] for arg in self.args)
        ))
    
    def generate_c(self):
        return "(!(" + " & ".join(arg.generate_c() for arg in self.args) + "))"
    def _generate_VHDL(self, cache, array_name):
        return "(" + " NAND ".join(cache[arg] for arg in self.args) + ")"
    def _generate_python(self, cache, array_name):
        return "(1-(" + " & ".join(cache[arg] for arg in self.args) + "))"


    def _binarize(self, new_nodes):
        return NAND(
            reduce(
                lambda a, b: AND(a,b),
                (new_nodes[arg] for arg in self.args[:-1])
            ), 
            new_nodes[self.args[-1]]
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (-a,-b,-c),
            (a,c),
            (b,c)
        ]

    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(gate_labels[0]),(-gate_labels[0])]
        
        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-gate_labels[0],arg_labels[0]),(gate_labels[0],-arg_labels[0])]
        
        elif len(arg_labels) == 2: # 2-arg => just use formula
            return self.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])
        
        else:
            # initial node uses 2 args
            clauses = AND.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])

            # afterward each node uses the previous + the next arg
            # using the associative operation and negating
            for i in range(1,len(gate_labels)-1):
                clauses += AND.tseytin_formula(
                    gate_labels[i-1], arg_labels[i+1], gate_labels[i]
                )

            # use negation for the final output
            idx = len(gate_labels)-1
            clauses += self.tseytin_formula(
                gate_labels[idx-1], arg_labels[idx+1], gate_labels[idx]
            )

            return clauses




class NOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, array, cache):
        return invert(reduce(
            lambda a, b: a | b,
            (cache[arg] for arg in self.args)
        ))
    def _eval_ANF(self, array, cache):
        return reduce(
            lambda a, b: a & b,
            (invert(cache[arg]) for arg in self.args)
        )
    
    def generate_c(self):
        return "(!(" + " | ".join(arg.generate_c() for arg in self.args) + "))"
    def _generate_VHDL(self, cache, array_name):
        return "(" + " NOR ".join(cache[arg] for arg in self.args) + ")"
    def _generate_python(self, cache, array_name):
        return "(1-(" + " | ".join(cache[arg] for arg in self.args) + "))"



    def _binarize(self, new_nodes):
        return NOR(
            reduce(
                lambda a, b: OR(a,b),
                (new_nodes[arg] for arg in self.args[:-1])
            ), 
            new_nodes[self.args[-1]]
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (a,b,c),
            (-a,-c),
            (-b,-c)
        ]

    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(gate_labels[0]),(-gate_labels[0])]
        
        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-gate_labels[0],arg_labels[0]),(gate_labels[0],-arg_labels[0])]
        
        elif len(arg_labels) == 2: # 2-arg => just use formula
            return self.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])
        
        else:
            # initial node uses 2 args
            clauses = OR.tseytin_formula(arg_labels[0],arg_labels[1],gate_labels[0])

            # afterward each node uses the previous + the next arg
            # using the associative operation and negating
            for i in range(1,len(gate_labels)-1):
                clauses += OR.tseytin_formula(
                    gate_labels[i-1], arg_labels[i+1], gate_labels[i]
                )

            # use negation for the final output
            idx = len(gate_labels)-1
            clauses += self.tseytin_formula(
                gate_labels[idx-1], arg_labels[idx+1], gate_labels[idx]
            )
            
            return clauses




class NOT(BooleanFunction):
    def __init__(self, *args):
        if len(args) != 1:
            raise ValueError("NOT takes only 1 argument")
        self.arg_limit = 1
        self.args = args

    @classmethod
    def _node_copy(self, fn, child_copies):
        return NOT(child_copies[fn.args[0]])
    
    def _merge_redundant(self, cache, subfunctions, in_place=False, p = False):
        return NOT(cache[self.args[0]])
    
    
    def _eval(self, array, cache):
        return invert(cache[self.args[0]])
    def _eval_ANF(self, array, cache):
        return invert(cache[self.args[0]])

    def generate_c(self):
        return "(!(" + f"{self.args[0].generate_c()}" + "))"
    def _generate_VHDL(self, cache, array_name):
        return "(NOT(" + f"{cache[self.args[0]]}" + "))"
    def _generate_python(self, cache, array_name):
        return "(1-(" + f"{cache[self.args[0]]}" + "))"


    def _binarize(self, new_nodes):
        return NOT(new_nodes[self.args[0]])
    
    @classmethod
    def tseytin_formula(self,a,c):
        return [
            (-a,-c),
            (a,c)
        ]

    @classmethod
    def tseytin_unroll(self,gate_labels,arg_labels):
        if len(arg_labels) == 1: # 1-arg => assert arg and output var are opposite
            return [(-gate_labels[0],-arg_labels[0]),(gate_labels[0],arg_labels[0])]

        if len(arg_labels) == 0: # should never happen, assert unsat
            return [(gate_labels[0]),(-gate_labels[0])]

    
    
