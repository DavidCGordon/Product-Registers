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

    def eval(self, array):
        return reduce(
            lambda a, b: a ^ b,
            (arg.eval(array) for arg in self.args)
        )
    def eval_ANF(self, array):
        return reduce(
            lambda a, b: a ^ b,
            (arg.eval_ANF(array) for arg in self.args)
        )
        
    def generate_c(self):
        return "(" + " ^ ".join(arg.generate_c() for arg in self.args) + ")"
    def generate_VHDL(self):
        return "(" + " XOR ".join(arg.generate_VHDL() for arg in self.args) + ")"
    def generate_python(self):
        return "(" + " ^ ".join(arg.generate_python() for arg in self.args) + ")"

    @cache
    def _binarize(self):
        return reduce(
            lambda a, b: XOR(a,b),
            (arg._binarize() for arg in self.args)
        )
 
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (-a,-b,-c),
            (a,b,-c),
            (a,-b,c),
            (-a,b,c)
        ]


class AND(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def eval(self, array):
        return reduce(
            lambda a, b: a & b,
            (arg.eval(array) for arg in self.args)
        )
    def eval_ANF(self, array):
        return reduce(
            lambda a, b: a & b,
            (arg.eval_ANF(array) for arg in self.args)
        )
        
    
    def generate_c(self):
        return "(" + " & ".join(arg.generate_c() for arg in self.args) + ")"
    def generate_VHDL(self):
        return "(" + " AND ".join(arg.generate_VHDL() for arg in self.args) + ")"
    def generate_python(self):
        return "(" + " & ".join(arg.generate_python() for arg in self.args) + ")"

    @cache
    def _binarize(self):
        return reduce(
            lambda a, b: AND(a,b),
            (arg._binarize() for arg in self.args)
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (-a,-b,c),
            (a,-c),
            (b,-c)
        ]
   
    
    

    


class OR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args
        
    def eval(self, array):
        return reduce(
            lambda a, b: a | b,
            (arg.eval(array) for arg in self.args)
        )
    def eval_ANF(self, array):
        return invert(reduce(
            lambda a, b: invert(a) & invert(b),
            (arg.eval_ANF(array) for arg in self.args)
        ))
    
    def generate_c(self):
        return "(" + " | ".join(arg.generate_c() for arg in self.args) + ")"
    def generate_VHDL(self):
        return "(" + " OR ".join(arg.generate_VHDL() for arg in self.args) + ")"
    def generate_python(self):
        return "(" + " | ".join(arg.generate_python() for arg in self.args) + ")"

    @cache
    def _binarize(self):
        return reduce(
            lambda a, b: OR(a,b),
            (arg._binarize() for arg in self.args)
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (a,b,-c),
            (-a,c),
            (-b,c)
        ]





class XNOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def eval(self, array):
        return invert(reduce(
            lambda a, b: a ^ b,
            (arg.eval(array) for arg in self.args)
        ))
    def eval_ANF(self, array):
        return invert(reduce(
            lambda a, b: a ^ b,
            (arg.eval_ANF(array) for arg in self.args)
        ))
        
    def generate_c(self):
        return "(!(" + " ^ ".join(arg.generate_c() for arg in self.args) + "))"
    def generate_VHDL(self):
        return "(" + " XNOR ".join(arg.generate_VHDL() for arg in self.args) + ")"
    def generate_python(self):
        return "(1-(" + " ^ ".join(arg.generate_python() for arg in self.args) + "))"

    @cache
    def _binarize(self):
        return XNOR(
            reduce(
                lambda a, b: XOR(a,b),
                (arg._binarize() for arg in self.args[:-1])
            ), 
            self.args[-1]._binarize()
        )
 
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (a,b,c),
            (-a,-b,c),
            (-a,b,-c),
            (a,-b,-c)
        ]
    



class NAND(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def eval(self, array):
        return invert(reduce(
            lambda a, b: a & b,
            (arg.eval(array) for arg in self.args)
        ))
    def eval_ANF(self, array):
        return invert(reduce(
            lambda a, b: a & b,
            (arg.eval_ANF(array) for arg in self.args)
        ))
    
    def generate_c(self):
        return "(!(" + " & ".join(arg.generate_c() for arg in self.args) + "))"
    def generate_VHDL(self):
        return "(" + " NAND ".join(arg.generate_VHDL() for arg in self.args) + ")"
    def generate_python(self):
        return "(1-(" + " & ".join(arg.generate_python() for arg in self.args) + "))"

    @cache
    def _binarize(self):
        return NAND(
            reduce(
                lambda a, b: AND(a,b),
                (arg._binarize() for arg in self.args[:-1])
            ), 
            self.args[-1]._binarize()
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (-a,-b,-c),
            (a,c),
            (b,c)
        ]





class NOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def eval(self, array):
        return invert(reduce(
            lambda a, b: a | b,
            (arg.eval(array) for arg in self.args)
        ))
    def eval_ANF(self, array):
        return reduce(
            lambda a, b: invert(a) & invert(b),
            (arg.eval_ANF(array) for arg in self.args)
        )
    
    def generate_c(self):
        return "(!(" + " | ".join(arg.generate_c() for arg in self.args) + "))"
    def generate_VHDL(self):
        return "(" + " NOR ".join(arg.generate_VHDL() for arg in self.args) + ")"
    def generate_python(self):
        return "(1-(" + " | ".join(arg.generate_python() for arg in self.args) + "))"


    @cache
    def _binarize(self):
        return NOR(
            reduce(
                lambda a, b: OR(a,b),
                (arg._binarize() for arg in self.args[:-1])
            ), 
            self.args[-1]._binarize()
        )
    
    @classmethod
    def tseytin_formula(self,a,b,c):
        return [
            (a,b,c),
            (-a,-c),
            (-b,-c)
        ]






class NOT(BooleanFunction):
    def __init__(self, *args):
        if len(args) != 1:
            raise ValueError("NOT takes only 1 argument")
        self.arg_limit = 1
        self.args = args
    
    def eval(self,array):
        return invert(self.args[0].eval(array))
    def eval_ANF(self,array):
        return invert(self.args[0].eval_ANF(array))

    def generate_c(self):
        return "(!(" + f"{self.args[0].generate_c()}" + "))"
    def generate_VHDL(self):
        return "(NOT(" + f"{self.args[0].generate_VHDL()}" + "))"
    def generate_python(self):
        return "(1-(" + f"{self.args[0].generate_python()}" + "))"

    @cache
    def _binarize(self):
        return NOT(self.args[0]._binarize())
    
    @classmethod
    def tseytin_formula(self,a,c):
        return [
            (-a,-c),
            (a,c)
        ]

            


    
    
