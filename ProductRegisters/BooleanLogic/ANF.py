from ProductRegisters.BooleanLogic.BooleanFunction import BooleanFunction
from ProductRegisters.BooleanLogic.Gates import XOR, AND
from ProductRegisters.BooleanLogic.Inputs import CONST, VAR

from itertools import product

# A container class which can hold an ANF_spec_repr of any hashable type
# Note that this class is unordered, because it uses sets. 
class ANF_spec_repr:
    @classmethod
    def _convert_iterable_term(self, term):
        # accept any iterable of hashable objects for terms
        # Accept True, 1, None, (), [], {} for the "1" element
        # do not accept False or 0.
        try:
            return frozenset(term)
        except TypeError:
            if (term == True or term == 1):
                return frozenset()
        return None

    def __init__(self, nested_iterable = None):
        self.terms = set()

        #if they want an empty ANF_spec_repr object.
        if nested_iterable == None: return

        #convert any nested iterable to proper set/frozenset format:
        for term in nested_iterable:
            new_term = ANF_spec_repr._convert_iterable_term(term)
            if not(new_term == None):
                self.terms ^= set([new_term])
        
    # ANF_spec_repr Operations:
    # ADD and XOR
    def __xor__(self,other): return self.__add__(other)
    def __add__(self,other): return ANF_spec_repr(self.terms ^ other.terms)
    

    # MUL and AND (can make an ANF_spec_repr very large, use w/ caution)
    def __and__(self,other): return self.__mul__(other)
    def __mul__(self,other):
        newANF_spec_repr = ANF_spec_repr()
        for a,b in product(self.terms,other.terms):
            newANF_spec_repr += ANF_spec_repr([a | b])
        return newANF_spec_repr
    
    # add an inverter operation
    def __invert__(self):
        return self ^ ANF_spec_repr([True])

    # use a pretty print for __str__, generic object for repr:
    def __str__(self):
        stringBeginning = ""
        termStrings = []
        for term in self.terms:
            if term == frozenset():
                stringBeginning = "True,"
            elif len(term) == 1:
                termStrings.append(repr(tuple(term)[0]) + ",")
            else:

                #if the internal objects have an order, sort terms.
                try: 
                    termStrings.append("(" + ",".join(repr(t) for t in sorted(term)) + "),")
                except TypeError:
                    termStrings.append("(" + ",".join(repr(t) for t in term) + "),")
        
        #sort string terms, and print vaguely by size 
        return stringBeginning + "".join(sorted(termStrings, key = lambda x: len(x)))[:-1]


    # Generic container methods
    def __len__(self): return len(self.terms)
    def __eq__(self,other): return self.terms == other.terms
    def __hash__(self): return self.terms.__hash__()
    def __iter__(self): return iter(self.terms)
    def add(self,term): self.terms.add(self._convert_iterable_term(term))
    def remove(self,term): self.terms.remove(self._convert_iterable_term(term))
    def __contains__(self,term): return (self._convert_iterable_term(term) in self.terms)

    # Conversion methods

    @classmethod
    def from_BooleanFunction(self,fn):
        var_list = {i: ANF_spec_repr([[i]]) for i in fn.idxs_used()}
        new_fn = fn.remap_constants([
            (0, ANF_spec_repr([0])),
            (1, ANF_spec_repr([1]))
        ])
        return new_fn.eval_ANF(var_list)
    
    def to_BooleanFunction(self):
        top_node = XOR()
        for term in self.terms:
            if type(term) == bool or ((not term) and type(term) != bool):
                new_arg = CONST(1)
            else:
                new_arg = AND(*(VAR(i) for i in term))
            top_node.add_arguments(new_arg)

        # don't return empty XORs:
        if not top_node.args:
            top_node.add_arguments(CONST(0))
            
        return top_node

# after loading, add ANF_spec_repr based methods to all boolean functions
def translate_ANF(fn): return ANF_spec_repr.from_BooleanFunction(fn).to_BooleanFunction()
BooleanFunction.translate_ANF = translate_ANF

def construct_ANF_spec_repr(cls, iterable): return ANF_spec_repr(iterable).to_BooleanFunction()
BooleanFunction.construct_ANF = classmethod(construct_ANF_spec_repr)

def anf_str(self):
    anf = self.translate_ANF()
    stringBeginning = ""
    termStrings = []
    for term in anf.args:
        if type(term) == CONST and term.value == 1:
            stringBeginning = "True,"
        elif type(term) == VAR:
            termStrings.append((str(term.index) + ","))
        else:
            termStrings.append("(" + ",".join(str(var.index) for var in term.args) + "),")
    return stringBeginning + "".join(sorted(termStrings, key = lambda x: len(x)))[:-1]
BooleanFunction.anf_str = anf_str

def degree(self, convert = True):
    if convert:
        return len(max(
            ANF_spec_repr.from_BooleanFunction(self),
            key = lambda x: len(x),
            default = tuple()
        ))
    
    degree = 0
    for term in self.args:
        if type(term) == CONST: degree = max(degree,0)
        elif type(term) == VAR: degree = max(degree,1)
        else:
            degree = max(degree, len(term.args))
        
    return degree
BooleanFunction.degree = degree

def monomial_count(self, convert = True):
    if convert:
        return len(
            ANF_spec_repr.from_BooleanFunction(self)
        )
    
    return len(self.args)
BooleanFunction.monomial_count = monomial_count