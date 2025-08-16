from typing import Optional, Any
from collections.abc import Iterable

from PyPR.BooleanLogic.BooleanFunction import BooleanFunction
from PyPR.BooleanLogic.Gates import XOR, AND
from PyPR.BooleanLogic.FunctionInputs import CONST, VAR

from itertools import product

# A container class which can hold an BooleanANF of any hashable type
# Note that this class is unordered, because it uses sets. 
class BooleanANF:
    terms: frozenset[frozenset[Any]]

    @classmethod
    def _convert_iterable_term(cls, 
        term: Iterable[Any] | bool | int
    ) -> frozenset[Any] | None:
        """Helper function converts a variety of iterables into a common from (frozensets)

        Any iterable of hashable objects will be converted using a `frozenset()` call. There are
        some checks for specific values which will be interpreted as either logical 1 (represented
        in the BooleanANF class as an empty frozenset) or a logical 0 (represented by returning
        None). These checks enable constructing a BooleanANF in a more flexible and readible way.
        The behaviour of the function is as follows:

        Will return `frozenset()` (logical True) if and only if `term` is:
         - An empty iterable, which converts naturally to an empty frozenset()
         - The bool `True`
         - The int 1
         
        Will return `None` if and only if `term` is:
         - The bool False
         - The int 0

        Anything else will be cast to frozenset, and an error will be raised if that conversion fails
          
        :param term: The term to be converted. Any iterable of hashable objects, as well as
            specific int and bool representations for the constant term.
        :type term: Iterable[Any] | bool | int
        :raises TypeError: if term is passed which does not match the special exceptions, and
            which is not able to be converted to frozenset,
        :return: a frozenset representing the term
        :rtype: frozenset[Any] | None
        """
        # exceptions for ease of use:
        if (term == False or term == 0):
            return None
        if (term == True or term == 1):
            return frozenset()

        try:
            # we can ignore type warning because we plan to catch error.
            return frozenset(term) # type: ignore
        except TypeError:
            raise TypeError( # thow better error.
                f"Unable to convert term {str(term)} (of type {type(term)})  to frozenset."
            )

    def __init__(self, 
        nested_iterable: Any = None, 
        fast_init: bool = False
    ):
        # fast init allows you to skip extra input handling
        if fast_init:
            self.terms = nested_iterable
            return 
        
        #if they want an empty BooleanANF object.
        if nested_iterable == None: 
            self.terms = frozenset()
            return
        
        #convert any nested iterable to proper set/frozenset format:
        terms = set()
        for term in nested_iterable:
            new_term = BooleanANF._convert_iterable_term(term)
            if not(new_term == None):
                terms ^= set([new_term])

        self.terms = frozenset(terms)

    def degree(self) -> int:
        """Compute the algebraic degree of the ANF expression.

        In the ANF representation, algebraic degree can be found by simply taking the
        maximum over the lengths of all terms. 

        :return: The algebraic degree of the function
        :rtype: int
        """
        return max([len(term) for term in self.terms], default=0)

    # BooleanANF Operations:
    # ADD and XOR
    def __xor__(self, other: "BooleanANF") -> "BooleanANF": 
        """Compute the BooleanANF corresponding to the XOR of two BooleanANFs

        This is also equivalent to `__add__(self, other)`, and is computed using
        the exclusive union of the two term frozensets.

        :param other: the other BooleanANF involved in the operation
        :type other: "BooleanANF"
        :return: the XOR of the two BooleanANFs
        :rtype: BooleanANF
        """
        return BooleanANF(self.terms ^ other.terms,fast_init=True)
    def __add__(self, other: "BooleanANF") -> "BooleanANF": 
        """An Alias for `__xor__`, which computes the BooleanANF corresponding to 
        the XOR of two BooleanANFs

        Internally, this is computed using the exclusive union of the two term frozensets.

        :param other: the other BooleanANF involved in the operation
        :type other: "BooleanANF"
        :return: the XOR of the two BooleanANFs
        :rtype: BooleanANF
        """
        return BooleanANF(self.terms ^ other.terms,fast_init=True)
    

    # MUL and AND (can make an BooleanANF very large, use w/ caution)
    def __and__(self,other: "BooleanANF") -> "BooleanANF":
        """Compute the BooleanANF corresponding to the AND of two BooleanANFs

        This is also equivalent to `__mul__(self, other)`. It is computed using the
        frozenset union to merge terms, and fast addition/removal of terms in a regular
        set, which is converted to a frozenset at the end of the multiplication.

        :param other: the other BooleanANF involved in the operation
        :type other: "BooleanANF"
        :return: the AND of the two BooleanANFs
        :rtype: BooleanANF
        """
        termset = set()
        for a,b in product(self.terms,other.terms):
            new_term = a|b
            if new_term in termset:
                termset.remove(new_term)
            else:
                termset.add(new_term)
        return BooleanANF(frozenset(termset),fast_init=True)
    def __mul__(self,other: "BooleanANF") -> "BooleanANF":
        """An Alias for `__and__`, which computes the BooleanANF corresponding to the
        AND of two BooleanANFs

        It is computed using the frozenset union to merge terms, and fast addition/removal
        of terms in a regular set, which is converted to a frozenset at the end of the
        multiplication.

        :param other: the other BooleanANF involved in the operation
        :type other: "BooleanANF"
        :return: the AND of the two BooleanANFs
        :rtype: BooleanANF
        """
        termset = set()
        for a,b in product(self.terms,other.terms):
            new_term = a|b
            if new_term in termset:
                termset.remove(new_term)
            else:
                termset.add(new_term)
        return BooleanANF(frozenset(termset),fast_init=True)
    
    # add an inverter operation
    def __invert__(self) -> "BooleanANF":
        """Invert a BooleanANF by XORing with logical one (empty frozenset)

        :return: The inverted ANF
        :rtype: BooleanANF
        """
        return self ^ BooleanANF(set([frozenset()]),fast_init=True)

    # use a pretty print for __str__, generic object for repr:
    def __str__(self) -> str:
        """Convert a booleanANF into a simplified string representation

        Each term is grouped by parentheses, and comma-separated. Terms with
        only one variable omit the parentheses to save space.

        :return: A string representation of the function
        :rtype: str
        """
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
    def __len__(self) -> int: 
        """Return the number of terms in the BooleanANF

        :return: the number of terms in the BooleanANF
        :rtype: int
        """
        return len(self.terms)
    def __eq__(self,other: "BooleanANF") -> bool: 
        """Determine if two BooleanANFs are equal.

        because frozenset equality depends on the the equality of the stored items (and is
        order independent), we can use the built-in equality check on the terms sets to check 
        equality between BooleanANF objects. Because ANF is a normal form, this also implies
        that the two BooleanANFs are equivalent as functions.

        :param other: the other BooleanANF involved in the operation
        :type other: "BooleanANF"
        :return: whether or not the two functions are equal.
        :rtype: bool
        """
        return self.terms == other.terms
    def __hash__(self) -> int: 
        """Return the hash of the underlying term set

        Because we are using the builtin term set equality, we also use the builtin 
        term set hash, which will always agree with equality. This also makes it 
        possible to store and hash BooleanANFs nicely.

        :return: The hash of the BooleanANF
        :rtype: int
        """
        return self.terms.__hash__()
    def __iter__(self) -> Iterable[frozenset[Any]]: 
        """Expose an iterator to the underlying term set.

        Because frozensets are unordered, there is no guarantee for order in the returned
        iterator. Still, this is convenient because it allows you to loop through terms or
        pass BooleanANFs to `sorted` (or other functions expecting an iterable)
        
        :return: An iterator for the underlying term set.
        :rtype: Iterable[frozenset[Any]]
        """
        return iter(self.terms)
    def __contains__(self, term: Iterable[Any] | int | bool) -> bool: 
        """Determines if a term is present in the ANF

        uses `_convert_iterable_term`, which allows the user to pass iterables 
        other than frozensets, and also use more conveniet shortcuts for the
        logical one (i.e. the bool True or the int 1). If a term corresponding 
        to an empty ANF (i.e.) is passed, True is returned by default.

        :param term: The term to check
        :type term: Iterable[Any] | int | bool
        :return: whether or not the term is contained in the BooleanANF
        :rtype: bool
        """
        converted_term = self._convert_iterable_term(term)
        if converted_term == None:
            return True
        else:
            return (self._convert_iterable_term(term) in self.terms)

    # Conversion methods
    @classmethod
    def from_BooleanFunction(cls,
        fn: BooleanFunction
    ) -> "BooleanANF":
        """Convert a `BooleanFunction` to an equivalent `BooleanANF`

        :param fn: the BooleanFunction to convert.
        :type fn: BooleanFunction
        :return: A BooleanANF which matches the ANF of the input function
        :rtype: BooleanANF
        """

        var_list = {i: BooleanANF([[i]]) for i in fn.idxs_used()}
        new_fn = fn.remap_constants([
            (0, BooleanANF([0])),
            (1, BooleanANF([1]))
        ])
        return new_fn.eval_ANF(var_list)
    
    def to_BooleanFunction(self: "BooleanANF") -> "BooleanFunction":
        """Convert to an equivalent `BooleanFunction`

        Every nonconstant term in the ANF will be converted to an AND (even
        terms with only one variable). Constant terms will be converted to
        CONST inputs. If the BooleanANF is empty, it will be filled with a CONST(0) 
        to prevent a gate with no arguments.

        :return: A BooleanFunction which matches the given BooleanANF
        :rtype: BooleanFunction
        """
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


# after loading, add BooleanANF based methods to all boolean functions
def translate_ANF(self: "BooleanFunction") -> "BooleanFunction":
    """Convert a BooleanFunction to its ANF representation

    This function translates a BooleanFunction to another BooleanFunction, but in
    ANF form. This means that it is top level XOR gate, the arguments of which are
    either CONST gates or an AND of VAR nodes. Note that this is the same as 
    `BooleanANF.from_BooleanFunction(fn).to_BooleanFunction`, but that the output
    not be confused with BooleanFunction - This is a "round-trip" back to BooleanFunctions
    
    :return: A BooleanFunction which models an ANF equation.
    :rtype: BooleanFunction
    """
    return BooleanANF.from_BooleanFunction(self).to_BooleanFunction()
BooleanFunction.translate_ANF = translate_ANF

def from_ANF(cls, nested_iterable: Any) -> "BooleanFunction": 
    """Generates a BooleanFunction from either a BooleanANF or a nested iterable.

    If the passed anf is a BooleanANF, convert it to a BooleanFunction (equivalent
    to `anf.to_BooleanFunction()`). Otherwise, anf is assumed to be a nested iterable,
    which can be converted to BooleanANF via the BooleanANF constructor. Then the
    returned function is equivalent to `BooleanANF(anf).to_BooleanFunction`. This
    enables forming a BooleanFunction from a nested list, or data structure which
    can represent ANF via literals in a convenient way.

    :param anf: The ANF object to be used to create the function. Because there are
        a variety of types available to be parsed by the BooleanANF constructor, the
        type is annotated as `Any` to avoid downstream type errors. 
    :type anf: Any
    :return: The boolean function corresponding to the input ANF.
    :rtype: BooleanFunction
    """
    if type(nested_iterable) == BooleanANF:
        return nested_iterable.to_BooleanFunction()
    return BooleanANF(nested_iterable).to_BooleanFunction()
BooleanFunction.from_ANF = classmethod(from_ANF) # type: ignore

def anf_str(self):
    """Print a BooleanFunction using the ANF string format 
        
    This is identical to `str(BooleanANF.from_BooleanFunction(fn))`, and is mostly
    included to make the code less verbose. A minor benefit is that it makes it slightly
    easier to reproduce some of the really old experiments

    :return: The string for the BooleanANF of the function
    :rtype: str
    """
    return str(BooleanANF.from_BooleanFunction(self))
BooleanFunction.anf_str = anf_str

def degree(self):
    """Calculate the algebraic degree of the function

    This is implemented by computing the ANF, and then taking the max of the
    degree of each monomial in the ANF. This is equivalent to
    `BooleanANF.from_BooleanFunction(fn).degree()`. As with all methods which 
    depend on the ANF, this has the potential to become computationally infeasible 
    (but is usually fine).

    :return: The algebraic degree of the function
    :rtype: int
    """
    return BooleanANF.from_BooleanFunction(self).degree()
BooleanFunction.degree = degree

def monomial_count(self, convert = True):
    """Calculate the number of monomials in the ANF of the function

    This is implemented by computing the ANF, and then counting the number of the
    monomials in the ANF. This is equivalent to 
    `len(BooleanANF.from_BooleanFunction(fn).degree().terms)`, but is more readible
    and far less verbose. As with all methods which depend on the ANF, this has the 
    potential to become computationally infeasible (but is usually fine).

    :return: The number of monomials in the ANF of the function
    :rtype: int
    """
    if convert:
        return len(
            BooleanANF.from_BooleanFunction(self)
        )
    
    return len(self.args)
BooleanFunction.monomial_count = monomial_count

# def anf_optimize(self):
#     # graded list of terms, in frozenset form
#     termlist = sorted(
#         BooleanANF.from_BooleanFunction(self).terms,
#         key = lambda x: len(x)
#     )

#     # map term -> its function
#     # pick divisors to reuse
#     # then sequentially add vars
#     fn_dict = {}
#     for term in termlist:
#         if term == frozenset():
#             fn_dict[term] = CONST(1)
#             continue

#         divisor = max(
#             [t for t in fn_dict.keys() if t <= term], 
#             key = lambda x: len(x),
#             default=None
#         )

#         if divisor == None:
#             fn_dict[term] = AND(*(VAR(x) for x in term))
#         else:
#             fn_dict[term] = AND(
#                 fn_dict[divisor], 
#                 *(VAR(x) for x in term - divisor)
#             )

#     top_node = XOR(*fn_dict.values())
#     if not top_node.args:
#         top_node.add_arguments(CONST(0))
#     return top_node
# BooleanFunction.anf_optimize = anf_optimize