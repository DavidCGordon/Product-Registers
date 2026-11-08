from functools import reduce

from PyPR.BooleanLogic import BooleanFunction
from PyPR.BooleanLogic.FunctionInputs import CONST


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


def _merge_associative(node, gate, cache, subfunctions, in_place, cancel):
    """`_merge_redundant` for an associative gate: one `gate` over the terms of every
    same-gate argument, nested to any depth.

    Terms are an insertion-ordered dict used as a set. With `cancel` (XOR) a term
    occurring twice is dropped, since x ^ x = 0; otherwise (AND, OR) it is kept
    once, since x & x = x | x = x. A subfunction is a single term, never flattened,
    because its other parents still need it whole.

    The result's arguments are not built here but deferred to the end of the walk
    (`cache.pending`), and that is what makes a nested chain linear rather than
    quadratic. A same-gate argument that is not a subfunction has this node as its
    only parent, so its result is discarded once absorbed: rather than copying its
    terms, this node takes them over, folding the smaller set into the larger.
    Building every intermediate gate in full instead copied the whole chain
    beneath it at each level. The order of the merged terms is unspecified.

    When merging leaves a single term the term itself is returned, and an XOR
    left with none returns `CONST(0)`.
    """
    if len(node.args) == 1:
        return cache[node.args[0]]

    pending = cache.pending
    terms = {}
    for arg in node.args:
        result = cache[arg]

        # a single term: a different gate, a leaf, or a subfunction
        if type(result) is not gate or arg in subfunctions:
            if result not in terms:
                terms[result] = None
            elif cancel:
                del terms[result]
            continue

        made_for, nested = pending[result]
        # Take the argument's terms over when they are the larger set. Only a
        # result made for `arg` itself: one passed through a one-argument gate
        # may belong to a node with other parents.
        if made_for is arg and len(nested) > len(terms):
            del pending[result]  # absorbed and unreachable: never written
            terms, nested = nested, terms
        for term in nested:
            if term not in terms:
                terms[term] = None
            elif cancel:
                del terms[term]

    # what is left may not need a gate at all: a single term stands for itself,
    # and an XOR whose terms all cancelled is the constant 0 (the empty sum)
    if len(terms) == 1:
        return next(iter(terms))
    if cancel and not terms:
        return CONST(0)

    result = node if in_place else gate(arg_limit = node.arg_limit)
    pending[result] = (node, terms)
    return result

class XOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, values, array):
        return reduce(
            lambda a, b: a ^ b,
            (values[arg] for arg in self.args)
        )
    def _eval_ANF(self, values, array):
        return reduce(
            lambda a, b: a ^ b,
            (values[arg] for arg in self.args)
        )

    def _generate_c(self, c_strings, array_name):
        return "(" + " ^ ".join(c_strings[arg] for arg in self.args) + ")"
    def _generate_VHDL(self, vhdl_strings, array_name):
        return "(" + " XOR ".join(vhdl_strings[arg] for arg in self.args) + ")"
    def _generate_python(self, python_strings, array_name):
        return "(" + " ^ ".join(python_strings[arg] for arg in self.args) + ")"


    def _merge_redundant(self, cache, subfunctions, in_place = False):
        return _merge_associative(self, XOR, cache, subfunctions, in_place, cancel = True)


    def _binarize(self,cache):
        return reduce(
            lambda a, b: XOR(a,b),
            (cache[arg] for arg in self.args)
        )

    @classmethod
    def tseytin_formula(cls, output, a, b):
        """The two-input clause kernel for XOR: four clauses pinning
        `output` to a XOR b.

        Named and shaped for this class alone -- the walker never calls it, and
        nothing outside `XOR.tseytin_unroll` depends on its arity.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the first argument wire.
        :type a: int
        :param b: Label of the second argument wire.
        :type b: int
        :return: Clauses forcing `output` to equal a XOR b.
        :rtype: list[tuple[int, ...]]
        """
        return [
            (-a,-b,-output),
            (a,b,-output),
            (a,-b,output),
            (-a,b,output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for an n-argument XOR.

        XOR is associative, so an n-input gate is a left-fold of the
        two-input exclusive-or: each intermediate wire combines the previous result
        with the next argument, and `tseytin_formula` supplies the clauses at
        every step. Degenerate arities are handled directly -- no arguments
        asserts a contradiction (a unit clause and its negation), one argument
        asserts equality with the output.

        :param label_map: Wires for every node encoded so far, including
            this one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the whole gate.
        :rtype: list[tuple[int, ...]]
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 0: # empty gate (e.g. XOR() ) => assert unsat using this gate
            return [(own_labels[0],),(-own_labels[0],)]

        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-own_labels[0],arg_labels[0]),(own_labels[0],-arg_labels[0])]

        else:
            # initial node uses 2 args
            clauses = self.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

            # afterward each node uses the previous + the next arg
            for i in range(1,len(own_labels)):
                clauses += self.tseytin_formula(
                    own_labels[i], own_labels[i-1], arg_labels[i+1]
                )
            return clauses




class AND(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, values, array):
        return reduce(
            lambda a, b: a & b,
            (values[arg] for arg in self.args)
        )
    def _eval_ANF(self, values, array):
        return reduce(
            lambda a, b: a & b,
            (values[arg] for arg in self.args)
        )


    def _generate_c(self, c_strings, array_name):
        return "(" + " & ".join(c_strings[arg] for arg in self.args) + ")"
    def _generate_VHDL(self, vhdl_strings, array_name):
        return "(" + " AND ".join(vhdl_strings[arg] for arg in self.args) + ")"
    def _generate_python(self, python_strings, array_name):
        return "(" + " & ".join(python_strings[arg] for arg in self.args) + ")"

    def _merge_redundant(self, cache, subfunctions, in_place = False):
        return _merge_associative(self, AND, cache, subfunctions, in_place, cancel = False)


    def _binarize(self, cache):
        return reduce(
            lambda a, b: AND(a,b),
            (cache[arg] for arg in self.args)
        )

    @classmethod
    def tseytin_formula(cls, output, a, b):
        """The two-input clause kernel for AND: four clauses pinning
        `output` to a AND b.

        Named and shaped for this class alone -- the walker never calls it, and
        nothing outside `AND.tseytin_unroll` depends on its arity.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the first argument wire.
        :type a: int
        :param b: Label of the second argument wire.
        :type b: int
        :return: Clauses forcing `output` to equal a AND b.
        :rtype: list[tuple[int, ...]]
        """
        return [
            (-a,-b,output),
            (a,-output),
            (b,-output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for an n-argument AND.

        AND is associative, so an n-input gate is a left-fold of the
        two-input conjunction: each intermediate wire combines the previous result
        with the next argument, and `tseytin_formula` supplies the clauses at
        every step. Degenerate arities are handled directly -- no arguments
        asserts a contradiction (a unit clause and its negation), one argument
        asserts equality with the output.

        :param label_map: Wires for every node encoded so far, including
            this one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the whole gate.
        :rtype: list[tuple[int, ...]]
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(own_labels[0],),(-own_labels[0],)]

        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-own_labels[0],arg_labels[0]),(own_labels[0],-arg_labels[0])]

        else:
            # initial node uses 2 args
            clauses = self.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

            # afterward each node uses the previous + the next arg
            for i in range(1,len(own_labels)):
                clauses += self.tseytin_formula(
                    own_labels[i], own_labels[i-1], arg_labels[i+1]
                )
            return clauses




class OR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _eval(self, values, array):
        return reduce(
            lambda a, b: a | b,
            (values[arg] for arg in self.args)
        )
    def _eval_ANF(self, values, array):
        return invert(reduce(
            lambda a, b: a & b,
            (invert(values[arg]) for arg in self.args)
        ))

    def _generate_c(self, c_strings, array_name):
        return "(" + " | ".join(c_strings[arg] for arg in self.args) + ")"
    def _generate_VHDL(self, vhdl_strings, array_name):
        return "(" + " OR ".join(vhdl_strings[arg] for arg in self.args) + ")"
    def _generate_python(self, python_strings, array_name):
        return "(" + " | ".join(python_strings[arg] for arg in self.args) + ")"


    def _merge_redundant(self, cache, subfunctions, in_place = False):
        return _merge_associative(self, OR, cache, subfunctions, in_place, cancel = False)


    def _binarize(self, cache):
        return reduce(
            lambda a, b: OR(a,b),
            (cache[arg] for arg in self.args)
        )

    @classmethod
    def tseytin_formula(cls, output, a, b):
        """The two-input clause kernel for OR: four clauses pinning
        `output` to a OR b.

        Named and shaped for this class alone -- the walker never calls it, and
        nothing outside `OR.tseytin_unroll` depends on its arity.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the first argument wire.
        :type a: int
        :param b: Label of the second argument wire.
        :type b: int
        :return: Clauses forcing `output` to equal a OR b.
        :rtype: list[tuple[int, ...]]
        """
        return [
            (a,b,-output),
            (-a,output),
            (-b,output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for an n-argument OR.

        OR is associative, so an n-input gate is a left-fold of the
        two-input disjunction: each intermediate wire combines the previous result
        with the next argument, and `tseytin_formula` supplies the clauses at
        every step. Degenerate arities are handled directly -- no arguments
        asserts a contradiction (a unit clause and its negation), one argument
        asserts equality with the output.

        :param label_map: Wires for every node encoded so far, including
            this one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the whole gate.
        :rtype: list[tuple[int, ...]]
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(own_labels[0],),(-own_labels[0],)]

        elif len(arg_labels) == 1: # 1-arg => assert arg and output var are equal
            return [(-own_labels[0],arg_labels[0]),(own_labels[0],-arg_labels[0])]

        else:
            # initial node uses 2 args
            clauses = self.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

            # afterward each node uses the previous + the next arg
            for i in range(1,len(own_labels)):
                clauses += self.tseytin_formula(
                    own_labels[i], own_labels[i-1], arg_labels[i+1]
                )
            return clauses




class XNOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _merge_redundant(self, cache, subfunctions, in_place = False):
        # with one argument this gate is a negation, so it cannot be dropped
        # the way a one-argument XOR, AND or OR can
        if len(self.args) == 1:
            return NOT(cache[self.args[0]])
        return super()._merge_redundant(cache, subfunctions, in_place)

    def _eval(self, values, array):
        return invert(reduce(
            lambda a, b: a ^ b,
            (values[arg] for arg in self.args)
        ))
    def _eval_ANF(self, values, array):
        return invert(reduce(
            lambda a, b: a ^ b,
            (values[arg] for arg in self.args)
        ))

    def _generate_c(self, c_strings, array_name):
        return "(!(" + " ^ ".join(c_strings[arg] for arg in self.args) + "))"
    def _generate_VHDL(self, vhdl_strings, array_name):
        # the native xnor is only right for exactly two arguments: VHDL reads a chained
        # xnor as (a xnor b) xnor c = a xor b xor c, and one argument
        # has to be negated -- otherwise negate the associative form
        if len(self.args) == 2:
            return "(" + " XNOR ".join(vhdl_strings[arg] for arg in self.args) + ")"
        return "(NOT(" + " XOR ".join(vhdl_strings[arg] for arg in self.args) + "))"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + " ^ ".join(python_strings[arg] for arg in self.args) + "))"

    def _binarize(self, cache):
        if len(self.args) == 1:
            return NOT(cache[self.args[0]])
        return XNOR(
            reduce(
                lambda a, b: XOR(a,b),
                (cache[arg] for arg in self.args[:-1])
            ),
            cache[self.args[-1]]
        )

    @classmethod
    def tseytin_formula(cls, output, a, b):
        """The two-input clause kernel for XNOR: four clauses pinning
        `output` to NOT(xor of its two inputs).

        Named and shaped for this class alone -- the walker never calls it, and
        nothing outside `XNOR.tseytin_unroll` depends on its arity.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the first argument wire.
        :type a: int
        :param b: Label of the second argument wire.
        :type b: int
        :return: Clauses forcing `output` to equal NOT(xor of its two inputs).
        :rtype: list[tuple[int, ...]]
        """
        return [
            (a,b,output),
            (-a,-b,output),
            (-a,b,-output),
            (a,-b,-output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for an n-argument XNOR.

        XNOR is **not** associative, so this is not a fold of two-input
        XNORs. An n-input XNOR is the negation of an n-input XOR, so the
        chain is built from `XOR.tseytin_formula` -- the associative core --
        and this class's own `tseytin_formula` is applied once, at the final
        wire, to negate the result. Degenerate arities are handled directly.

        :param label_map: Wires for every node encoded so far, including
            this one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the whole gate.
        :rtype: list[tuple[int, ...]]
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 0: # empty gate (e.g. XOR() ) => assert unsat using this gate
            return [(own_labels[0],),(-own_labels[0],)]

        elif len(arg_labels) == 1: # 1-arg => the negation of its one argument, as eval and codegen have it
            return [(-own_labels[0],-arg_labels[0]),(own_labels[0],arg_labels[0])]

        elif len(arg_labels) == 2: # 2-arg => just use formula
            return self.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

        else:
            # initial node uses 2 args
            clauses = XOR.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

            # afterward each node uses the previous + the next arg
            # using the associative operation and negating
            for i in range(1,len(own_labels)-1):
                clauses += XOR.tseytin_formula(
                    own_labels[i], own_labels[i-1], arg_labels[i+1]
                )

            # use negation for the final output
            idx = len(own_labels)-1
            clauses += self.tseytin_formula(
                own_labels[idx], own_labels[idx-1], arg_labels[idx+1]
            )

            return clauses




class NAND(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _merge_redundant(self, cache, subfunctions, in_place = False):
        # with one argument this gate is a negation, so it cannot be dropped
        # the way a one-argument XOR, AND or OR can
        if len(self.args) == 1:
            return NOT(cache[self.args[0]])
        return super()._merge_redundant(cache, subfunctions, in_place)

    def _eval(self, values, array):
        return invert(reduce(
            lambda a, b: a & b,
            (values[arg] for arg in self.args)
        ))
    def _eval_ANF(self, values, array):
        return invert(reduce(
            lambda a, b: a & b,
            (values[arg] for arg in self.args)
        ))

    def _generate_c(self, c_strings, array_name):
        return "(!(" + " & ".join(c_strings[arg] for arg in self.args) + "))"
    def _generate_VHDL(self, vhdl_strings, array_name):
        # the native nand is only right for exactly two arguments: VHDL forbids
        # chaining nand without parentheses, and one argument
        # has to be negated -- otherwise negate the associative form
        if len(self.args) == 2:
            return "(" + " NAND ".join(vhdl_strings[arg] for arg in self.args) + ")"
        return "(NOT(" + " AND ".join(vhdl_strings[arg] for arg in self.args) + "))"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + " & ".join(python_strings[arg] for arg in self.args) + "))"


    def _binarize(self, cache):
        if len(self.args) == 1:
            return NOT(cache[self.args[0]])
        return NAND(
            reduce(
                lambda a, b: AND(a,b),
                (cache[arg] for arg in self.args[:-1])
            ),
            cache[self.args[-1]]
        )

    @classmethod
    def tseytin_formula(cls, output, a, b):
        """The two-input clause kernel for NAND: four clauses pinning
        `output` to NOT(and of its two inputs).

        Named and shaped for this class alone -- the walker never calls it, and
        nothing outside `NAND.tseytin_unroll` depends on its arity.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the first argument wire.
        :type a: int
        :param b: Label of the second argument wire.
        :type b: int
        :return: Clauses forcing `output` to equal NOT(and of its two inputs).
        :rtype: list[tuple[int, ...]]
        """
        return [
            (-a,-b,-output),
            (a,output),
            (b,output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for an n-argument NAND.

        NAND is **not** associative, so this is not a fold of two-input
        NANDs. An n-input NAND is the negation of an n-input AND, so the
        chain is built from `AND.tseytin_formula` -- the associative core --
        and this class's own `tseytin_formula` is applied once, at the final
        wire, to negate the result. Degenerate arities are handled directly.

        :param label_map: Wires for every node encoded so far, including
            this one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the whole gate.
        :rtype: list[tuple[int, ...]]
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(own_labels[0],),(-own_labels[0],)]

        elif len(arg_labels) == 1: # 1-arg => the negation of its one argument, as eval and codegen have it
            return [(-own_labels[0],-arg_labels[0]),(own_labels[0],arg_labels[0])]

        elif len(arg_labels) == 2: # 2-arg => just use formula
            return self.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

        else:
            # initial node uses 2 args
            clauses = AND.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

            # afterward each node uses the previous + the next arg
            # using the associative operation and negating
            for i in range(1,len(own_labels)-1):
                clauses += AND.tseytin_formula(
                    own_labels[i], own_labels[i-1], arg_labels[i+1]
                )

            # use negation for the final output
            idx = len(own_labels)-1
            clauses += self.tseytin_formula(
                own_labels[idx], own_labels[idx-1], arg_labels[idx+1]
            )

            return clauses




class NOR(BooleanFunction):
    def __init__(self, *args, arg_limit = None):
        self.arg_limit = arg_limit
        self.args = args

    def _merge_redundant(self, cache, subfunctions, in_place = False):
        # with one argument this gate is a negation, so it cannot be dropped
        # the way a one-argument XOR, AND or OR can
        if len(self.args) == 1:
            return NOT(cache[self.args[0]])
        return super()._merge_redundant(cache, subfunctions, in_place)

    def _eval(self, values, array):
        return invert(reduce(
            lambda a, b: a | b,
            (values[arg] for arg in self.args)
        ))
    def _eval_ANF(self, values, array):
        return reduce(
            lambda a, b: a & b,
            (invert(values[arg]) for arg in self.args)
        )

    def _generate_c(self, c_strings, array_name):
        return "(!(" + " | ".join(c_strings[arg] for arg in self.args) + "))"
    def _generate_VHDL(self, vhdl_strings, array_name):
        # the native nor is only right for exactly two arguments: VHDL forbids
        # chaining nor without parentheses, and one argument
        # has to be negated -- otherwise negate the associative form
        if len(self.args) == 2:
            return "(" + " NOR ".join(vhdl_strings[arg] for arg in self.args) + ")"
        return "(NOT(" + " OR ".join(vhdl_strings[arg] for arg in self.args) + "))"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + " | ".join(python_strings[arg] for arg in self.args) + "))"



    def _binarize(self, cache):
        if len(self.args) == 1:
            return NOT(cache[self.args[0]])
        return NOR(
            reduce(
                lambda a, b: OR(a,b),
                (cache[arg] for arg in self.args[:-1])
            ),
            cache[self.args[-1]]
        )

    @classmethod
    def tseytin_formula(cls, output, a, b):
        """The two-input clause kernel for NOR: four clauses pinning
        `output` to NOT(or of its two inputs).

        Named and shaped for this class alone -- the walker never calls it, and
        nothing outside `NOR.tseytin_unroll` depends on its arity.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the first argument wire.
        :type a: int
        :param b: Label of the second argument wire.
        :type b: int
        :return: Clauses forcing `output` to equal NOT(or of its two inputs).
        :rtype: list[tuple[int, ...]]
        """
        return [
            (a,b,output),
            (-a,-output),
            (-b,-output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for an n-argument NOR.

        NOR is **not** associative, so this is not a fold of two-input
        NORs. An n-input NOR is the negation of an n-input OR, so the
        chain is built from `OR.tseytin_formula` -- the associative core --
        and this class's own `tseytin_formula` is applied once, at the final
        wire, to negate the result. Degenerate arities are handled directly.

        :param label_map: Wires for every node encoded so far, including
            this one and its arguments.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the whole gate.
        :rtype: list[tuple[int, ...]]
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 0: # empty gate (e.g. AND() ) => assert unsat using this gate
            return [(own_labels[0],),(-own_labels[0],)]

        elif len(arg_labels) == 1: # 1-arg => the negation of its one argument, as eval and codegen have it
            return [(-own_labels[0],-arg_labels[0]),(own_labels[0],arg_labels[0])]

        elif len(arg_labels) == 2: # 2-arg => just use formula
            return self.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

        else:
            # initial node uses 2 args
            clauses = OR.tseytin_formula(own_labels[0], arg_labels[0], arg_labels[1])

            # afterward each node uses the previous + the next arg
            # using the associative operation and negating
            for i in range(1,len(own_labels)-1):
                clauses += OR.tseytin_formula(
                    own_labels[i], own_labels[i-1], arg_labels[i+1]
                )

            # use negation for the final output
            idx = len(own_labels)-1
            clauses += self.tseytin_formula(
                own_labels[idx], own_labels[idx-1], arg_labels[idx+1]
            )

            return clauses




class NOT(BooleanFunction):
    def __init__(self, *args):
        if len(args) != 1:
            raise ValueError("NOT takes only 1 argument")
        self.arg_limit = 1
        self.args = args

    def _copy(self, child_copies):
        # NOT's constructor fixes arg_limit itself and takes no keyword for it
        return NOT(child_copies[self.args[0]])

    def _merge_redundant(self, cache, subfunctions, in_place=False):
        return NOT(cache[self.args[0]])


    def _eval(self, values, array):
        return invert(values[self.args[0]])
    def _eval_ANF(self, values, array):
        return invert(values[self.args[0]])

    def _generate_c(self, c_strings, array_name):
        return "(!(" + f"{c_strings[self.args[0]]}" + "))"
    def _generate_VHDL(self, vhdl_strings, array_name):
        return "(NOT(" + f"{vhdl_strings[self.args[0]]}" + "))"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + f"{python_strings[self.args[0]]}" + "))"


    def _binarize(self, cache):
        return NOT(cache[self.args[0]])

    @classmethod
    def tseytin_formula(cls,output,a):
        """The one-input clause kernel for NOT: two clauses pinning `output`
        opposite to `a`.

        Unary, unlike the two-input kernels the chaining gates carry under this
        name. NOT has `arg_limit = 1` and never chains, so there is no fold for
        a two-input kernel to serve.

        :param output: Label of the wire carrying the gate's result.
        :type output: int
        :param a: Label of the single argument wire.
        :type a: int
        :return: Clauses forcing `output` to equal NOT a.
        :rtype: list[tuple[int, ...]]
        """
        return [
            (-a,-output),
            (a,output)
        ]

    def _tseytin_clauses(self,label_map):
        """Clauses for a NOT.

        There is nothing to unroll: `arg_limit = 1` means the only real case is
        a single argument, whose clauses are emitted directly. The empty case is
        kept as a contradiction for symmetry with the chaining gates, and any
        larger arity is rejected rather than silently returning nothing.

        :param label_map: Wires for every node encoded so far, including
            this one and its argument.
        :type label_map: dict[BooleanFunction, list[int]]
        :return: Clauses encoding the negation.
        :rtype: list[tuple[int, ...]]
        :raises ValueError: If given more than one argument wire.
        """
        own_labels = label_map[self]
        arg_labels = [label_map[arg][-1] for arg in self.args]
        if len(arg_labels) == 1: # 1-arg => assert arg and output var are opposite
            return [(-own_labels[0],-arg_labels[0]),(own_labels[0],arg_labels[0])]

        if len(arg_labels) == 0: # should never happen, assert unsat
            return [(own_labels[0],),(-own_labels[0],)]

        # NOT has arg_limit = 1, so add_arguments rejects a second argument
        # before this can be reached. Raising rather than falling through keeps
        # the declared return type honest: every path returns clauses or fails.
        raise ValueError(
            f"{type(self).__name__} encodes a single argument, "
            f"but {len(arg_labels)} argument wires were given"
        )



