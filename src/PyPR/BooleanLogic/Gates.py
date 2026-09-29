from functools import reduce

from PyPR.BooleanLogic import BooleanFunction


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
    def _generate_tex(self, cache, array_name):
        return " \\oplus \\,".join(cache[arg] for arg in self.args)


    def _merge_redundant(self,
        cache,
        subfunctions,
        in_place = False
    )->"BooleanFunction":
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
            self.args = tuple(new_args)
            return self
        else:
            return XOR(*new_args, arg_limit = self.arg_limit)


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
    def _generate_tex(self, cache, array_name):
        return "".join(cache for arg in self.args)

    def _merge_redundant(self, cache, subfunctions, in_place = False):
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
            self.args = tuple(new_args)
            return self
        else:
            return AND(*new_args, arg_limit = self.arg_limit)


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
    def _generate_tex(self, cache, array_name):
        return " \\vee ".join(cache[arg] for arg in self.args)


    def _merge_redundant(self, cache, subfunctions, in_place = False):
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
            self.args = tuple(new_args)
            return self
        else:
            return OR(*new_args, arg_limit = self.arg_limit)


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
        return "(" + " XNOR ".join(vhdl_strings[arg] for arg in self.args) + ")"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + " ^ ".join(python_strings[arg] for arg in self.args) + "))"

    def _binarize(self, cache):
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
        return "(" + " NAND ".join(vhdl_strings[arg] for arg in self.args) + ")"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + " & ".join(python_strings[arg] for arg in self.args) + "))"


    def _binarize(self, cache):
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
        return "(" + " NOR ".join(vhdl_strings[arg] for arg in self.args) + ")"
    def _generate_python(self, python_strings, array_name):
        return "(1-(" + " | ".join(python_strings[arg] for arg in self.args) + "))"



    def _binarize(self, cache):
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



