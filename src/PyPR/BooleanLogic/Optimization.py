"""Logic optimization for several Boolean functions at once.

`optimize` rewrites a list of functions -- the bits of a feedback function, say
-- as one DAG that computes the same outputs with fewer two-input operations,
sharing whatever logic the outputs have in common. `generate_python(optimize=True)`
writes the result, on a single function or a whole feedback function, and
`compile` passes that flag, since compiled code evaluates exactly that DAG.

Finding a smallest program is NP-hard even for XOR alone (the shortest linear
program problem), so every step here is a heuristic that never increases the
operation count. See docs/architecture/Logic Optimization.md for the network
representation, each step's rule and its source, and the cost model.
"""
import copy
import heapq
from collections.abc import Callable, Sequence
from numbers import Integral
from typing import Any

from PyPR.BooleanLogic.BooleanFunction import BooleanFunction
from PyPR.BooleanLogic.FunctionInputs import CONST, VAR
from PyPR.BooleanLogic.Gates import AND, NAND, NOR, NOT, OR, XNOR, XOR

# A literal is a node of the network or its complement: `2 * node + negated`.
# Node 0 is the constant 0, so its two literals are the constants.
_FALSE = 0
_TRUE = 1

# Node kinds. A gate's terms are a sorted tuple of literals; XOR terms are
# always positive, since a complemented input only flips the output.
_CONST = "const"
_VAR = "var"
_ATOM = "atom"
_XOR = "xor"
_AND = "and"
_OR = "or"

# How many term pairs pair extraction will count, per gate kind. Counting is
# quadratic in a gate's size, so this is what bounds the pass on huge gates;
# gates that do not fit are left as they are.
PAIR_BUDGET = 400_000


class _Network:
    """A hash-consed network of XOR, AND and OR gates over variables and atoms.

    Building a gate through `make` simplifies it first and then returns the
    existing node if an identical gate was built before, so equal logic is one
    node however many times it is asked for.
    """

    def __init__(self) -> None:
        self.kind: list[str] = [_CONST]
        self.terms: list[tuple[int, ...]] = [()]
        # a variable's index, or the original node an atom stands for
        self.payload: list[Any] = [None]
        self._table: dict[tuple, int] = {}

    def _add(self, kind: str, terms: tuple[int, ...], payload: Any = None) -> int:
        self.kind.append(kind)
        self.terms.append(terms)
        self.payload.append(payload)
        return len(self.kind) - 1

    def _unique(self, kind: str, terms: tuple[int, ...]) -> int:
        key = (kind, terms)
        if key not in self._table:
            self._table[key] = self._add(kind, terms)
        return self._table[key]

    def var(self, index: int) -> int:
        key = (_VAR, index)
        if key not in self._table:
            self._table[key] = self._add(_VAR, (), index)
        return self._table[key] << 1

    def atom(self, node: BooleanFunction, terms: Sequence[int]) -> int:
        """A node this module does not understand: kept as it is, over optimized arguments."""
        return self._add(_ATOM, tuple(terms), node) << 1

    def make(self, kind: str, literals: Sequence[int]) -> int:
        if kind == _XOR:
            return self.make_xor(literals)
        if kind == _AND:
            return self.make_and(literals)
        return self.make_or(literals)

    def make_xor(self, literals: Sequence[int]) -> int:
        # x ^ ~y = ~(x ^ y), so complements move to the output; x ^ x = 0, so a
        # node occurring twice cancels; x ^ 0 = x
        negated = 0
        nodes: set[int] = set()
        for literal in literals:
            negated ^= literal & 1
            node = literal >> 1
            if node == 0:
                continue
            if node in nodes:
                nodes.remove(node)
            else:
                nodes.add(node)
        if not nodes:
            return negated
        if len(nodes) == 1:
            return (next(iter(nodes)) << 1) | negated
        terms = tuple(sorted(node << 1 for node in nodes))
        return (self._unique(_XOR, terms) << 1) | negated

    def make_and(self, literals: Sequence[int]) -> int:
        # x & 0 = 0, x & 1 = x, x & x = x, x & ~x = 0
        chosen: set[int] = set()
        for literal in literals:
            if literal == _FALSE:
                return _FALSE
            if literal == _TRUE:
                continue
            if literal ^ 1 in chosen:
                return _FALSE
            chosen.add(literal)
        if not chosen:
            return _TRUE
        if len(chosen) == 1:
            return next(iter(chosen))
        # De Morgan: ~a & ~b = ~(a | b), one complement instead of one per input
        if all(literal & 1 for literal in chosen):
            return self.make_or([literal ^ 1 for literal in chosen]) ^ 1
        return self._unique(_AND, tuple(sorted(chosen))) << 1

    def make_or(self, literals: Sequence[int]) -> int:
        # x | 1 = 1, x | 0 = x, x | x = x, x | ~x = 1
        chosen: set[int] = set()
        for literal in literals:
            if literal == _TRUE:
                return _TRUE
            if literal == _FALSE:
                continue
            if literal ^ 1 in chosen:
                return _TRUE
            chosen.add(literal)
        if not chosen:
            return _FALSE
        if len(chosen) == 1:
            return next(iter(chosen))
        # De Morgan: ~a | ~b = ~(a & b)
        if all(literal & 1 for literal in chosen):
            return self.make_and([literal ^ 1 for literal in chosen]) ^ 1
        return self._unique(_OR, tuple(sorted(chosen))) << 1


def _reachable(net: _Network, roots: Sequence[int]) -> list[int]:
    """The nodes the root literals depend on, every node after the nodes it uses."""
    order: list[int] = []
    seen: set[int] = set()
    stack: list[int] = [literal >> 1 for literal in reversed(roots)]
    while stack:
        node = stack.pop()
        if node < 0:
            order.append(~node)
            continue
        if node in seen:
            continue
        seen.add(node)
        stack.append(~node)
        for literal in reversed(net.terms[node]):
            if literal >> 1 not in seen:
                stack.append(literal >> 1)
    return order


def _fanout(net: _Network, roots: Sequence[int], order: Sequence[int]) -> dict[int, int]:
    """How many times each node is used: by a gate's term, or as an output."""
    uses = dict.fromkeys(order, 0)
    for node in order:
        for literal in net.terms[node]:
            uses[literal >> 1] += 1
    for literal in roots:
        uses[literal >> 1] += 1
    return uses


def _used_literals(net: _Network, roots: Sequence[int], order: Sequence[int]) -> set[int]:
    """Every literal some gate takes as a term, or some output is."""
    used = set(roots)
    for node in order:
        used.update(net.terms[node])
    return used


def _free_complements(net: _Network, order: Sequence[int], used: set[int]) -> dict[int, int]:
    """XOR gates whose complement costs nothing, each with the term that makes it so.

    ~(a ^ b) = ~a ^ b. So when an XOR is only ever used complemented, and ~a is
    computed anyway for some term a, the gate can be built over ~a directly and
    needs no complement of its own.
    """
    free: dict[int, int] = {}
    for node in order:
        if net.kind[node] != _XOR:
            continue
        if (node << 1) in used or (node << 1 | 1) not in used:
            continue
        for term in net.terms[node]:
            if term ^ 1 in used:
                free[node] = term
                break
    return free


def _cost(net: _Network, roots: Sequence[int]) -> int:
    """Two-input operations needed to evaluate the roots, each node computed once.

    A gate of n terms is n - 1 operations. A complement is one operation per
    node that is used complemented anywhere, since the complement is shared,
    except where `_free_complements` shows it costs nothing.
    """
    order = _reachable(net, roots)
    used = _used_literals(net, roots, order)
    operations = 0
    for node in order:
        if net.kind[node] in (_XOR, _AND, _OR):
            operations += len(net.terms[node]) - 1
    complemented = {literal >> 1 for literal in used if literal & 1 and literal >> 1 != 0}
    return operations + len(complemented) - len(_free_complements(net, order, used))


def _build(functions: Sequence[BooleanFunction]) -> tuple[_Network, list[int]]:
    """The network computing `functions`, with a root literal for each."""
    net = _Network()
    literal_of: dict[BooleanFunction, int] = {}
    makers: dict[type, tuple[str, int]] = {
        XOR: (_XOR, 0), AND: (_AND, 0), OR: (_OR, 0),
        XNOR: (_XOR, 1), NAND: (_AND, 1), NOR: (_OR, 1),
    }
    roots = []
    for function in functions:
        for node in function.postorder(stop=literal_of.__contains__):
            # exact types only: a subclass may mean something else, and is an atom
            node_type = type(node)
            terms = [literal_of[arg] for arg in node.args]
            if isinstance(node, VAR) and node_type is VAR:
                literal_of[node] = net.var(node.index)
            elif (isinstance(node, CONST) and node_type is CONST
                    and isinstance(node.value, Integral) and node.value in (0, 1)):
                literal_of[node] = int(node.value)
            elif node_type is NOT:
                literal_of[node] = terms[0] ^ 1
            elif node_type in makers and terms:
                kind, negated = makers[node_type]
                literal_of[node] = net.make(kind, terms) ^ negated
            else:
                literal_of[node] = net.atom(node, terms)
        roots.append(literal_of[function])
    return net, roots


def _rebuild(net: _Network, roots: Sequence[int]) -> tuple[_Network, list[int]]:
    """Rebuild a network through `make`, merging gates into the gates that use them.

    A gate used exactly once, by a gate of the same kind, is absorbed into it:
    (a ^ b) ^ c becomes a ^ b ^ c. A gate with other users is kept whole, since
    they share it. Rebuilding also re-simplifies and re-hashes every gate, which
    the other passes rely on after editing gates in place.
    """
    order = _reachable(net, roots)
    uses = _fanout(net, roots, order)
    new = _Network()
    moved: dict[int, int] = {0: _FALSE}
    for node in order:
        kind = net.kind[node]
        if kind == _CONST:
            continue
        if kind == _VAR:
            moved[node] = new.var(net.payload[node])
            continue
        if kind == _ATOM:
            terms = [moved[literal >> 1] ^ (literal & 1) for literal in net.terms[node]]
            moved[node] = new.atom(net.payload[node], terms)
            continue

        literals: list[int] = []
        for literal in net.terms[node]:
            child = literal >> 1
            new_literal = moved[child] ^ (literal & 1)
            absorb = (
                uses[child] == 1
                and net.kind[child] == kind
                and new.kind[new_literal >> 1] == kind
                # an XOR's complement moves to the output; AND and OR need the
                # child uncomplemented to be the same operation
                and (kind == _XOR or new_literal & 1 == 0)
            )
            if absorb:
                literals.extend(new.terms[new_literal >> 1])
                if new_literal & 1:
                    literals.append(_TRUE)
            else:
                literals.append(new_literal)
        moved[node] = new.make(kind, literals)
    return new, [moved[literal >> 1] ^ (literal & 1) for literal in roots]


def _extract_pairs(net: _Network, roots: Sequence[int], kind: str) -> None:
    """Share pairs of terms between gates of one kind, greedily, most frequent first.

    If terms a and b occur together in several gates, a new gate p = a . b is
    made once and each of those gates uses p in their place: k gates sharing
    the pair save k - 1 operations. This is Paar's greedy algorithm for XOR
    networks, applied to AND and OR as well, where a shared pair is a common
    two-literal cube. Edits `net` in place; the caller rebuilds it.
    """
    order = _reachable(net, roots)
    gates = [node for node in order if net.kind[node] == kind and len(net.terms[node]) >= 2]
    terms = {gate: set(net.terms[gate]) for gate in gates}

    occurrences: dict[int, set[int]] = {}
    for gate in gates:
        for term in terms[gate]:
            if term not in occurrences:
                occurrences[term] = set()
            occurrences[term].add(gate)

    # only a term in two or more gates can be part of a shared pair
    budget = PAIR_BUDGET
    counts: dict[tuple[int, int], int] = {}
    counted: set[int] = set()
    for gate in sorted(gates, key=lambda gate: (len(terms[gate]), gate)):
        shared = sorted(term for term in terms[gate] if len(occurrences[term]) >= 2)
        pairs = len(shared) * (len(shared) - 1) // 2
        if pairs > budget:
            continue
        budget -= pairs
        counted.add(gate)
        for i, first in enumerate(shared):
            for second in shared[i + 1:]:
                pair = (first, second)
                if pair in counts:
                    counts[pair] += 1
                else:
                    counts[pair] = 1
    if len(counted) < len(gates):
        for members in occurrences.values():
            members &= counted

    # a gate that is exactly a pair already computes it
    existing = {tuple(sorted(terms[gate])): gate for gate in gates if len(terms[gate]) == 2}

    heap = [(-count, pair) for pair, count in counts.items() if count >= 2]
    heapq.heapify(heap)

    def change(first: int, second: int, amount: int) -> None:
        pair = (first, second) if first < second else (second, first)
        if pair in counts:
            count = counts[pair] + amount
        elif amount > 0:
            count = amount
        else:
            # never counted (one of its terms was in a single gate): nothing to lower
            return
        counts[pair] = count
        if count >= 2:
            heapq.heappush(heap, (-count, pair))

    while heap and budget > 0:
        negative_count, pair = heapq.heappop(heap)
        # entries are not updated in place, so skip any that are out of date
        if pair not in counts or counts[pair] != -negative_count:
            continue
        first, second = pair
        del counts[pair]
        users = sorted(occurrences[first] & occurrences[second])
        if len(users) < 2:
            continue
        if pair in existing:
            shared_gate = existing[pair]
        else:
            shared_gate = net._add(kind, pair)
            existing[pair] = shared_gate
        new_term = shared_gate << 1

        for gate in users:
            if gate == shared_gate:
                continue
            gate_terms = terms[gate]
            gate_terms.remove(first)
            gate_terms.remove(second)
            occurrences[first].discard(gate)
            occurrences[second].discard(gate)
            budget -= len(gate_terms)
            for term in gate_terms:
                change(first, term, -1)
                change(second, term, -1)

            if new_term in gate_terms:
                # the gate already used the shared pair as a term: a second
                # copy cancels it in an XOR and changes nothing in AND or OR
                if kind == _XOR:
                    gate_terms.remove(new_term)
                    occurrences[new_term].discard(gate)
                    for term in gate_terms:
                        change(new_term, term, -1)
            else:
                for term in gate_terms:
                    change(new_term, term, 1)
                gate_terms.add(new_term)
                if new_term not in occurrences:
                    occurrences[new_term] = set()
                occurrences[new_term].add(gate)

            if len(gate_terms) == 2:
                key = tuple(sorted(gate_terms))
                if key not in existing:
                    existing[key] = gate

    for gate in gates:
        net.terms[gate] = tuple(sorted(terms[gate]))


def _factor_terms(
    net: _Network,
    outer: str,
    inner: str,
    cubes: list[frozenset[int]],
    others: list[int],
) -> int:
    """Build `outer` over `others` and the `inner` gates in `cubes`, factoring out common literals.

    While some literal x occurs in two or more cubes, those cubes are replaced
    by one: with outer = XOR and inner = AND, x&a ^ x&b becomes x & (a ^ b),
    saving an AND for each cube after the first. The cofactors a, b are
    factored in turn. The same identity holds for (OR, AND) and (AND, OR).
    """
    cubes = list(cubes)
    others = list(others)
    while True:
        frequency: dict[int, int] = {}
        for cube in cubes:
            for literal in cube:
                if literal in frequency:
                    frequency[literal] += 1
                else:
                    frequency[literal] = 1
        best = None
        for literal in sorted(frequency):
            if frequency[literal] >= 2 and (best is None or frequency[literal] > frequency[best]):
                best = literal
        if best is None:
            break

        inner_cubes: list[frozenset[int]] = []
        inner_others: list[int] = []
        remaining: list[frozenset[int]] = []
        for cube in cubes:
            if best not in cube:
                remaining.append(cube)
                continue
            cofactor = cube - {best}
            if len(cofactor) == 1:
                inner_others.append(next(iter(cofactor)))
            else:
                inner_cubes.append(cofactor)
        cubes = remaining
        cofactors = _factor_terms(net, outer, inner, inner_cubes, inner_others)
        others.append(net.make(inner, [best, cofactors]))

    return net.make(outer, others + [net.make(inner, sorted(cube)) for cube in cubes])


def _factor(net: _Network, roots: Sequence[int]) -> tuple[_Network, list[int]]:
    """Factor common literals out of sums of products (and products of sums).

    Only an inner gate used once, by the gate being factored, takes part: one
    with other users is shared logic, and factoring it would compute it twice.
    """
    order = _reachable(net, roots)
    uses = _fanout(net, roots, order)
    inner_of = {_XOR: _AND, _OR: _AND, _AND: _OR}
    new = _Network()
    moved: dict[int, int] = {0: _FALSE}
    for node in order:
        kind = net.kind[node]
        if kind == _CONST:
            continue
        if kind == _VAR:
            moved[node] = new.var(net.payload[node])
            continue
        terms = [moved[literal >> 1] ^ (literal & 1) for literal in net.terms[node]]
        if kind == _ATOM:
            moved[node] = new.atom(net.payload[node], terms)
            continue

        inner = inner_of[kind]
        cubes: list[frozenset[int]] = []
        others: list[int] = []
        for literal, new_literal in zip(net.terms[node], terms, strict=True):
            child = literal >> 1
            is_cube = (
                literal & 1 == 0
                and net.kind[child] == inner
                and uses[child] == 1
                and new_literal & 1 == 0
                and new.kind[new_literal >> 1] == inner
            )
            if is_cube:
                cubes.append(frozenset(new.terms[new_literal >> 1]))
            else:
                others.append(new_literal)
        moved[node] = _factor_terms(new, kind, inner, cubes, others)
    return new, [moved[literal >> 1] ^ (literal & 1) for literal in roots]


def _to_functions(net: _Network, roots: Sequence[int]) -> list[BooleanFunction]:
    """The network as BooleanFunctions: one object per node, shared by everything that uses it."""
    gate_types: dict[str, Callable[..., BooleanFunction]] = {_XOR: XOR, _AND: AND, _OR: OR}
    order = _reachable(net, roots)
    free = _free_complements(net, order, _used_literals(net, roots, order))
    positive: dict[int, BooleanFunction] = {}
    complements: dict[int, BooleanFunction] = {}

    def function_of(literal: int) -> BooleanFunction:
        if literal >> 1 == 0:
            return CONST(literal)
        if literal & 1 == 0:
            return positive[literal >> 1]
        if literal >> 1 not in complements:
            complements[literal >> 1] = NOT(positive[literal >> 1])
        return complements[literal >> 1]

    for node in order:
        kind = net.kind[node]
        if kind == _CONST:
            continue
        if kind == _VAR:
            positive[node] = VAR(net.payload[node])
        elif kind == _ATOM:
            original: BooleanFunction = net.payload[node]
            if original.is_leaf():
                positive[node] = original.copy()
            else:
                args = [function_of(literal) for literal in net.terms[node]]
                positive[node] = original._copy(dict(zip(original.args, args, strict=True)))
        elif node in free:
            # only its complement is used, and it is built over a complemented term
            args = [function_of(literal ^ 1 if literal == free[node] else literal)
                    for literal in net.terms[node]]
            complements[node] = XOR(*args)
        else:
            positive[node] = gate_types[kind](*(function_of(literal) for literal in net.terms[node]))
    return [function_of(literal) for literal in roots]


def optimize(functions: Sequence[BooleanFunction]) -> list[BooleanFunction]:
    """Rewrite several functions as one DAG computing the same outputs in fewer operations.

    The functions are treated as Boolean: every variable takes the value 0 or
    1, and a `CONST` of 0 or 1 is that bit. Under that reading the result
    computes exactly the same outputs. Identical logic is built once, wherever
    it occurs; constants, cancelling pairs (x ^ x), repeats (x & x) and
    contradictions (x & ~x) are simplified away; pairs of terms that several
    gates share are computed once; and common literals are factored out of
    sums of products. Nodes of any other class, and constants that are not
    bits, are kept as they are, over optimized arguments.

    The result is a heuristic: it never has more operations than the input
    (as counted by `operation_count`), but it is not guaranteed to be the
    fewest possible.

    :param functions: The functions to optimize together, such as a feedback
        function's `fn_list`.
    :type functions: Sequence[BooleanFunction]
    :return: One function per input, in the same order. They share nodes: a
        node used by several outputs is one object. The inputs are unchanged.
    :rtype: list[BooleanFunction]
    """
    net, roots = _rebuild(*_build(functions))
    best = _cost(net, roots)

    # Each pass is kept only if it lowers the count, so the result is never
    # worse than plain sharing. Extraction runs again after factoring, since
    # the cofactor sums factoring creates can share pairs of their own.
    def extract(net: _Network, roots: list[int]) -> tuple[_Network, list[int]]:
        edited = _Network()
        edited.kind = list(net.kind)
        edited.terms = list(net.terms)
        edited.payload = list(net.payload)
        for kind in (_XOR, _AND, _OR):
            _extract_pairs(edited, roots, kind)
        return _rebuild(edited, roots)

    for step in (extract, _factor, extract):
        candidate, candidate_roots = step(net, roots)
        cost = _cost(candidate, candidate_roots)
        if cost < best:
            net, roots, best = candidate, candidate_roots, cost

    result = _to_functions(net, roots)

    # The network's cost model and `operation_count` agree except in corners
    # (a complement the input happened to place better), so check the promise
    # itself and fall back to the input's own structure if it would be broken.
    if operation_count(result) > operation_count(functions):
        return copy.deepcopy(list(functions))
    return result


def operation_count(functions: Sequence[BooleanFunction]) -> int:
    """Count the two-input operations a list of functions takes to evaluate together.

    Each node object is counted once however many of the functions use it, so
    this measures the DAG as it is shared, not as it would be written out. A
    gate with n arguments is n - 1 operations, a NOT is one, and a negated gate
    (XNOR, NAND, NOR) is one more than its plain form.

    :param functions: The functions to count over.
    :type functions: Sequence[BooleanFunction]
    :return: The number of operations.
    :rtype: int
    """
    seen: set[BooleanFunction] = set()
    operations = 0
    for function in functions:
        for node in function.postorder(stop=seen.__contains__):
            seen.add(node)
            node_type = type(node)
            if node_type is NOT:
                operations += 1
            elif node_type in (XOR, AND, OR):
                operations += max(len(node.args) - 1, 0)
            elif node_type in (XNOR, NAND, NOR):
                operations += max(len(node.args) - 1, 0) + 1
    return operations
