from math import gcd, lcm

from PyPR.Reporting import get_logger

from PyPR.Tools.RootCounting.Combinatorics import binsum, powerset
from PyPR.Tools.RootCounting.EvaluationFunctions import pessimistic_expected_value
from PyPR.Tools.RootCounting.JordanSet import JordanSet
from PyPR.Tools.RootCounting.OverlappingRectangle import rectangle_solve

log = get_logger(__name__)


def _coprime_groups(basis):
    """Group the field degrees of a basis into classes connected by a common factor.

    Two degrees are linked when their gcd exceeds 1; the groups are the connected
    components of that relation. Fields in different groups have pairwise coprime
    degrees, so products across groups are injective (Proposition 3); within a
    group they need not be.

    :param basis: Field degrees, as filed in the basis table.
    :type basis: tuple[int, ...]
    :return: Lists of positions into `basis`, one list per group.
    :rtype: list[list[int]]
    """
    parent = list(range(len(basis)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(basis)):
        for j in range(i + 1, len(basis)):
            if gcd(basis[i], basis[j]) > 1:
                parent[find(i)] = find(j)
    groups = {}
    for i in range(len(basis)):
        root = find(i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)
    return list(groups.values())


def _group_room(basis, group, pessimistic):
    """How many roots the products over one group of fields can take at most.

    The multiplicative groups GF(2^e)^* of the group's fields generate, inside
    the cyclic group GF(2^L)^* (L the lcm of the degrees), the subgroup of order
    lcm(2^e - 1), and every product lies in it. This holds whatever primitive
    polynomials the fields come from. A bound through exponent weights (a ball)
    would need the fields' generators to be compatible across the embeddings,
    which a CMPR's independently chosen blocks are not.

    The pessimistic version applies the same statistical model as
    `pessimistic_expected_value`: the root 1 is left out, and the count is
    scaled by (1 - 2^-L)^2.

    :param basis: Field degrees of the term.
    :type basis: tuple[int, ...]
    :param group: Positions into `basis` of one group (`_coprime_groups`).
    :type group: list[int]
    :param pessimistic: Whether to apply the statistical model.
    :type pessimistic: bool
    :return: The bound for the group.
    :rtype: int | float
    """
    order = lcm(*(2**basis[i] - 1 for i in group))
    if not pessimistic:
        return order
    degree = lcm(*(basis[i] for i in group))
    return (order - 1) * (2**degree - 1)**2 / 2**(2 * degree)


def _capped_basis_count(basis, rectangles, pessimistic):
    """Bound the count of one filed basis whose degrees are not pairwise coprime.

    `rectangle_solve` measures the union of the terms as boxes, which treats
    every tuple of roots as a different product. When degrees share a factor,
    different tuples can give the same product. Two further bounds hold
    regardless:
    * the sum over terms of each term's Jordan length times, per group of
      fields sharing a factor, the smaller of the group's box and its room
      (`_group_room`), groups multiplying since their degrees are coprime;
    * the largest Jordan length times the room of every group.

    The smallest of the three is returned.

    :param basis: Field degrees of the basis.
    :type basis: tuple[int, ...]
    :param rectangles: The filed terms: Jordan length first, then one count per field.
    :type rectangles: list[tuple]
    :param pessimistic: Whether the counts follow the statistical model (`lower`).
    :type pessimistic: bool
    :return: The capped count.
    :rtype: int | float
    """
    groups = _coprime_groups(basis)
    room = [_group_room(basis, group, pessimistic) for group in groups]
    per_term = 0
    for counts in rectangles:
        size = counts[0]
        for group, group_room in zip(groups, room):
            box = 1
            for i in group:
                box *= counts[i + 1]
            size *= min(box, group_room)
        per_term += size
    whole = max(counts[0] for counts in rectangles)
    for group_room in room:
        whole *= group_room
    return min(rectangle_solve(rectangles), per_term, whole)


def _warn_not_coprime(bases):
    log.warning(
        "Root expression has products over fields whose degrees are not "
        "pairwise coprime (%s); they are counted with min(box, the subgroup the "
        "fields generate), which is safe but may be loose.",
        ", ".join(str(b) for b in sorted(set(bases))),
    )


def is_embedded_subset(a,b):
    a_keyset = set(a.roots.keys())
    b_keyset = set(b.roots.keys())

    # there should be a strict subset relationship on keys:
    if not a_keyset.issubset(b_keyset):
        return False

    # all shared keys should have A <= B .
    for k in b_keyset.intersection(a_keyset):
        if a.roots[k] > b.roots[k]:
            return False

    # all non shared keys should be maximum (for each pair (b,c), b should equal c).
    for k in (b_keyset - a_keyset):
        if b.roots[k] != k:
            return False

    # if all conditions pass, return true
    return True

def is_root_subset(a,b):
    # all counts in A must be <= B to be a subset.
    for k in a.roots.keys():
        if a.roots[k] > b.roots[k]:
            return False
    return True


class RootExpression:
    def __init__(self,root_table):
        self.root_table = root_table

    def extend(self, extension_jordan_set):
        max_mult = 0
        for field_tuple in self.root_table:
            for term in self.root_table[field_tuple]:
                if is_embedded_subset(extension_jordan_set, term):
                    max_mult = max(max_mult,*term.mults)

        field_tuple = tuple(sorted(extension_jordan_set.roots.keys()))
        return self + RootExpression({
            field_tuple: {JordanSet(extension_jordan_set.roots, {max_mult + 1})}
        })

    def __str__(self):
        outstr = "{\n"
        for field_tuple in self.root_table:
            outstr += f"{field_tuple}: " + " + ".join(str(term) for term in self.root_table[field_tuple])
            outstr += "\n"
        outstr += "}"
        return outstr

    def __copy__(self):
        new_table = {}
        for field_tuple in self.root_table:
            new_table[field_tuple] = set()
            for JSet in self.root_table[field_tuple]:
                new_table[field_tuple].add(JSet.__copy__())
        return RootExpression(new_table)



    def __xor__(self, other): return self.__add__(other)
    def __add__(self, other):
        # clean out redundant subsets and merge.

        # add selfs terms
        output = self.__copy__()
        for field_tuple in other.root_table:
            # if there is no conflict, just write directly:
            if field_tuple not in self.root_table:
                output.root_table[field_tuple] = set()
                for JSet in other.root_table[field_tuple]:
                    output.root_table[field_tuple].add(JSet)
                continue

            # otherwise, split mults:
            # problem there could be smaller things, which were kept for mult reasons
            for JSet in other.root_table[field_tuple]:
                new_term = JSet.__copy__()

                # add in the new term:
                to_remove = []
                merge_location = None
                for already_added in output.root_table[field_tuple]:
                    smaller = is_root_subset(new_term, already_added)
                    bigger = is_root_subset(already_added,new_term)

                    if smaller and bigger:
                        merge_location = already_added
                    elif smaller:
                        new_term.mults -= already_added.mults
                        # can break early in this case
                        if not new_term.mults:
                            break
                    elif bigger:
                        already_added.mults -= new_term.mults
                        if not already_added.mults:
                            to_remove.append(already_added)

                # update the table entry accordingly:
                if new_term.mults and merge_location is not None:
                    merge_location.mults |= new_term.mults
                elif new_term.mults and merge_location is None:
                    output.root_table[field_tuple].add(new_term)
                for term in to_remove:
                    output.root_table[field_tuple].remove(term)

        #print("\n")
        return output


    def __and__(self, other): return self.__mul__(other)
    def __mul__(self, other):
        # clean out redundant subsets and merge.

        # add selfs terms
        output = RootExpression({})
        for field_tuple1 in self.root_table:
            for field_tuple2 in other.root_table:
                target_field_tuple = tuple(sorted(set(field_tuple1) | set(field_tuple2)))

                # add if needed
                if target_field_tuple not in output.root_table:
                    output.root_table[target_field_tuple] = set()

                for JSet1 in self.root_table[field_tuple1]:
                    for JSet2 in other.root_table[field_tuple2]:
                        prod = JSet1 * JSet2

                        # add in the new term:
                        to_remove = []
                        merge_location = None
                        for already_added in output.root_table[target_field_tuple]:
                            smaller = is_root_subset(prod, already_added)
                            bigger = is_root_subset(already_added,prod)

                            if smaller and bigger:
                                merge_location = already_added
                            elif smaller:
                                prod.mults -= already_added.mults
                                # can break early in this case
                                if not prod.mults:
                                    break
                            elif bigger:
                                already_added.mults -= prod.mults
                                if not already_added.mults:
                                    to_remove.append(already_added)

                        # update the table entry accordingly:
                        if prod.mults and merge_location is not None:
                            merge_location.mults |= prod.mults
                        elif prod.mults and merge_location is None:
                            output.root_table[target_field_tuple].add(prod)
                        for term in to_remove:
                            output.root_table[target_field_tuple].remove(term)
        return output

    @classmethod
    #def logical_one(self): return RootExpression({(): set([JordanSet({}, {1})])})
    def logical_one(cls): return RootExpression({(): {JordanSet({}, {1})}})

    @classmethod
    def logical_zero(cls): return RootExpression({})

    def __invert__(self): return self ^ RootExpression.logical_one()


    """
    Evaluation
    """

    def upper(self):
        """Upper bound on the linear complexity: the number of roots the expression allows.

        When the field degrees of a basis are pairwise coprime this is the size of
        the union of the terms' root sets (see docs/theory/Root Expressions and LC
        Estimation.md, Proposition 3 and the paragraph after Proposition 5). When
        they are not, distinct tuples can give the same product, and the box count
        overstates the term. For such a basis a warning is logged and the count
        is capped by `_capped_basis_count`: per term at the smaller of its box and
        the subgroup its fields generate. The caps hold whatever polynomials the
        fields come from, so the result stays an upper bound; how tight it is in
        between is not known.

        :return: The upper bound.
        :rtype: int
        """
        # initialize values
        linear_complexity = 0
        basis_table = {}

        # Assume any embedded sets are present:
        #  - this means no cleaning is necessary
        #  - but each embedded term might needs to be added to multiple entries in the basis table.
        #       - because the roots span across several subfields.
        for field_tuple in self.root_table:
            for term in self.root_table[field_tuple]:
                full = []           # - All bases where b == c (contains the 1 coset)
                partial = []        # - All bases where b != c (no embedded sets)

                # Add each (b,c) pair to the right list:
                for b,c in term.roots.items():
                    if b == c: full.append((b,binsum(b,c-1)))
                    else: partial.append((b,binsum(b,c)))

                # Use the full list to distribute embedded sets to appropriate entries
                #   - every subset of bases in the full list represents a valid subfield
                for modifier in powerset(full):
                    pairs = sorted(list(modifier) + partial)
                    basis = tuple([x[0] for x in pairs])
                    counts = tuple([max(term.mults)] + [x[1] for x in pairs])

                    if basis in basis_table:
                        basis_table[basis].append(counts)
                    else:
                        basis_table[basis] = [counts]

        # evaluate the basis table using hyperrec algorithm
        not_coprime = []
        for basis, x in basis_table.items():
            if len(_coprime_groups(basis)) < len(basis):
                not_coprime.append(basis)
                linear_complexity += _capped_basis_count(basis, x, pessimistic=False)
            else:
                linear_complexity += rectangle_solve(x)
        if not_coprime:
            _warn_not_coprime(not_coprime)
        return linear_complexity


    def lower(self):
        """Statistical lower estimate of the linear complexity (not a bound).

        Each class is counted by `pessimistic_expected_value`, which leaves out
        the root 1 and scales by the chance that a coset survives. For a basis
        whose degrees are not pairwise coprime, the count is capped as in
        `upper` (`_capped_basis_count`), with the subgroup's room given the same
        statistical treatment. The two estimates then use the same caps, so the
        lower estimate does not exceed the upper bound on such a basis. The
        warning is left to `upper`, which is called alongside.

        :return: The estimate.
        :rtype: float
        """
        # initalize values:
        linear_complexity = 0
        basis_table = {}

        # Evaluate only the maximum coset?
        # cleaned_terms = maximal_elements(
        #     leq_ordering = is_embedded_set,
        #     inputs = [self.terms]
        # )

        # add each term into the basis table.
        for field_tuple in self.root_table:
            for term in self.root_table[field_tuple]:
                pairs = sorted([
                    (b,pessimistic_expected_value(b,c))
                    for b,c in term.roots.items()
                ])

                basis = tuple([x[0] for x in pairs])
                counts = tuple([max(term.mults)] + [x[1] for x in pairs])

                if basis in basis_table:
                    basis_table[basis].append(counts)
                else:
                    basis_table[basis] = [counts]

        #solve basis table using hyperrec algorithm:
        for basis, rectangle_list in basis_table.items():
            if len(_coprime_groups(basis)) < len(basis):
                linear_complexity += _capped_basis_count(basis, rectangle_list, pessimistic=True)
            else:
                linear_complexity += rectangle_solve(rectangle_list)
        return linear_complexity


    """    
    ALTERNATE MULTIPLICATION (test speeds/profile as future work).

    def __mul__(self, other):
        new_root_expression = RootExpression()
        for term_a, term_b in product(self.terms, other.terms):
            new_root_expresion += (term_a * term_b)
        return new_root_expression
    """
