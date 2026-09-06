"""Gröbner basis equation store with Gebauer-Möller pair management.

Extends Buchberger's algorithm with:

- **Deferred S-polynomial computation** via (lcm, i, j) pair triples.
  S-polynomials are only formed when their pair is popped from the
  priority queue, so pairs pruned by Gebauer-Möller never incur the
  cost of polynomial multiplication.

- **Gebauer-Möller B, F, M criteria** for pair-set management.
  When a new basis element is inserted:

    - *Product criterion*: coprime leading terms are skipped
      (same as GroebnerEqStore).
    - *F criterion*: among new pairs with the same lcm, only one
      is kept.
    - *M criterion*: a new pair whose lcm is properly divisible by
      a surviving pair's lcm is dropped.
    - *B criterion*: old pairs in the queue whose lcm is divisible
      by the new element's leading term (with strictly smaller lcms
      on both sides) are pruned.

- **Generation-counter staleness detection** replaces the ``seen``
  set.  Each basis element carries a generation counter bumped on
  any modification; pairs store the generation at enqueue time and
  are discarded at pop time if either counter has advanced.

- **Interreduction**: when a new basis element's leading term
  divides an existing element's leading term, the existing element
  is reduced by the new one.  Its lead changes, stale pairs are
  invalidated via the generation counter, and fresh pairs are
  regenerated.

- **Non-constant linear-lead propagation**: when a degree-1
  polynomial x_i + p(others) enters the basis with a short tail
  (≤ ``linear_sub_threshold`` monomials), x_i is substituted
  throughout the basis — eliminating the variable even before its
  concrete value is known.
"""

import heapq
from functools import cmp_to_key

from PyPR.BooleanLogic.BooleanANF import BooleanANF
from PyPR.BooleanLogic.FunctionInputs import CONST
from PyPR.Cryptanalysis.Components.EquationStores.BaseEqStore import BaseEqStore


# ── monomial helpers (same order as GrobnerEqStore) ──────────────────

def monomial_compare(term_1, term_2):
    """Graded reverse-lex comparison on monomials (frozensets).

    :return: 0 if equal, 1 if term_1 > term_2, -1 otherwise.
    :rtype: int
    """
    if term_1 == term_2:
        return 0
    if len(term_1) > len(term_2):
        return 1
    if len(term_1) < len(term_2):
        return -1
    if sorted(term_1, reverse=True) >= sorted(term_2, reverse=True):
        return 1
    return -1

monomial_order = cmp_to_key(monomial_compare)


def lead_term(f):
    """Leading term of a BooleanANF under graded reverse-lex.

    :param f: A polynomial over GF(2).
    :type f: BooleanANF
    :return: The leading monomial (frozenset of variable indices),
        or ``None`` for the zero polynomial.
    :rtype: frozenset[int] | None
    """
    if f.terms:
        return max(f.terms, key=monomial_order)
    return None


# ── pair-queue entry ─────────────────────────────────────────────────

class _PairEntry:
    """Min-heap entry for a deferred S-polynomial pair.

    Ordered by graded monomial order on the lcm, matching the
    monomial order used for basis elements.  Generation counters
    let :meth:`GroebnerEqStore2._pop_fresh_pair` detect and skip
    stale entries without a global ``seen`` set.
    """
    __slots__ = ('lcm', 'i', 'j', 'gen_i', 'gen_j')

    def __init__(self, lcm, i, j, gen_i, gen_j):
        self.lcm = lcm
        self.i = i
        self.j = j
        self.gen_i = gen_i
        self.gen_j = gen_j

    def __lt__(self, other):
        if len(self.lcm) != len(other.lcm):
            return len(self.lcm) < len(other.lcm)
        return (
            sorted(self.lcm, reverse=True)
            <= sorted(other.lcm, reverse=True)
        )


# ── variable substitution ───────────────────────────────────────────

def _substitute_var(poly, var, tail):
    """Replace *var* with *tail* in *poly*, returning the result.

    Splits poly = cofactor·var + rest (where cofactor collects every
    monomial containing *var* with *var* factored out), then returns
    rest + cofactor·tail.

    :param poly: The polynomial to transform.
    :type poly: BooleanANF
    :param var: Variable index to replace.
    :type var: int
    :param tail: Replacement polynomial (must not contain *var*).
    :type tail: BooleanANF
    :return: The substituted polynomial.
    :rtype: BooleanANF
    """
    cofactor_terms = []
    rest_terms = []

    for term in poly.terms:
        if var in term:
            cofactor_terms.append(term - frozenset([var]))
        else:
            rest_terms.append(term)

    if not cofactor_terms:
        return poly

    cofactor = BooleanANF(frozenset(cofactor_terms), fast_init=True)
    rest = BooleanANF(frozenset(rest_terms), fast_init=True)
    return rest + cofactor * tail


# ── store ────────────────────────────────────────────────────────────

_DEFAULT_LINEAR_THRESHOLD = 4


class GroebnerEqStore2(BaseEqStore):
    """GF(2) polynomial equation store with Gebauer-Möller pair management.

    Drop-in replacement for :class:`GroebnerEqStore` with the same
    public interface (``enqueue_equation``, ``process_pending``,
    ``solved_vars``, ``unknown_vars``, ``is_determined``, ``queue``,
    ``num_eqs``, ``reduce``, ``_simplify``).

    :param simplify_mode: Simplification strategy (currently unused;
        reserved for future back-reduction modes).
    :type simplify_mode: str | None
    :param linear_sub_threshold: Maximum number of monomials in a
        linear polynomial's tail for non-constant propagation to
        fire.  Higher values eliminate more variables but risk
        polynomial blow-up.
    :type linear_sub_threshold: int
    """

    def __init__(self, simplify_mode=None,
                 linear_sub_threshold=_DEFAULT_LINEAR_THRESHOLD):
        super().__init__(consistent=False)
        self.eager = False
        self.filtering = True

        self._simplify_mode = simplify_mode
        self._linear_threshold = linear_sub_threshold

        # Basis: append-only parallel arrays.  Inactive slots stay
        # for index stability; _active[i] marks liveness.
        self._polys: list[BooleanANF] = []
        self._leads: list[frozenset[int]] = []
        self._active: list[bool] = []
        self._gen: list[int] = []

        self.num_eqs: int = 0

        # Deferred pair queue
        self._pairs: list[_PairEntry] = []

        # Input buffer: raw ANFs not yet reduced / inserted
        self._inputs: list[BooleanANF] = []

        # Variable tracking
        self.unknown_vars: set[int] = set()
        self.solved_vars: dict[int, int] = {}

    # ── compatibility properties ─────────────────────────────────

    @property
    def queue(self):
        """List whose length equals the combined input + pair queue size."""
        return [None] * (len(self._inputs) + len(self._pairs))

    @property
    def is_determined(self):
        return not (self.unknown_vars - set(self.solved_vars.keys()))

    @property
    def equations(self):
        """Active basis polynomials (read-only snapshot)."""
        return [p for p, a in zip(self._polys, self._active) if a]

    # ── public interface ─────────────────────────────────────────

    def enqueue_equation(self, equation):
        """Stage an equation for later processing.

        Solved-variable values are substituted immediately so that
        the input is as simple as possible before it enters the
        reduction pipeline.

        :param equation: A BooleanFunction to add.
        :type equation: BooleanFunction
        """
        equation = equation.compose({
            var: CONST(val) for var, val in self.solved_vars.items()
        })
        self.unknown_vars |= set(equation.idxs_used())
        anf = BooleanANF.from_BooleanFunction(equation)
        if anf.terms:
            self._inputs.append(anf)

    def insert_equation(self, equation, identifier=None,
                        translate_ANF=True):
        self.enqueue_equation(equation)
        self.process_pending()

    def queue_equation(self, equation, identifier=None,
                       translate_ANF=True):
        self.enqueue_equation(equation)

    def process_pending(self, *, verbose=False, batch_size=None,
                        _print_depth=0):
        """Process the input and pair queues.

        Pops entries from the input queue (priority) or pair queue,
        reduces them, and inserts non-zero remainders into the basis
        — which triggers Gebauer-Möller pair generation, interreduction,
        and propagation.

        :param verbose: Print live progress.
        :type verbose: bool
        :param batch_size: Maximum number of queue pops before
            returning.  ``None`` processes everything.
        :type batch_size: int | None
        :param _print_depth: Indentation level for verbose output.
        :type _print_depth: int
        :return: Number of queue entries consumed.
        :rtype: int
        """
        if batch_size is None:
            batch_size = -1          # will decrement but never reach 0

        indent1 = '|   ' * (_print_depth + 1)
        indent2 = '|   ' * (_print_depth + 2)
        printed_header = False

        if verbose and (self._inputs or self._pairs):
            print(f"{indent1}Running Groebner basis reduction (GM):")
            printed_header = True

        consumed = 0
        while (self._inputs or self._pairs) and batch_size:
            # ── pop next polynomial ──────────────────────────────
            if self._inputs:
                poly = self._inputs.pop()
            else:
                poly = self._pop_fresh_pair()
                if poly is None:
                    break

            batch_size -= 1
            consumed += 1

            # ── reduce against current basis ─────────────────────
            reduced = self.reduce(poly)
            reduced_lead = lead_term(reduced)

            if reduced_lead is None:
                if verbose:
                    self._print_status(indent2, consumed)
                continue

            if reduced.terms == frozenset([frozenset()]):
                if verbose:
                    print(
                        f"\n{indent2}Contradiction found!"
                        f"  (processed: {consumed}"
                        f"  --  basis: {self.num_eqs}"
                        f"  --  solved: {len(self.solved_vars)})"
                    )
                raise ValueError("Inconsistent")

            # ── insert, generate pairs, interreduce, propagate ───
            idx = self._insert(reduced, reduced_lead)
            self._generate_pairs_gm(idx)
            self._interreduce_by(idx)

            if len(reduced_lead) <= 1:
                self._propagate()

            if verbose:
                self._print_status(indent2, consumed)

        if printed_header:
            print()

        return consumed

    def reduce(self, poly):
        """Top-reduce *poly* against all active basis elements.

        :param poly: Polynomial to reduce.
        :type poly: BooleanANF
        :return: The remainder after iterated leading-term reduction.
        :rtype: BooleanANF
        """
        curr_lead = lead_term(poly)
        while curr_lead is not None:
            selected = None
            for i in range(len(self._polys)):
                if self._active[i] and self._leads[i] <= curr_lead:
                    selected = i
                    break
            if selected is None:
                return poly
            quotient = curr_lead - self._leads[selected]
            poly = poly + BooleanANF([quotient]) * self._polys[selected]
            curr_lead = lead_term(poly)
        return poly

    def _simplify(self):
        """Substitute solved variables into the basis and propagate.

        Called externally by SplitGrob_Solver after confirming a
        variable's value.
        """
        self._apply_solved()
        self._propagate()

    # ── basis management ─────────────────────────────────────────

    def _insert(self, poly, lt):
        """Append a reduced, non-zero polynomial to the basis.

        :return: The new element's index.
        :rtype: int
        """
        idx = len(self._polys)
        self._polys.append(poly)
        self._leads.append(lt)
        self._active.append(True)
        self._gen.append(0)
        self.num_eqs += 1
        return idx

    def _pop_fresh_pair(self):
        """Pop pair entries until a fresh one is found, compute its
        S-polynomial, and return it.

        :return: The S-polynomial, or ``None`` if the queue is
            exhausted (all remaining entries were stale).
        :rtype: BooleanANF | None
        """
        while self._pairs:
            e = heapq.heappop(self._pairs)
            if (self._active[e.i]
                    and self._active[e.j]
                    and self._gen[e.i] == e.gen_i
                    and self._gen[e.j] == e.gen_j):
                return (
                    self._polys[e.i]
                    * BooleanANF([self._leads[e.j] - self._leads[e.i]])
                    + self._polys[e.j]
                    * BooleanANF([self._leads[e.i] - self._leads[e.j]])
                )
        return None

    # ── Gebauer-Möller pair generation ───────────────────────────

    def _generate_pairs_gm(self, new_idx):
        """Generate pairs for a newly inserted element, applying
        product / F / M criteria to the new pairs and the B
        criterion to prune stale old pairs.

        :param new_idx: Basis index of the newly inserted element.
        :type new_idx: int
        """
        new_lt = self._leads[new_idx]

        # ---- candidate pairs (product criterion filters coprime) ----
        candidates = []
        for i in range(len(self._polys)):
            if i == new_idx or not self._active[i]:
                continue
            if not (new_lt & self._leads[i]):        # coprime
                continue
            candidates.append((new_lt | self._leads[i], i))

        # ---- F criterion: one pair per distinct lcm ----
        by_lcm: dict[frozenset[int], int] = {}
        for lcm, i in candidates:
            if lcm not in by_lcm:
                by_lcm[lcm] = i

        # ---- M criterion: drop pairs whose lcm is properly
        #      divisible by a surviving pair's lcm ----
        lcms_asc = sorted(by_lcm, key=len)
        surviving_lcms: list[frozenset[int]] = []
        for lcm in lcms_asc:
            if any(prev < lcm for prev in surviving_lcms):
                continue
            surviving_lcms.append(lcm)

        # ---- B criterion: prune old pairs ----
        # Remove old pair (b_i, b_j) when lt(new) divides the pair's
        # lcm and both "shortcut" lcms — lcm(new, b_i) and
        # lcm(new, b_j) — are strictly smaller than the pair's lcm.
        pruned: list[_PairEntry] = []
        for entry in self._pairs:
            if (new_lt <= entry.lcm
                    and self._active[entry.i]
                    and self._active[entry.j]):
                lcm_ni = new_lt | self._leads[entry.i]
                lcm_nj = new_lt | self._leads[entry.j]
                if lcm_ni != entry.lcm and lcm_nj != entry.lcm:
                    continue                         # pruned
            pruned.append(entry)

        # ---- push surviving new pairs, rebuild heap ----
        for lcm in surviving_lcms:
            i = by_lcm[lcm]
            pruned.append(_PairEntry(
                lcm, new_idx, i,
                self._gen[new_idx], self._gen[i],
            ))
        heapq.heapify(pruned)
        self._pairs = pruned

    def _regenerate_pairs_for(self, idx, exclude=None):
        """Push fresh pairs for a modified basis element.

        :param idx: Index of the modified element.
        :type idx: int
        :param exclude: Skip this index (avoids duplicating pairs
            already generated by :meth:`_generate_pairs_gm`).
        :type exclude: int | None
        """
        lt = self._leads[idx]
        for j in range(len(self._polys)):
            if j == idx or j == exclude or not self._active[j]:
                continue
            if not (lt & self._leads[j]):
                continue
            heapq.heappush(self._pairs, _PairEntry(
                lt | self._leads[j], idx, j,
                self._gen[idx], self._gen[j],
            ))

    # ── interreduction ───────────────────────────────────────────

    def _interreduce_by(self, new_idx):
        """Reduce existing basis elements whose leading terms are
        divisible by the new element's leading term.

        Modified elements get their generation counter bumped
        (invalidating all their queued pairs) and fresh pairs
        regenerated.

        :param new_idx: Basis index of the newly inserted element.
        :type new_idx: int
        """
        new_lt = self._leads[new_idx]
        for i in range(len(self._polys)):
            if i == new_idx or not self._active[i]:
                continue
            if not (new_lt <= self._leads[i]):
                continue

            quotient = self._leads[i] - new_lt
            reduced = (self._polys[i]
                       + BooleanANF([quotient]) * self._polys[new_idx])
            new_lead = lead_term(reduced)

            if new_lead is None:
                self._active[i] = False
                self.num_eqs -= 1
                self._gen[i] += 1
                continue

            if reduced.terms == frozenset([frozenset()]):
                raise ValueError("Inconsistent")

            old_lead = self._leads[i]
            self._polys[i] = reduced
            self._leads[i] = new_lead

            if new_lead != old_lead:
                self._gen[i] += 1
                self._regenerate_pairs_for(i, exclude=new_idx)

    # ── propagation ──────────────────────────────────────────────

    def _propagate(self):
        """Iterate constant unit-propagation and non-constant
        linear-lead propagation until a fixed point.
        """
        changed = True
        while changed:
            changed = False

            # ── constant unit propagation ────────────────────────
            new_solved = []
            for v in list(self.unknown_vars):
                reduced = self.reduce(BooleanANF([[v]]))
                if reduced == BooleanANF([True]):
                    self.solved_vars[v] = 1
                    new_solved.append(v)
                    changed = True
                elif reduced == BooleanANF([]):
                    self.solved_vars[v] = 0
                    new_solved.append(v)
                    changed = True

            for v in new_solved:
                self.unknown_vars.discard(v)

            if new_solved:
                self._apply_solved()

            # ── non-constant linear-lead propagation ─────────────
            for i in range(len(self._polys)):
                if not self._active[i] or len(self._leads[i]) != 1:
                    continue
                var = next(iter(self._leads[i]))
                if var in self.solved_vars:
                    continue

                # tail = poly - lead monomial
                tail = self._polys[i] + BooleanANF([self._leads[i]])
                if len(tail.terms) > self._linear_threshold:
                    continue
                if any(var in term for term in tail.terms):
                    continue

                for j in range(len(self._polys)):
                    if j == i or not self._active[j]:
                        continue
                    if not any(var in t for t in self._polys[j].terms):
                        continue

                    new_poly = _substitute_var(
                        self._polys[j], var, tail,
                    )
                    if new_poly.terms == self._polys[j].terms:
                        continue

                    if new_poly.terms == frozenset([frozenset()]):
                        raise ValueError("Inconsistent")

                    new_lead = lead_term(new_poly)
                    if new_lead is None:
                        self._active[j] = False
                        self.num_eqs -= 1
                        self._gen[j] += 1
                        changed = True
                        continue

                    old_lead = self._leads[j]
                    self._polys[j] = new_poly
                    self._leads[j] = new_lead

                    if new_lead != old_lead:
                        self._gen[j] += 1
                        self._regenerate_pairs_for(j)
                        changed = True

    def _apply_solved(self):
        """Substitute solved-variable values into all active basis
        elements.  Bumps generation counters for elements whose
        leading terms change and regenerates their pairs.
        """
        zero_set = frozenset(
            v for v, val in self.solved_vars.items() if val == 0
        )
        ones_set = frozenset(
            v for v, val in self.solved_vars.items() if val == 1
        )

        for i in range(len(self._polys)):
            if not self._active[i]:
                continue

            poly = self._polys[i]
            reduced_terms = [
                term - ones_set for term in poly.terms
                if not (term & zero_set)
            ]

            filtered: set[frozenset[int]] = set()
            for term in reduced_terms:
                if term in filtered:
                    filtered.remove(term)
                else:
                    filtered.add(term)

            new_poly = BooleanANF(frozenset(filtered), fast_init=True)

            if new_poly.terms == frozenset([frozenset()]):
                raise ValueError("Inconsistent")

            new_lead = lead_term(new_poly)
            if new_lead is None:
                self._active[i] = False
                self.num_eqs -= 1
                self._gen[i] += 1
                continue

            old_lead = self._leads[i]
            self._polys[i] = new_poly
            self._leads[i] = new_lead

            if new_lead != old_lead:
                self._gen[i] += 1
                self._regenerate_pairs_for(i)

    # ── verbose output ───────────────────────────────────────────

    def _print_status(self, indent, consumed):
        pending = len(self._inputs) + len(self._pairs)
        print(
            f"\r\033[K{indent}Processed: {consumed}"
            f"  --  Basis: {self.num_eqs}"
            f"  --  Queue: {pending}"
            f"  --  Solved: {len(self.solved_vars)}",
            end=''
        )
