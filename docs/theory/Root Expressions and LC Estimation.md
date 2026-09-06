# Root Expressions and Linear Complexity Estimation

## Why We Need Estimation

For a CMPR of size $n$ bits, the linear complexity can be as large as $\prod_i (2^{n_i} - 1)$ — an exponentially large number. Computing it exactly via Berlekamp-Massey requires observing at least $2L$ output bits, which is infeasible when $L \approx 2^n$. When the register is highly nonlinear, berlekamp massey (and other exact LC algorithms) are exponential in the size of the register.

The estimation algorithm computes a pair (lower, upper) bounding the true linear complexity, without running the register at all. The upper bound is **tight** (meaning it equals the true linear complexity with high probability). The lower bound is a statistical estimate that holds with high probability. Both are computed in polynomial time (in the size of the register).

## What is a Root Expression?

The linear complexity of a signal is determined by which roots appear in its exponential representation. A signal $s[t]$ over GF(2) that satisfies a linear recurrence can be written as:

$$s[t] = \sum_\alpha c_\alpha \cdot \text{Tr}(\alpha^t)$$

where the sum is over roots $\alpha$ in the algebraic closure $\overline{\mathbb{F}_2}$, $c_\alpha \in \{0,1\}$, and $\text{Tr}$ is the field trace -- the GF(2) analogue of how Binet's formula decomposes the Fibonacci numbers into a sum of exponentials at the roots of $x^2 - x - 1$ (see [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md) for the root/field background). The linear complexity equals the number of roots $\alpha$ for which $c_\alpha \neq 0$ (counted with appropriate multiplicity).

Rather than computing the exact set of roots in a signal, a **root expression** bounds which roots can possibly appear, expressed in terms of **coset classes** — sets of roots grouped by field and coset weight.

### Coset Classes

In $\mathbb{F}_{2^e}$, squaring an element $\alpha^k$ gives $\alpha^{2k}$ — a cyclic left shift of $k$'s $e$-bit binary representation (since $\alpha^{2^e} = \alpha$). The **cyclotomic coset** of $k$ is the orbit $\{k, 2k, 4k, \ldots \bmod 2^e - 1\}$ under this shift. Because cyclic left shift preserves Hamming weight, all elements of a cyclotomic coset share the same binary weight, called the **coset weight**.

The **coset class** $\langle e \cdot w \rangle$ is the set of all roots in $\mathbb{F}_{2^e}$ whose exponent belongs to a cyclotomic coset with weight at most $w$. A **root expression** is a formal algebraic expression of coset classes — a sum of products $\sum_j \prod_i \langle e_{ji} \cdot w_{ji} \rangle$, where each product represents roots from independent component fields and the sum represents the union of possibilities across terms. The upper bound on linear complexity is computed by counting the maximum number of roots the expression could represent (using the same combinatorial rules as [monomial profiles](Monomial%20Profile%20Theory.md)).

Cross-field products of coset classes have a natural geometric interpretation as rectangular **regions** in a grid whose axes are the component fields and whose coordinates are coset weights. See [Coset Regions](../conventions/Coset%20Regions.md) for diagrams and examples of this picture.

## The Five Root-Expression Propositions

The propositions below state the counting algebra used by root expressions. They are propositions about the upper-bound sets, not claims that every listed root survives in every concrete sequence.

**Proposition 1** (Counting):

$$|\langle e \cdot w \rangle| = \sum_{i=1}^{w} \binom{e}{i}$$

*Proof:* Each exponent in $\mathbb{F}_{2^e}$ can be written as an $e$-bit binary string. There are $\binom{e}{i}$ roots with exponents of coset weight exactly $i$. Summing over all weights up to $w$ gives the count. $\square$

**Proposition 2** (Same-field product): Let $\langle e \cdot w_1 \rangle, \langle e \cdot w_2 \rangle, \ldots, \langle e \cdot w_k \rangle$ be coset classes with common exponent $e$. Then:

$$\prod_{i=1}^{k} \langle e \cdot w_i \rangle = \left\langle e \cdot \sum_{i=1}^{k} w_i \right\rangle$$

where the weight saturates at $e$ (since exponents are $e$-bit strings, coset weight cannot exceed $e$).

*Proof:* When two sequences with roots from the same field are multiplied termwise, the product's root representation is:

$$c[t] = a[t]\, b[t] = \left(\sum_{i} A_i (\alpha^i)^t\right)\left(\sum_{j} B_j (\alpha^j)^t\right) = \sum_{i}\sum_{j} A_i B_j\, (\alpha^{i+j})^t$$

New exponents $i + j$ (and new cosets) may be generated, but no resulting coset has weight greater than the sum of the weights of the two input cosets [Key76], since adding two $e$-bit binary representations can increase Hamming weight by at most the sum of the individual weights. $\square$

**Proposition 3** (Cross-field product): Let $\langle e_1 \cdot w_1 \rangle, \langle e_2 \cdot w_2 \rangle, \ldots, \langle e_k \cdot w_k \rangle$ be coset classes with distinct exponents $e_i \neq e_j$ for any $i, j$. Then:

$$\left| \prod_{i=1}^{k} \langle e_i \cdot w_i \rangle \right| = \prod_{i=1}^{k} |\langle e_i \cdot w_i \rangle|$$

*Proof:* This follows from induction on $k$. A root choice from one independent component field determines no choice in any other, so the cardinalities multiply. $\square$

**Proposition 4** (Addition / XOR): For any two root expressions $E_1$ and $E_2$:

$$|E_1 + E_2| \leq |E_1 \cup E_2| = |E_1| + |E_2| - |E_1 \cap E_2|$$

*Proof:* Consider termwise addition in the root representation. For any root $\alpha^i$, the coefficient $(A_i \oplus B_i)$ is nonzero if at least one of $A_i$ or $B_i$ is nonzero and $A_i \neq B_i$. Thus a root can only be present in the sum's root expression if it already existed in one of the original sequences. The second equality follows from inclusion-exclusion. $\square$

**Proposition 5** (Intersection of single products): Let $E_1$ and $E_2$ be root expressions, each a single product of coset classes. If they have different sets of exponents, then $E_1 \cap E_2$ is empty. Otherwise, if both products share the same exponents $\{e_1, \ldots, e_k\}$, then:

$$\left( \prod_{i=1}^{k} \langle e_i \cdot w_{1,i} \rangle \right) \cap \left( \prod_{i=1}^{k} \langle e_i \cdot w_{2,i} \rangle \right) = \prod_{i=1}^{k} \langle e_i \cdot \min(w_{1,i}, w_{2,i}) \rangle$$

*Proof:* A root in the intersection must have coset weight at most $w_{1,i}$ and at most $w_{2,i}$ in each component $i$, so the weight bound is the minimum. If one expression has an extra component field, no root of that field can appear in the other expression's product, so the intersection is empty. $\square$

These propositions parallel the five core propositions of [Monomial Profile Theory](Monomial%20Profile%20Theory.md). However, root expressions additionally track multiplicities via Jordan decomposition of generalized eigenspaces — see [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md) for how Propositions 2 and 5 extend when eigenvalues have nontrivial Jordan blocks.

## How Root Expressions Are Propagated

Root expressions propagate through a CMPR block by block: each block's root expression combines its own roots with the roots inherited from upstream chaining. This section justifies why the propagation is valid and gives the algorithm.

### Propagation Through Boolean Functions

A chaining function $\mathcal{C}$ is a boolean function of the lower blocks' state bits. Since its ANF decomposes into XOR and AND operations, the root expression of $\mathcal{C}$'s output can be computed symbolically from the root expressions of its inputs:

- **XOR** combines root expressions via union (Proposition 4).
- **AND** multiplies root expressions — same-field products add coset weights (Proposition 2), cross-field products multiply cardinalities (Proposition 3).
- **CONST(1)** contributes no roots (linear complexity 0).

In code, this is `BooleanFunction.eval_ANF(root_expressions)`.

### Why the Resolvent Justifies Propagation

A block with $s$-bit update matrix $U$ driven by chaining input $\mathcal{C}(A[t])$ satisfies the recurrence $B[t+1] = UB[t] \oplus \mathcal{C}(A[t])$. Taking the D-transform (see [Resolvent Analysis](Resolvent%20Analysis.md) for the full derivation):

$$B(D) = (I \oplus UD)^{-1}\bigl(D\,\mathcal{C}(D) \oplus B[0]\bigr)$$

Expanding the inverse via the adjugate:

$$B(D) = \frac{1}{\chi_U(D)}\,\operatorname{Adj}(I \oplus UD)\bigl(D\,\mathcal{C}(D) \oplus B[0]\bigr)$$

where $\chi_U(D) = \det(I \oplus UD)$ is the characteristic polynomial of $U$. Each entry of $B(D)$ is a linear combination of entries of $\mathcal{C}(D)$, filtered through rational functions with denominator $\chi_U(D)$.

The roots of each bit's output sequence therefore come from two sources:

1. **The block's own roots** — from the denominator $\chi_U(D)$, which contributes roots at period $2^s - 1$. For a primitive polynomial of degree $s$, these form the weight-1 cyclotomic coset, giving the coset class $\langle s \cdot 1 \rangle$.
2. **The chaining input's roots** — from $\mathcal{C}(D)$ in the numerator, bounded by the root expression of $\mathcal{C}$ evaluated symbolically on the upstream blocks' root expressions.

The adjugate entries are polynomials in $D$ — they select and combine entries of $\mathcal{C}(D)$ but cannot introduce new periodicities. The denominator $\chi_U$ is fixed by the block's update matrix. So the resolvent filters cannot create roots beyond these two sources, and the root expression for each bit of the driven block is bounded by $\text{RE}(\mathcal{C}) + \langle s \cdot 1 \rangle$.

For the full resolvent derivation in both D-transform and Z-transform conventions, and the relationship between the two resolvent matrices, see [Resolvent Analysis](Resolvent%20Analysis.md). The eigenstructure of the resolvent entries feeds the root/multiplicity calculation described in [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md).

### The CMPR Root Expression Algorithm

Given a CMPR with $q$ MPR blocks $\{M_1, \ldots, M_q\}$ of sizes $\{s_1, \ldots, s_q\}$ and chaining functions $\mathcal{C}_1, \ldots, \mathcal{C}_{q-1}$:

```
procedure UPPER_BOUND(CMPR):
    RootExpressions ← table indexed by bit position

    # Block 1 (source) — no chaining input
    for b in bits(M_1):
        RootExpressions[b] ← ⟨s_1 · 1⟩

    # Each subsequent block: chaining composition + own roots
    for i ← 2 to q:
        composition ← C_{i-1}.eval_ANF(RootExpressions)
        for b in bits(M_i):
            RootExpressions[b] ← composition + ⟨s_i · 1⟩

    # Evaluate the output function on the full table
    return Evaluate(RootExpressions[output_bit])
```

Block 1 (the source) has no chaining input, so each of its bits gets $\langle s_1 \cdot 1 \rangle$. For each subsequent block $i$, the chaining function $\mathcal{C}_{i-1}$ is evaluated symbolically on the current root expression table — applying Propositions 1–5 gate by gate — and the result is summed with the block's own $\langle s_i \cdot 1 \rangle$. The algorithm processes blocks from upstream to downstream, matching CMPR's chaining order (see [Notation and Terminology](../conventions/Notation%20and%20Terminology.md)). The final `Evaluate` step counts the maximum number of roots in the output bit's root expression, giving the upper bound $\Lambda$.

**Implementation note:** Both `RootExpression` and `MonomialProfile` are computed from the update matrix's eigenstructure. For root expressions, generalized eigenspaces require the multiplicity-aware Jordan handling described in [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md); the companion matrix determines which root cosets can appear in each bit's signal.

## The estimate_LC Algorithm

`C.estimate_LC(output_bit)` works as follows:

1. Call `C.root_expressions()` to get a root expression for every bit of the register.
2. Evaluate the output function symbolically with these root expressions to get the root expression for the output bit.
3. Compute `upper = re.upper()` — the maximum number of roots the expression could represent.
4. Compute `lower = max(blockLen, re.lower())` — a statistical lower estimate, floored at the size of the block the output bit lives in (since the block's own roots always contribute).

The lower bound is a float (a statistical expectation), not an integer. The upper bound is an integer.

### When the Upper Bound is Not Tight

The upper bound is the maximum possible linear complexity if every root in the expression is actually present. In practice, some roots cancel, making the true LC smaller. The cases where the bound is loose:

- **Random cancellation:** Under the Key-Rueppel model, each cyclotomic coset of size $n$ is absent with probability approximately $1/2^n$. For most cosets this probability is negligible, but the weight-$e$ coset in $\mathbb{F}_{2^e}$ contains only 1 element (since $\binom{e}{e} = 1$), giving a cancellation probability of $1/2$ — a coin flip. The pessimistic (lower-bound) estimate handles this by capping each coset class's weight at $e - 1$, always assuming the worst case where these weight-$e$ roots degenerate. The optimistic (upper-bound) estimate includes them.

Note also that these bounds estimate *trajectory* LC, not *column* LC in an equation matrix -- the two have a counterintuitive [inverse relationship across CMPR block structure](Trajectory%20and%20Column%20LC.md) that means RE upper bounds should not be used to predict column LC.

In the 10-trial diagnostic on `[M7, M5, M3, M2]` with `old_ANF_template()`, the upper bound was 13023 and the lower was ~10090. 8/10 trials landed in range; 2 trials had actual LC ~3930, well below the lower bound. These failures correspond to unlucky random chaining functions that cause significant cancellation — a degenerate case rather than a bug in the bounds.

## Locked Blocks

The `locked_list` parameter controls contributions from blocks. This is used in cryptanalytic scenarios where some component registers are already known to the attacker ("locked") and do not contribute to the actual security. Pass one truthy/falsey flag per block; a truthy flag leaves that block's roots active, while a falsey flag suppresses its root contribution and gives tighter bounds for the remaining unknown components.

This feature is noted as "slightly under-tested" in the notebook — use with caution.

## T-Functions

Because `TFunction` is a subclass of `CMPR`, `estimate_LC` technically works on T-Functions too. However, there are two problems:
- T-Functions have many components (one per bit), so the $O(nk^2 N^4)$ complexity becomes impractical quickly.
- Repeated components (many bits share the same underlying polynomial) disable the [mesh optimization](../architecture/Mesh%20Optimization.md), which requires distinct block sizes.

The function handles small T-Functions embedded in larger CMPRs reasonably well, but was not designed for them. For T-functions specifically, the period and complexity structure is better understood through the multiplicity bracket framework — see [Root Multiplicities — Multiplicity Brackets and the Period Staircase](Root%20Multiplicities%20and%20Jordan%20Decomposition.md#multiplicity-brackets-and-the-period-staircase).

## Runtime Complexity of Root Counting

Let $n$ = total bits, $k$ = number of MPR components, $S$ = set of component sizes. Define:

$$N \leq \frac{\prod_{s_i \in S} s_i}{\max(S)}$$

as the maximum number of terms in any root expression. This is the key quantity — it bounds how large the intermediate objects get.

**Cost per addition of root expressions:** $O(kN^2)$
**Cost per multiplication of root expressions:** $O(kN^3)$

*Derivation:* Building the intermediate term list costs $O(kN)$ for addition, $O(kN^2)$ for multiplication (pairs of terms, each costing $O(k)$). Pruning to maximal elements requires comparing each intermediate term against the current maximal list (size $\leq N$), each comparison costing $O(k)$ — giving $O(kN^2)$ for the pruning pass after addition and $O(kN^3)$ for multiplication.

**Number of operations per block:** The chaining function's monomial profile (which has at most $N$ terms, each involving at most $n$ variables) guides the computation. There are at most $O(N)$ additions and $O(nN)$ multiplications. Since multiplications dominate:

$$\text{cost per block} = O(nN) \cdot O(kN^3) = O(knN^4)$$

**Total cost** (summed over $k$ blocks):

$$O(nk^2 N^4)$$

Substituting the bound $N \leq \frac{1}{\max(S)} \left(\frac{n}{k}\right)^k$ and using $\max(S) \geq n/k$:

$$O\!\left(\frac{k^6}{n^3} \left(\frac{n}{k}\right)^{4k}\right)$$

**Comparison to BM:** Berlekamp-Massey needs $O(L \log L)$ operations where $L = $ linear complexity $\approx 2^n$. So BM costs roughly $O(n \cdot 2^n)$.

The root counting algorithm is **polynomial in $n$** (for fixed $k$), while BM is exponential. For large $n$, root counting wins by an enormous margin — even though the polynomial can be severe (e.g., $O(n^{32})$ for $k=8$). The dominant cost driver is $k$: adding more components to a CMPR makes estimation massively more expensive (the $4k$ in the exponent).

When block sizes are all distinct and chaining is simple, the [mesh optimization](../architecture/Mesh%20Optimization.md) bypasses gate-by-gate composition entirely by directly enumerating maximal weight distribution vectors (a multinomial expansion over the component fields). This avoids building intermediate expression objects and can be orders of magnitude faster than the default path.

## API

```python
lower, upper = C.estimate_LC(output_bit=0)
REs = C.root_expressions()                    # one RootExpression per bit
filtered_re = filter_fn.eval_ANF(REs)         # propagate through an output function
```

Berlekamp-Massey (for small enough sequences):
```python
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey, berlekamp_massey_iterator
lc, poly = berlekamp_massey(seq)              # full BM, requires len(seq) >= 2*lc
for lc, poly in berlekamp_massey_iterator(seq_generator, yield_rate=1000):
    pass  # get running lower bound from a streaming iterator
```

Note: `berlekamp_massey` outputs the **dual** polynomial (i.e. `convolve(poly, seq) = 0`), consistent with the Fibonacci/Galois LFSR constructor convention. See [Polynomial Conventions](../conventions/Polynomial%20Conventions.md) for the primal/dual distinction.
