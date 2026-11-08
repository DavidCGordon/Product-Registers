# Root Expressions and Linear Complexity Estimation

This document describes how the library bounds the linear complexity of a CMPR output sequence without running the register. A **CMPR** is a composite MPR, built by chaining MPRs (Mersenne product registers) together. Each MPR block is a linear register over $\mathbb{F}_2$; every block after the first is also driven by a **chaining function**, a Boolean function of upstream blocks' state bits that is added into its update. The **linear complexity** (LC) of a sequence is the length of the shortest LFSR that generates it.

The document defines root expressions and coset classes, states the five propositions that govern their arithmetic, shows how root expressions propagate through the blocks of a CMPR, and then covers the `estimate_LC` method, its limits, and its cost. It assumes the root and field background of [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md): the algebraic closure $\overline{\mathbb{F}_2}$, Frobenius orbits, and the field trace.

## Why We Need Estimation

For a CMPR of size $n$ bits with $k$ blocks of sizes $s_1, \ldots, s_k$, the linear complexity can be as large as $\prod_i (2^{s_i} - 1)$ — an exponentially large number. Computing it exactly via Berlekamp-Massey requires observing at least $2L$ output bits, which is infeasible when $L \approx 2^n$. When the register is highly nonlinear, Berlekamp-Massey (and other exact LC algorithms) take time exponential in the size of the register.

The estimation algorithm computes a pair (lower, upper) without running the register at all. The upper value counts every root the analysis allows, so it is the linear complexity the sequence would have if none of those roots cancelled. The lower value is a statistical estimate, not a bound. The one recorded comparison with measured values, described under [When the Upper Bound is Not Tight](#when-the-upper-bound-is-not-tight), had eight of ten measurements between the two values and two below the lower one. Both values are computed in time polynomial in $n$ for a fixed number of blocks $k$.

## What is a Root Expression?

The linear complexity of a signal is determined by which roots appear in its exponential representation. Let $s[t]$ be a signal over $\mathbb{F}_2$ that satisfies a linear recurrence whose characteristic polynomial has no repeated roots. Then it can be written as

$$s[t] = \sum_\alpha A_\alpha\, \alpha^t$$

where the sum is over finitely many distinct roots $\alpha$ in the algebraic closure $\overline{\mathbb{F}_2}$, with coefficients $A_\alpha \in \overline{\mathbb{F}_2}$. This is the $\mathbb{F}_2$ counterpart of Binet's formula, which writes the Fibonacci numbers as a sum of exponentials at the roots of $x^2 - x - 1$: the roots $\alpha$ correspond to the two roots of that polynomial, and the $A_\alpha$ to its two constants.

Because $s$ takes values in $\mathbb{F}_2$, the coefficients satisfy $A_{\alpha^2} = A_\alpha^2$ (derived in [Roots — The Frobenius Automorphism and Conjugates](Roots%20and%20the%20Algebraic%20Closure.md#the-frobenius-automorphism-and-conjugates)). If the Frobenius orbit of $\alpha$ has $d$ elements, then $\alpha^{2^d} = \alpha$ gives $A_\alpha^{2^d} = A_\alpha$, so $A_\alpha \in \mathbb{F}_{2^d}$, and the terms of one orbit add up to a trace:

$$\sum_{j=0}^{d-1} A_\alpha^{2^j} \left(\alpha^{2^j}\right)^t = \text{Tr}_{\mathbb{F}_{2^d}/\mathbb{F}_2}\!\left(A_\alpha\, \alpha^t\right).$$

So $s[t]$ is a sum, over Frobenius orbits, of $\text{Tr}_{\mathbb{F}_{2^d}/\mathbb{F}_2}(A_\alpha \alpha^t)$ with one representative $\alpha$ per orbit and the coefficient $A_\alpha \in \mathbb{F}_{2^d}$ inside the trace. The coefficient cannot be pulled out as a bit: for $\alpha$ primitive in $\mathbb{F}_{2^3}$, the sequence $\text{Tr}(\alpha \cdot \alpha^t) = \text{Tr}(\alpha^{t+1})$ is a shift of $\text{Tr}(\alpha^t)$, neither zero nor equal to it.

The linear complexity equals the number of roots $\alpha$ with $A_\alpha \neq 0$. To see this, write $E$ for the shift $(Es)[t] = s[t+1]$. For a polynomial $g \in \mathbb{F}_2[E]$,

$$g(E)\,s = \sum_\alpha A_\alpha\, g(\alpha)\, \alpha^t,$$

and by the linear independence of the sequences $t \mapsto \alpha^t$ (the Vandermonde argument used in that derivation) this is zero if and only if $g(\alpha) = 0$ for every $\alpha$ with $A_\alpha \neq 0$. So the minimal polynomial annihilating $s$ is the product of $(E - \alpha)$ over those $\alpha$, and its degree is their number. The formula above is the multiplicity-free case. When the recurrence has repeated roots, terms of the form $\binom{t}{r}\alpha^t$ appear as well; those are the subject of [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md).

Rather than computing the exact set of roots in a signal, a **root expression** bounds which roots can possibly appear, expressed in terms of **coset classes** — sets of roots grouped by field and coset weight.

### Coset Classes

Fix a primitive element $\alpha$ of $\mathbb{F}_{2^e}$, so that every nonzero element is $\alpha^k$ for exactly one exponent $k$ in $\{1, \ldots, 2^e - 1\}$. Exponents are always represented in this range. In particular the identity is $1 = \alpha^{2^e - 1}$, whose exponent is the all-ones $e$-bit string, of Hamming weight $e$; the exponent 0 is not used.

Squaring an element $\alpha^k$ gives $\alpha^{2k}$ — a cyclic left shift of $k$'s $e$-bit binary representation (since $\alpha^{2^e} = \alpha$). The **cyclotomic coset** of $k$ is the orbit $\{k, 2k, 4k, \ldots \bmod 2^e - 1\}$ under this shift. Because cyclic left shift preserves Hamming weight, all elements of a cyclotomic coset share the same binary weight, called the **coset weight**. With exponents in $\{1, \ldots, 2^e - 1\}$, coset weights run from 1 to $e$.

The **coset class** $\langle e \cdot w \rangle$ is the set of all roots in $\mathbb{F}_{2^e}$ whose exponent belongs to a cyclotomic coset with weight from 1 to $w$. A **root expression** is a formal algebraic expression of coset classes — a sum of products $\sum_j \prod_i \langle e_{ji} \cdot w_{ji} \rangle$, where each product represents roots from independent component fields and the sum represents the union of possibilities across terms. The upper bound on linear complexity is computed by counting the maximum number of roots the expression could represent (using the same combinatorial rules as [monomial profiles](Monomial%20Profile%20Theory.md)).

Cross-field products of coset classes have a geometric interpretation as rectangular **regions** in a grid whose axes are the component fields and whose coordinates are coset weights. See [Coset Regions](../conventions/Coset%20Regions.md) for diagrams and examples of this picture.

## The Five Root-Expression Propositions

The propositions below state the counting algebra used by root expressions. They are propositions about the upper-bound sets, not claims that every listed root survives in every concrete sequence.

**Proposition 1** (Counting):

$$|\langle e \cdot w \rangle| = \sum_{i=1}^{w} \binom{e}{i}$$

*Proof:* Each exponent in $\mathbb{F}_{2^e}$ can be written as an $e$-bit binary string. There are $\binom{e}{i}$ roots with exponents of coset weight exactly $i$. Summing over all weights up to $w$ gives the count. $\square$

**Proposition 2** (Same-field product): Let $\langle e \cdot w_1 \rangle, \langle e \cdot w_2 \rangle, \ldots, \langle e \cdot w_k \rangle$ be coset classes with common exponent $e$. Then every product of one root from each class lies in the class whose weight is the sum:

$$\prod_{i=1}^{k} \langle e \cdot w_i \rangle \subseteq \left\langle e \cdot \sum_{i=1}^{k} w_i \right\rangle$$

where the weight saturates at $e$ (since exponents are $e$-bit strings, coset weight cannot exceed $e$). Root-expression arithmetic uses the right-hand side as the product of the classes. It is an upper bound on the set of products, like every root expression.

*Proof:* When two sequences with roots from the same field are multiplied termwise, the product's root representation is:

$$c[t] = a[t]\, b[t] = \left(\sum_{i} A_i (\alpha^i)^t\right)\left(\sum_{j} B_j (\alpha^j)^t\right) = \sum_{i}\sum_{j} A_i B_j\, (\alpha^{i+j})^t$$

New exponents $i + j$ (and new cosets) may be generated, but no resulting coset has weight greater than the sum of the weights of the two input cosets [Key76]. Writing $\text{wt}$ for Hamming weight, the integer sum satisfies $\text{wt}(i + j) \leq \text{wt}(i) + \text{wt}(j)$, because carries only remove ones. If $i + j \geq 2^e$, write $i + j = 2^e + r$; reduction modulo $2^e - 1$ replaces it by $r + 1$, and $\text{wt}(r + 1) \leq \text{wt}(r) + 1 = \text{wt}(2^e + r)$. So the reduced exponent has weight at most $\text{wt}(i) + \text{wt}(j)$. This is the case $k = 2$; the general case follows by applying it $k - 1$ times. $\square$

**Proposition 3** (Cross-field product): Let $\langle e_1 \cdot w_1 \rangle, \langle e_2 \cdot w_2 \rangle, \ldots, \langle e_k \cdot w_k \rangle$ be coset classes whose exponents are pairwise coprime, $\gcd(e_i, e_j) = 1$ for $i \neq j$. Then:

$$\left| \prod_{i=1}^{k} \langle e_i \cdot w_i \rangle \right| = \prod_{i=1}^{k} |\langle e_i \cdot w_i \rangle|$$

*Proof:* Take $k = 2$ first. Since $\gcd(e_1, e_2) = 1$, the fields meet only in the base field: $\mathbb{F}_{2^{e_1}} \cap \mathbb{F}_{2^{e_2}} = \mathbb{F}_2$. Suppose $\beta_1 \beta_2 = \beta_1' \beta_2'$ with $\beta_1, \beta_1' \in \mathbb{F}_{2^{e_1}}$ and $\beta_2, \beta_2' \in \mathbb{F}_{2^{e_2}}$ nonzero. Then $\beta_1/\beta_1' = \beta_2'/\beta_2$ lies in both fields, hence in $\mathbb{F}_2$, and being nonzero it is 1. So $\beta_1 = \beta_1'$ and $\beta_2 = \beta_2'$: the product map is injective, and the cardinalities multiply. For general $k$, induct: the products of the first $k - 1$ factors lie in $\mathbb{F}_{2^L}$ with $L = \text{lcm}(e_1, \ldots, e_{k-1})$, and $\gcd(e_k, L) = 1$, so the same argument applies to that set and the last class. $\square$

Coprimality is needed; distinct exponents are not enough. Take $e = (2, 4)$ and $w = (1, 2)$, with exponents written relative to a primitive element of $\mathbb{F}_{2^4}$. Then $\langle 2 \cdot 1 \rangle = \{5, 10\}$ and $\langle 4 \cdot 2 \rangle = S = \{1, 2, 4, 8, 3, 6, 12, 9, 5, 10\}$. Products add exponents modulo 15, and $5 + 6 = 11 = 10 + 1$, so two different pairs give the same product. The product set $(S + 5) \cup (S + 10)$ has 15 elements, not $2 \cdot 10 = 20$. Without coprimality the product of the counts is still an upper bound on the size of the product set, since every element of the product set comes from at least one tuple.

**Proposition 4** (Addition / XOR): For any two root expressions $E_1$ and $E_2$:

$$|E_1 + E_2| \leq |E_1 \cup E_2| = |E_1| + |E_2| - |E_1 \cap E_2|$$

*Proof:* Consider termwise addition in the root representation. For any root $\alpha^i$, the coefficient $(A_i \oplus B_i)$ is nonzero if at least one of $A_i$ or $B_i$ is nonzero and $A_i \neq B_i$. Thus a root can only be present in the sum's root expression if it already existed in one of the original sequences. The second equality follows from inclusion-exclusion. $\square$

**Proposition 5** (Intersection of single products): Let $E_1 = \prod_{i \in A} \langle e_i \cdot w_i \rangle$ and $E_2 = \prod_{i \in B} \langle e_i \cdot v_i \rangle$ be single products of coset classes, where the exponents $e_i$ for $i \in A \cup B$ are pairwise coprime. Then $E_1 \cap E_2$ is nonempty exactly when $w_i = e_i$ for every $i \in A \setminus B$ and $v_i = e_i$ for every $i \in B \setminus A$, and in that case

$$E_1 \cap E_2 = \prod_{i \in A \cap B} \langle e_i \cdot \min(w_i, v_i) \rangle ,$$

where the empty product is $\{1\}$.

*Proof:* By the argument of Proposition 3, applied to the fields indexed by $A \cup B$, the map $(\beta_i)_{i \in A \cup B} \mapsto \prod_i \beta_i$ on tuples of nonzero elements $\beta_i \in \mathbb{F}_{2^{e_i}}$ is injective. An element of $E_1$ is the image of a tuple with $\beta_i \in \langle e_i \cdot w_i \rangle$ for $i \in A$ and $\beta_i = 1$ for $i \in B \setminus A$; an element of $E_2$ is the image of a tuple with $\beta_i \in \langle e_i \cdot v_i \rangle$ for $i \in B$ and $\beta_i = 1$ for $i \in A \setminus B$. By injectivity, a common element is the image of a single tuple satisfying both descriptions. For $i \in A \setminus B$ that requires $1 \in \langle e_i \cdot w_i \rangle$; the exponent of $1$ is the all-ones string, of weight $e_i$, and $w_i \le e_i$, so this holds exactly when $w_i = e_i$. Symmetrically, $v_i = e_i$ for $i \in B \setminus A$. For $i \in A \cap B$ it requires $\beta_i$ in both classes, and since $\langle e \cdot w \rangle \subseteq \langle e \cdot w' \rangle$ for $w \le w'$, that is $\beta_i \in \langle e_i \cdot \min(w_i, v_i) \rangle$. Conversely, under these conditions every tuple with $\beta_i \in \langle e_i \cdot \min(w_i, v_i) \rangle$ for $i \in A \cap B$ and $\beta_i = 1$ elsewhere has its image in both products. $\square$

The full-class condition in the conclusion is needed: products over different fields can meet, because every full class $\langle e \cdot e \rangle$ contains $1$, so $\langle 3 \cdot 3 \rangle \cap \langle 5 \cdot 5 \rangle = \{1\}$. The coprimality hypothesis is needed: without it the minimum formula fails. With exponents written modulo 15 relative to a primitive element $\alpha$ of $\mathbb{F}_{2^4}$, the product $\langle 2 \cdot 2 \rangle \times \langle 4 \cdot 1 \rangle$ contains $1 \cdot \alpha = \alpha$, and $\langle 2 \cdot 1 \rangle \times \langle 4 \cdot 2 \rangle$ contains $\alpha^{10} \cdot \alpha^{6} = \alpha$, but $\langle 2 \cdot 1 \rangle \times \langle 4 \cdot 1 \rangle$ has exponents $\{3, 6, 7, 9, 11, 12, 13, 14\}$ and does not contain $\alpha$.

The monomial-profile Proposition 5 ([Monomial Profile Theory](Monomial%20Profile%20Theory.md)) needs neither hypothesis: every monomial-profile pair draws at least one variable from its register, so no pair contains the empty monomial, which plays the role that $1$ plays here.

`RootExpression.upper()` does not rely on this proposition. It splits each full class $\langle e \cdot e \rangle$ of a term into its elements other than $1$ and $\{1\}$, over every subset of the term's full fields (the basis-table loop in `RootExpression.upper()`), and files each piece under the set of fields where its elements differ from $1$. When the field degrees are pairwise coprime, pieces filed under different sets are disjoint, by the injectivity above; under one set, every piece is a box anchored at the origin, and `rectangle_solve` measures their union. So `upper()` is the size of the union of the terms' root sets, each root counted once for every Jordan length up to the largest that a term covering it records. When the degrees filed under one set are not pairwise coprime, the box overstates the term, as the example after Proposition 3 shows. For such a set `upper()` logs a warning and caps the count. Fields that share a factor are grouped, and the products over a group lie in the subgroup of $\mathbb{F}_{2^L}^\times$ ($L$ the lcm of the degrees) generated by the fields' multiplicative groups. That subgroup has order $\operatorname{lcm}_i(2^{e_i} - 1)$, so each term is counted as the product, over groups, of the smaller of its box and that order. The count of the set is the least of this sum over terms, the box union, and the largest Jordan length times the product of the orders. These caps hold whatever polynomials the fields come from. A tighter cap through exponent weights would need the fields' generators to be compatible across the embeddings, which a CMPR's independently chosen blocks are not. On the example above the cap gives the exact count, 15. `lower()` applies the same caps, with each group's order given the statistical treatment of `pessimistic_expected_value`: the root $1$ is left out and the count is scaled by $(1 - 2^{-L})^2$. So the lower estimate does not exceed the upper bound on these sets. The cap is checked against Berlekamp–Massey on CMPRs with non-coprime block sizes (`test/Tools/test_root_expression_noncoprime.py`). Such products arise only from the third block on, since a block's chaining multiplies bits of the block before it. In those cases the true linear complexity is often far below both values, so `lower()` is not a bound there either.

These propositions parallel the five core propositions of [Monomial Profile Theory](Monomial%20Profile%20Theory.md). However, root expressions additionally track multiplicities via Jordan decomposition of generalized eigenspaces — see [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md) for how Propositions 2 and 5 extend when eigenvalues have nontrivial Jordan blocks.

## How Root Expressions Are Propagated

Root expressions propagate through a CMPR block by block: each block's root expression combines its own roots with the roots inherited from upstream chaining. This section justifies why the propagation is valid and gives the algorithm.

### Propagation Through Boolean Functions

A chaining function $\mathcal{C}$ is a Boolean function of the upstream blocks' state bits. Since its ANF decomposes into XOR and AND operations, the root expression of $\mathcal{C}$'s output can be computed symbolically from the root expressions of its inputs:

- **XOR** combines root expressions via union (Proposition 4).
- **AND** multiplies root expressions — same-field products add coset weights (Proposition 2), cross-field products multiply cardinalities (Proposition 3).
- **CONST(1)** contributes the single root $\alpha = 1$. The all-ones sequence is nonzero and satisfies $s[t+1] = s[t]$, so its linear complexity is 1. In code it is `RootExpression.logical_one()`, a term with no field component and Jordan length set $\{1\}$.

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

Given a CMPR with $k$ MPR blocks $\{M_1, \ldots, M_k\}$ of sizes $\{s_1, \ldots, s_k\}$ and chaining functions $\mathcal{C}_1, \ldots, \mathcal{C}_{k-1}$:

```
procedure UPPER_BOUND(CMPR):
    RootExpressions ← table indexed by bit position

    # Block 1 (source) — no chaining input
    for b in bits(M_1):
        RootExpressions[b] ← ⟨s_1 · 1⟩

    # Each subsequent block: chaining composition + own roots
    for i ← 2 to k:
        composition ← C_{i-1}.eval_ANF(RootExpressions)
        for b in bits(M_i):
            RootExpressions[b] ← composition + ⟨s_i · 1⟩

    # Evaluate the output function on the full table
    return Evaluate(RootExpressions[output_bit])
```

Block 1 (the source) has no chaining input, so each of its bits gets $\langle s_1 \cdot 1 \rangle$. For each subsequent block $i$, the chaining function $\mathcal{C}_{i-1}$ is evaluated symbolically on the current root expression table — applying Propositions 1–5 gate by gate — and the result is summed with the block's own $\langle s_i \cdot 1 \rangle$. The algorithm processes blocks from upstream to downstream, matching CMPR's chaining order (see [Notation and Terminology](../conventions/Notation%20and%20Terminology.md)). The final `Evaluate` step counts the maximum number of roots in the output bit's root expression, giving the upper bound $\Lambda$.

In the implementation (`CMPR.root_expressions` in `src/PyPR/FeedbackFunctions/CMPR.py`), the step "composition + $\langle s_i \cdot 1 \rangle$" is `RootExpression.extend`. Besides adding the block's own class, `extend` records a Jordan length one greater than the largest Jordan length of any term already present that covers the block's roots. Jordan lengths, and the multiplicity-aware handling of same-field products, are described in [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md).

## The estimate_LC Algorithm

`C.estimate_LC(output_bit)` works as follows:

1. Call `C.root_expressions()` to get a root expression for every bit of the register.
2. Evaluate the output function symbolically with these root expressions to get the root expression for the output bit.
3. Compute `upper = re.upper()` — the maximum number of roots the expression could represent.
4. Compute `lower = max(block_len, re.lower())` — a statistical lower estimate, floored at the size of the block the output bit lives in (since the block's own roots contribute). When that block is locked (see below) its own roots are not added and there is no floor.

The upper value is an integer. The lower value is in general not an integer, because of how `re.lower()` weights each class. For a class $\langle b \cdot c \rangle$, `lower()` uses `pessimistic_expected_value(b, c)` from `src/PyPR/Tools/RootCounting/EvaluationFunctions.py` in place of the count of Proposition 1:

$$\left(\sum_{i=1}^{c'} \binom{b}{i}\right) \cdot \left(\frac{2^b - 1}{2^b}\right)^{2}, \qquad c' = \begin{cases} c - 1 & \text{if } c = b, \\ c & \text{otherwise.} \end{cases}$$

The first factor is the count of Proposition 1 with the weight capped at $b - 1$, which excludes the weight-$b$ root discussed below.

> **OPEN:** What does the factor $\left((2^b - 1)/2^b\right)^2$ in `pessimistic_expected_value` model? — needs author input.
> The code multiplies each class count by it (`EvaluationFunctions.py`, lines 15–20). No document in the repository derives it, and the statistical model described below gives a per-coset absence probability, not this factor.

### When the Upper Bound is Not Tight

The upper bound is the maximum possible linear complexity if every root in the expression is actually present. In practice, some roots cancel, making the true LC smaller.

**Random cancellation.** The statistical model used here is the one described in [Monomial Profile Theory — Statistical Model](Monomial%20Profile%20Theory.md#statistical-model-for-expected-behavior): roots cancel a whole Frobenius coset at a time, and, assuming statistical independence, a coset of size $n$ is absent with probability approximately $1/2^n$. For most cosets this probability is negligible. The exception is the weight-$e$ coset in $\mathbb{F}_{2^e}$, which is the single root $1$ (exponent $2^e - 1$, and $\binom{e}{e} = 1$). Its coefficient lies in $\mathbb{F}_2$, so under the model it is absent with probability $1/2$. The pessimistic (lower) estimate handles this by capping each coset class's weight at $e - 1$, always assuming the worst case where the weight-$e$ root is absent. The optimistic (upper-bound) estimate includes it.

**A recorded comparison with measured values.** In a ten-trial diagnostic on a CMPR with blocks of sizes 7, 5, 3, and 2, the chaining functions were drawn at random with `old_ANF_template()` (`src/PyPR/BooleanLogic/ChainingGeneration/Templates/old_ANF_template.py`), a chaining template. For blocks of more than one bit, as here, that template samples a function for each bit of every block after the first in `cmpr.blocks`: an XOR of between 2 and `max_xor` terms, each term an AND of variables drawn from the preceding block, with the number of variables per term drawn from 2 to `max_and` and capped below that block's size (both parameters default to 4). The upper bound was 13023 and the lower estimate was approximately 10090. Eight of the ten trials had measured LC inside this range; the other two had measured LC of approximately 3930, well below the lower estimate. This is the only comparison recorded in this documentation, and it does not show how often the measured value equals the upper bound.

> **OPEN:** Can this diagnostic be reproduced, and what explains the two low trials? — needs author input.
> No script or seed for it exists in the repository, and neither the output bit that was measured nor the template parameters used are stated. An earlier version of this passage attributed the two low trials to randomly drawn chaining functions that cause significant cancellation; nothing in the repository supports or refutes that attribution.

**Trajectory LC, not column LC.** These bounds estimate *trajectory* LC, not *column* LC in an equation matrix. The two have an [inverse relationship across CMPR block structure](Trajectory%20and%20Column%20LC.md), so root-expression upper bounds should not be used to predict column LC.

## Locked Blocks

The `locked_list` parameter of `root_expressions` and `estimate_LC` selects which blocks contribute their own roots. This is used in cryptanalytic scenarios where some component registers are already known to the attacker ("locked") and do not contribute to the actual security. Pass one truthy/falsey flag per block; a truthy flag leaves that block's roots active, while a falsey flag suppresses its root contribution and gives tighter bounds for the remaining unknown components.

The tests in `test/Tools/linear_complexity_test.py` check that suppressing every block does not raise the upper bound, that the mesh and default algorithms agree for every pattern of locked blocks on a three-block register, and that locking a bit's own block removes the block-size floor. None of them compares a bound computed with a nontrivial `locked_list` against a measured linear complexity, so such bounds should be checked against measured values before they are relied on.

## T-Functions

Because `TFunction` is a subclass of `CMPR`, `estimate_LC` can be called on T-Functions too. There are two obstacles:
- T-Functions have many components (one per bit), so the $O(nk^2 N^4)$ complexity (derived under [Runtime Complexity of Root Counting](#runtime-complexity-of-root-counting) below) becomes impractical quickly.
- The [mesh optimization](../architecture/Mesh%20Optimization.md) for root expressions requires pairwise coprime block sizes and no 1-bit block. Every T-function block has one bit, so the 1-bit condition excludes them (repeated blocks of size 1 do not violate coprimality, since $\gcd(1, 1) = 1$).

The method was not designed for T-Functions, though small T-Function segments embedded in a larger CMPR remain within its practical range. For T-functions specifically, the period and complexity structure is better understood through the multiplicity bracket framework — see [Root Multiplicities — Multiplicity Brackets and the Period Staircase](Root%20Multiplicities%20and%20Jordan%20Decomposition.md#multiplicity-brackets-and-the-period-staircase).

## Runtime Complexity of Root Counting

Let $n$ be the total number of bits, $k$ the number of MPR blocks, and $S = \{s_1, \ldots, s_k\}$ the set of block sizes. Let $N$ be the maximum number of terms in any root expression arising in the computation. Every cost below is expressed in terms of $N$, because it bounds how large the intermediate objects get. The analysis uses the bound

$$N \leq \frac{\prod_{s_i \in S} s_i}{\max(S)}.$$

> **OPEN:** What is the proof of the bound $N \leq \prod_{s_i \in S} s_i / \max(S)$? — needs author input.
> The bound is used below to reach the closed-form cost. No derivation of it appears in this document or elsewhere in the repository.

**Cost per operation.** One addition of root expressions costs $O(kN^2)$, and one multiplication costs $O(kN^3)$.

*Derivation:* Building the intermediate term list costs $O(kN)$ for addition, $O(kN^2)$ for multiplication (pairs of terms, each costing $O(k)$). Pruning to maximal elements requires comparing each intermediate term against the current maximal list (size $\leq N$), each comparison costing $O(k)$ — giving $O(kN^2)$ for the pruning pass after addition and $O(kN^3)$ for multiplication.

**Number of operations per block:** The chaining function's monomial profile (which has at most $N$ terms, each involving at most $n$ variables) guides the computation. There are at most $O(N)$ additions and $O(nN)$ multiplications. Since multiplications dominate:

$$\text{cost per block} = O(nN) \cdot O(kN^3) = O(knN^4)$$

**Total cost** (summed over $k$ blocks):

$$O(nk^2 N^4)$$

Substituting the bound $N \leq \frac{1}{\max(S)} \left(\frac{n}{k}\right)^k$ and using $\max(S) \geq n/k$:

$$O\!\left(\frac{k^6}{n^3} \left(\frac{n}{k}\right)^{4k}\right)$$

**Comparison to Berlekamp-Massey (BM).** BM needs at least $2L$ output bits, where $L$ is the linear complexity. With $L \approx 2^n$ its running time is therefore exponential in $n$.

The root counting algorithm is **polynomial in $n$** (for fixed $k$), while BM is exponential, so for large $n$ root counting is the cheaper of the two — even though the degree of the polynomial can be high (e.g., $O(n^{29})$ for $k=8$, from $n^{4k - 3}$). The dominant cost driver is $k$: each component added to a CMPR raises that degree, through the $4k$ in the exponent.

When block sizes are pairwise coprime and chaining is simple, the [mesh optimization](../architecture/Mesh%20Optimization.md) bypasses gate-by-gate composition entirely by directly enumerating maximal weight distribution vectors: the terms of the nested expansion $\text{RE}_i = (\text{RE}_{i-1})^{d_i} + \langle s_i \cdot 1 \rangle$, where $d_i$ is block $i$'s chaining degree. This avoids building intermediate expression objects and can be orders of magnitude faster than the default path.

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
