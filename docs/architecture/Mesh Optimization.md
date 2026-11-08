# Mesh Optimization

The **mesh optimization** is an alternative algorithm for computing the root expressions and monomial profiles of a CMPR, a chain of MPR blocks in which every block after the first (the **source**) also receives in its update a **chaining function** of upstream blocks' bits. It applies when the chaining is simple, meaning that each block's chaining function refers only to the immediately preceding block ([The Model Hypothesis](#the-model-hypothesis) states this precisely); root expressions additionally require that no block have one bit and that the block sizes be pairwise coprime (see [Eligibility Conditions](#eligibility-conditions)). Under simple chaining, the root expression of each block is a power of the preceding block's root expression plus the block's own roots, so the expressions nest one inside the next ([The Nested Power Model](#the-nested-power-model)). The default path composes root expressions gate by gate through each chaining function's ANF. The mesh instead enumerates the terms of the nested expansion directly, as vectors of weights over the blocks, and keeps only the maximal ones ([Correctness of the Iterator](#correctness-of-the-iterator)).

The optimization lives in `src/PyPR/Tools/RootCounting/MeshOptimization.py` and is invoked automatically by `CMPR.root_expressions()` and `CMPR.monomial_profiles()` when the eligibility conditions are met. Pass `force_default=True` to bypass it.

## The Nested Power Model

Number the blocks $0, 1, \ldots, m$ in chaining order: block 0 is the source, and block $m$ is the **driven** block, the one whose expression is being computed. Block $j$ has $s_j$ bits. This is the order of `CMPR.blocks`, upstream to downstream ([Notation and Terminology](../conventions/Notation%20and%20Terminology.md)). Write $\text{RE}_j$ for the root expression that the propagation rule of [Root Expressions — How Root Expressions Are Propagated](../theory/Root%20Expressions%20and%20LC%20Estimation.md#how-root-expressions-are-propagated) gives the bits of block $j$, computed through the true ANFs of the chaining functions; step 1 of [The Recursion](#the-recursion) shows that it is the same for every bit of the block. This section derives a recursion for $\text{RE}_j$ and the weight-vector description the iterator works with. Neither path evaluates the true ANFs: both compute the same recursion with an exponent $D_i$ read from a monomial profile, and $D_i \ge d_i$, so what they compute dominates $\text{RE}_j$ ([Degree Extraction](#degree-extraction)). The section is written for root expressions; the last paragraph of [Terms and Dominance](#terms-and-dominance) states the changes for monomial profiles.

### The Model Hypothesis

The derivations below assume:

> **(H)** Every bit of block $i \ge 1$ updates as block $i$'s MPR linear update XOR a Boolean function of block $i - 1$'s bits only, and every bit of the source updates by its MPR update alone.

The code enforces part of (H). `root_expressions()` and `monomial_profiles()`, and therefore `estimate_LC`, first call `_require_chaining_linear_in_own_bits()`, on both paths. It raises `ValueError` when a block's merged update profile ([Degree Extraction](#degree-extraction)) contains a product with at least one factor from the block itself, that is, a term with more than one factor containing `VAR(block_id)`; the source block is checked too. The eligibility check for the mesh then requires that each block's update functions read no bits other than the block's own and those of the preceding block (`used_bits <= allowed_bits`); the source may read only its own bits. Two kinds of functions that read a block's own bits beyond its MPR update are outside (H), and the code treats them differently.

- **Nonlinear in the block's own bits: rejected.** On such chaining neither path returns a bound, which is the reason for the rejection. On a register of sizes $(7, 4, 3)$ whose middle block's chaining has terms `AND(own bit, source bit)`, the upper values for a bit of block 1 and a bit of block 2 were 32 and 35 on the mesh path and 4 and 7 on the default path, while Berlekamp–Massey measured linear complexities 90 and 81 (recorded in the docstring of `test_chaining_nonlinear_in_own_bits_is_rejected`). The check reads the as-written profile, so it also rejects a product that does not survive in the true ANF, such as an `AND` of an own bit with itself. It also misses no register whose true ANF has such a monomial. Let $\mu$ be a monomial of a bit's update in block $i$, with at least 2 variables, at least one of them in block $i$. By the [profile coverage](#degree-extraction) claim proved under Degree Extraction, the bit's profile has a term with the same block set as $\mu$ and count at least $|\mu \cap b|$ at each block $b$, and the merged XOR keeps that term or one covering it. That term therefore has at least 2 factors, one of them `VAR(i)`, so the check raises. Before the rejection was added, the two paths failed differently. The default path evaluates the own-bit factor against the block's own table entry, which is still empty while the block's expression is being computed (`block_table[block_id]` receives the result of the same `eval_ANF` call in `_re_default` and `_mp_default`), so the term contributes nothing. The mesh counts the own-bit factors in its rate ([Degree Extraction](#degree-extraction)). Neither has a term combining block $i$'s field with an upstream field: the default because the term vanishes, the mesh because block $i$'s own unit stands alone (Lemma 1(c) below).
- **Linear in the block's own bits: not detected.** An extra term linear in the block's own bits changes the block's update matrix $U$, and with it the characteristic polynomial $\chi_U$ whose roots are the block's own roots in the propagation bound. The recursion's own term $\langle s_i \cdot 1 \rangle$ assumes the MPR's polynomial. This kind is also outside (H), but `_require_chaining_linear_in_own_bits()` does not detect it, since a single own-bit variable is a term with one factor, like the MPR update's own terms.

> **OPEN:** Should chaining that reads a block's own bits linearly be supported, and if so, what own term should the recursion use for it? — needs author input.
> (H) excludes it because it changes $\chi_U$; neither the rejection nor the eligibility check detects it. What the two paths return on such a register has not been measured.

### Terms and Dominance

Two terms from [Root Expressions and LC Estimation](../theory/Root%20Expressions%20and%20LC%20Estimation.md) are used throughout. The **coset weight** of a root $\alpha^k$ of $\mathbb{F}_{2^s}$ is the Hamming weight of the $s$-bit binary representation of its exponent $k$; squaring the root cyclically shifts those bits, so every exponent in the same cyclotomic coset has the same weight. The **coset class** $\langle s \cdot w \rangle$ is the set of all roots in $\mathbb{F}_{2^s}$ whose coset weight is at most $w$.

A **term** of a root expression is a product of coset classes from different fields, $\prod_j \langle s_j \cdot w_j \rangle$, and the expression is the union of its terms. With pairwise coprime sizes, each block has its own field $\mathbb{F}_{2^{s_j}}$, any two of these fields meet only in $\mathbb{F}_2$, and [Proposition 3](../theory/Root%20Expressions%20and%20LC%20Estimation.md#the-five-root-expression-propositions) (cross-field product) makes the root count of a term the product of its factors' counts. A term is therefore described by a **weight vector** $k = (k_0, \ldots, k_m)$ of nonnegative integers: block $j$ contributes the factor $\langle s_j \cdot \min(k_j, s_j) \rangle$ when $k_j > 0$, and no factor when $k_j = 0$. Write

$$\bar k_j = \min(k_j, s_j), \qquad \mathrm{supp}(k) = \{\, j : k_j > 0 \,\}$$

for the **saturated weight** at block $j$ and the set of blocks the term uses. Weight above $s_j$ changes nothing: coset weights in $\mathbb{F}_{2^{s_j}}$ run only up to $s_j$, and Proposition 2 caps a product's weight there. So two vectors with the same support and the same saturated weights describe the same term.

Multiplying two terms adds their weight vectors. A field used by only one factor is carried over, and in a field used by both, Proposition 2 (same-field product) adds the coset weights, capped at $s_j$. Adding first and capping later gives the same result, because for nonnegative integers $y_1, \ldots, y_n$

$$\min\Bigl(s,\ \sum_r y_r\Bigr) = \min\Bigl(s,\ \sum_r \min(s, y_r)\Bigr):$$

when every $y_r < s$ the inner minima are the $y_r$ themselves and the two sides coincide, and when some $y_r \ge s$ both sides equal $s$. So the term of $k + k'$ is the product of the terms of $k$ and $k'$, and the iterator may carry weights above $s_j$, capping them only when it outputs a term.

A term $k$ is **dominated** by a term $k'$ when the two use the same blocks and $\bar k_j \le \bar k'_j$ for every $j$. Since $\langle s \cdot w \rangle$ grows with $w$, every root of $k$'s term is then a root of $k'$'s, and $k$ adds nothing to the bound. It is **strictly dominated** when, in addition, $\bar k_j < \bar k'_j$ for some $j$. A term of an expression is **maximal** when no term of the expression strictly dominates it. Dominance compares only terms over the same fields, and it is the comparison both paths prune with: `is_root_subset` is applied within one field tuple in `RootExpression.__add__` and `RootExpression.__mul__`, and `is_monomial_subset` requires equal block sets. It is weaker than containment of root sets. The class $\langle 7 \cdot 7 \rangle$ contains the root $1$, so $\langle 7 \cdot 7 \rangle \times \langle 5 \cdot 3 \rangle$ contains every root of $\langle 5 \cdot 3 \rangle$, yet neither term dominates the other. "Maximal" below always means maximal for dominance.

Two facts about dominance are used below.

- **(a) Dominance survives multiplication.** If $k$ is dominated by $k'$, then $k + x$ is dominated by $k' + x$ for every $x$. The two sums use the same blocks, $\mathrm{supp}(k) \cup \mathrm{supp}(x)$, and by the identity above the saturated weight of a sum at block $j$ is $\min(s_j, \bar k_j + \bar x_j)$, which does not decrease when $\bar k_j$ is replaced by the larger $\bar k'_j$.
- **(b) Pruning before multiplying does not change the maximal terms.** By (a), a product that uses a dominated term is dominated by the product that uses the dominating term in its place. So removing dominated terms from an expression before multiplying it leaves the maximal terms of the product unchanged. This is what allows the default path, which prunes after every operation, to be compared with the mesh, which builds no intermediate expressions.

**Monomial profiles.** The definitions and facts above, Lemma 1 below, and [Correctness of the Iterator](#correctness-of-the-iterator) use only this arithmetic of weight vectors, and they hold for monomial profiles with three changes. The pair $\langle s_j \cdot w \rangle$ denotes the monomials with 1 to $w$ variables from the $s_j$-bit register ([Monomial Profile Theory](../theory/Monomial%20Profile%20Theory.md)). Coordinates are registers rather than fields, so two blocks of the same size are still two coordinates. And Propositions 2 and 3 of that document replace those of root expressions: weights from one register add, capped at the register's $s_j$ variables (`TermSet.__mul__` adds counts per block and caps them at the block size), and counts from different registers multiply whatever their sizes.

### The Recursion

Assume (H). For a block $i \ge 1$ and a bit $b$ of it, write $\mathcal{C}_{i,b}$ for the chaining function of bit $b$: the function of block $i-1$'s bits that (H) adds to the bit's MPR update. Its **algebraic degree** is the largest number of variables in a monomial of its [ANF](../theory/Algebraic%20Normal%20Form.md). The **chaining degree** $d_i$ is the largest of these degrees over the bits $b$ of block $i$. The code reads an exponent $D_i \ge d_i$ off a monomial profile instead; [Degree Extraction](#degree-extraction) says how, and proves the inequality.

Suppose for now that $d_i \ge 1$ for every $i \ge 1$ and that no chaining function has a constant term; [Constants and Blocks Without Chaining Input](#constants-and-blocks-without-chaining-input) treats both exceptions. Then

$$\text{RE}_0 = \langle s_0 \cdot 1 \rangle, \qquad \text{RE}_i = (\text{RE}_{i-1})^{d_i} + \langle s_i \cdot 1 \rangle \quad (1 \le i \le m),$$

where $+$ is the union of terms and the power is the product of $d_i$ copies. The power is expanded by distributing: a term of $(\text{RE}_{i-1})^{d_i}$ is the product of one term of $\text{RE}_{i-1}$ chosen for each of the $d_i$ factors, which is how `RootExpression.__mul__` forms a product (every pair of terms). The equation holds up to dominated terms: its two sides have the same maximal terms. The derivation has three steps.

1. **Each bit of block $i$ is bounded by its own roots and the roots of the block's whole chaining input.** [Root Expressions — Why the Resolvent Justifies Propagation](../theory/Root%20Expressions%20and%20LC%20Estimation.md#why-the-resolvent-justifies-propagation) shows that every entry of the block's output is a linear combination of the entries of its chaining input $\mathcal{C}$, filtered through rational functions with denominator $\chi_U$, so the root expression of every bit of block $i$ is bounded by $\text{RE}(\mathcal{C}) + \langle s_i \cdot 1 \rangle$. A linear combination of sequences has its roots among theirs (Proposition 4), so $\text{RE}(\mathcal{C})$ is the union of $\text{RE}(\mathcal{C}_{i,b})$ over the bits $b$ of the block, and the bound is the same for every bit.
2. **A monomial of degree $a$ contributes $(\text{RE}_{i-1})^a$.** Each $\text{RE}(\mathcal{C}_{i,b})$ is computed from the ANF ([Propagation Through Boolean Functions](../theory/Root%20Expressions%20and%20LC%20Estimation.md#propagation-through-boolean-functions)): AND multiplies root expressions and XOR takes their union. A monomial of degree $a$ is an AND of $a$ bits of block $i-1$, each with root expression $\text{RE}_{i-1}$. So $\text{RE}(\mathcal{C})$ is the union of $(\text{RE}_{i-1})^a$ over the monomials of all the block's chaining functions, and the largest exponent $a$ occurring is $d_i$.
3. **Lower powers are dominated.** Let $E$ be a root expression and $1 \le a < d$. A term of $E^a$ is a sum $x_1 + \cdots + x_a$ of term vectors of $E$. Adding $d - a$ further copies of $x_1$ gives a term of $E^d$ that uses the same blocks (those of $x_1$ are already among them) and whose saturated weights are at least as large. So every term of $E^a$ is dominated by a term of $E^d$, and the union in step 2 has the same maximal terms as $(\text{RE}_{i-1})^{d_i}$ alone.

For monomial profiles, the default path builds the same recursion directly: `_mp_default` evaluates block $i$'s chaining profile on the table of profiles already computed and adds the term with count 1 at block $i$ (`TermSet({block_id: size}, {block_id: 1})`), the profile counterpart of $\langle s_i \cdot 1 \rangle$.

When `locked_list` locks blocks, the recursion drops the own term $\langle s_j \cdot 1 \rangle$ of each locked block; [Locked Blocks](#locked-blocks) states that recursion and proves what the mesh computes for it.

### Constants and Blocks Without Chaining Input

The constant sequence 1 has the single root $1$ ([Propagation Through Boolean Functions](../theory/Root%20Expressions%20and%20LC%20Estimation.md#propagation-through-boolean-functions)). In a root expression it is the **constant term**, the term with no fields (`RootExpression.logical_one()`), and its weight vector is $0$. The constant term is not the same as the root $1$: a full class $\langle s \cdot s \rangle$ also contains the root $1$, but it is a term over the field $\mathbb{F}_{2^s}$. The constant term enters $\text{RE}_i$ in two ways.

- **A constant term in block $i$'s chaining.** If some $\mathcal{C}_{i,b}$ has the constant 1 in its ANF, step 2 puts the constant term into $\text{RE}(\mathcal{C})$, so $\text{RE}_i$ contains it.
- **A constant term inherited from upstream.** Suppose $\text{RE}_{i-1} = E + 1$, where $E$ does not contain the constant term. Choosing $1$ for $d_i - a$ of the factors and a term of $E$ for the other $a$ shows that $(E + 1)^{d_i}$ is the union of $E^a$ for $0 \le a \le d_i$, where $E^0 = 1$. By step 3, the powers $E^a$ with $1 \le a < d_i$ are dominated by $E^{d_i}$. So the maximal terms of $(\text{RE}_{i-1})^{d_i}$ are those of $E^{d_i}$ together with the constant term.

In both cases the constant term stands by itself: no other term has empty support, so it neither dominates nor is dominated by any other. The remaining terms are those of the recursion without constants. The constant term can therefore be handled separately. The mesh computes the recursion without constants, and the driver adds the constant term to $\text{RE}_i$ when block $i$'s chaining has a constant term or $\text{RE}_{i-1}$ contains the constant term ([Output Construction](#output-construction)).

**A block with no chaining input.** If no $\mathcal{C}_{i,b}$ reads a bit of block $i-1$, each of them is a constant, 0 or 1, and $d_i = 0$. The power $(\text{RE}_{i-1})^0$ is then the empty product, and it stands for that constant: it contributes nothing when every chaining function is 0, and the constant term when some chaining function is the constant 1. So $\text{RE}_i = \langle s_i \cdot 1 \rangle$, together with the constant term in the second case. The code does not run the mesh with $d_i = 0$. A register in which some block after the source takes no input from upstream is handled entirely by the default algorithm ([Degree Extraction](#degree-extraction)).

### Weight Vectors and Their Value

In weight-vector form, $\langle s_j \cdot 1 \rangle$ is the unit vector $e_j$ (weight 1 at block $j$ and 0 elsewhere), and a product of terms is a sum of vectors. Let $T_j$ be the set of term vectors of $\text{RE}_j$ expanded without removing any term, constants excluded. The recursion says

$$T_0 = \{ e_0 \}, \qquad T_j = \{ e_j \} \cup \{\, x_1 + \cdots + x_{d_j} : x_1, \ldots, x_{d_j} \in T_{j-1} \,\} \quad (1 \le j \le m).$$

By fact (b) of [Terms and Dominance](#terms-and-dominance), $T_m$ has the same maximal terms as $\text{RE}_m$ computed with pruning at every step.

Set $P_0 = 1$ and $P_j = d_1 d_2 \cdots d_j$, and define the **value** of a vector as

$$V(k) = \sum_{j=0}^{m} k_j P_j .$$

$P_j$ is the worth of one unit of weight at block $j$, measured in units at block 0: one unit at block $j$ is worth $d_j$ units at block $j - 1$, and so $d_j d_{j-1} \cdots d_1$ units at block 0. A **trade** at block $j$, for $1 \le j \le m$, exchanges one unit at block $j$ for $d_j$ units at block $j - 1$:

$$k \;\longmapsto\; k - e_j + d_j\, e_{j-1}, \qquad \text{allowed when } k_j \ge 1.$$

A trade preserves value, since $-P_j + d_j P_{j-1} = 0$. In terms of the recursion, a unit at block $j$ stands for choosing the summand $\langle s_j \cdot 1 \rangle$ of $\text{RE}_j$. Trading it chooses the other summand, $(\text{RE}_{j-1})^{d_j}$, and the $d_j$ units it produces at block $j - 1$ stand for that power's $d_j$ factors, each of which is in turn kept or traded.

**Lemma 1.** Let $R_j$ be the set of vectors reachable from $e_j$ by trades at blocks $1, \ldots, j$. Then:

- (a) $R_j = T_j$;
- (b) every vector of $R_j$ has value $P_j$;
- (c) $e_j$ is the only vector of $R_j$ with $k_j \ne 0$.

*Proof.* (b) $V(e_j) = P_j$, and trades preserve value.

(c) Entries above $j$ are 0 in every vector of $R_j$, since a trade at a block $i \le j$ changes only blocks $i$ and $i - 1$. So a vector $k \in R_j$ with $k_j \ge 1$ has $P_j = V(k) = k_j P_j + \sum_{i < j} k_i P_i$. Every $P_i$ is positive because every $d_i \ge 1$, so $k_j = 1$ and $k_i = 0$ for $i < j$.

(a) By induction on $j$. For $j = 0$ no trade is available, since trades happen at blocks $1$ and above, so $R_0 = \{ e_0 \} = T_0$. Let $j \ge 1$. From $e_j$ the only allowed trade is at block $j$, and it gives $d_j e_{j-1}$. After it, block $j$ is empty, and no trade refills it, since every allowed trade adds weight at a block below $j$. So $R_j = \{ e_j \} \cup Q$, where $Q$ is the set reachable from $d_j e_{j-1}$ by trades at blocks $1, \ldots, j - 1$. It remains to show that $Q$ is the set of sums of $d_j$ vectors of $R_{j-1}$; the induction hypothesis $R_{j-1} = T_{j-1}$ then gives $R_j = T_j$.

- *Every such sum lies in $Q$.* Write $d_j e_{j-1}$ as the sum of $d_j$ copies of $e_{j-1}$, apply the trades that turn the first copy into the first summand, then those for the second, and so on. Each trade is allowed: the entry it draws from is at least the corresponding entry of the summand being built, because the other summands are nonnegative, and that entry is at least 1.
- *Every vector of $Q$ is such a sum.* By induction on the number of trades. The vector $d_j e_{j-1}$ is the sum of $d_j$ copies of $e_{j-1} \in R_{j-1}$. Suppose $y = x_1 + \cdots + x_{d_j}$ with every $x_r \in R_{j-1}$, and a trade at block $i$ is allowed on $y$. Then $y_i \ge 1$, so some $x_r$ has $(x_r)_i \ge 1$. The trade is allowed on that $x_r$ and keeps it in $R_{j-1}$, and the traded $y$ is the sum with $x_r$ replaced by its traded form. $\square$

Read in terms of the recursion, the inductive step says two things. The power multiplies the value by the new exchange rate: a sum of $d_j$ vectors of value $P_{j-1}$ has value $d_j P_{j-1} = P_j$. And the new summand $\langle s_j \cdot 1 \rangle$ is the single unit $e_j$, whose value is the whole budget $P_j$, so it never appears together with any other weight.

### Worked Example

Take a CMPR with blocks of sizes $[7, 5, 3]$ and chaining degrees $d_1 = d_2 = 3$. The code's degree list is $[1, 3, 3]$; its first entry is the degree extracted for the source block, which is never used as a rate ([Tree Traversal](#tree-traversal)). The recursion gives

$$\text{RE}_1 = \langle 7 \cdot 1 \rangle^3 + \langle 5 \cdot 1 \rangle = \langle 7 \cdot 3 \rangle + \langle 5 \cdot 1 \rangle, \qquad \text{RE}_2 = \bigl( \langle 7 \cdot 3 \rangle + \langle 5 \cdot 1 \rangle \bigr)^3 + \langle 3 \cdot 1 \rangle .$$

Choosing $\langle 7 \cdot 3 \rangle$ for $a$ of the three factors and $\langle 5 \cdot 1 \rangle$ for the other $3 - a$ gives the term with vector $(3a,\ 3 - a,\ 0)$. This is the binomial expansion of the cube with its coefficients dropped: a root expression records which products can occur, not how many times. The mesh for block 2 enumerates:

| Weight vector | Root expression term |
|---|---|
| $[0,\; 0,\; 1]$ | $\langle 3 \cdot 1 \rangle$ |
| $[0,\; 3,\; 0]$ | $\langle 5 \cdot 3 \rangle$ |
| $[3,\; 2,\; 0]$ | $\langle 7 \cdot 3 \rangle \times \langle 5 \cdot 2 \rangle$ |
| $[6,\; 1,\; 0]$ | $\langle 7 \cdot 6 \rangle \times \langle 5 \cdot 1 \rangle$ |
| $[9,\; 0,\; 0]$ | $\langle 7 \cdot 7 \rangle$ |

The values are $P = (1, 3, 9)$, and every row satisfies $k_0 + 3k_1 + 9k_2 = 9$. The first row is $e_2$, the driven block's own $\langle 3 \cdot 1 \rangle$. The second is the trade at block 2, the case $a = 0$. The remaining rows trade one, two, and three of the units at block 1 to block 0, each for three units there: the cases $a = 1, 2, 3$. Weight in $\mathbb{F}_{2^7}$ saturates at 7, so $[9, 0, 0]$ maps to $\langle 7 \cdot 7 \rangle$, not $\langle 7 \cdot 9 \rangle$. The only two rows over the same fields, $[3, 2, 0]$ and $[6, 1, 0]$, have saturated weights $(3, 2)$ and $(6, 1)$, and neither dominates the other, so all five terms are maximal.

The weights in a row do not sum to at most $d_2 = 3$: the row $[9, 0, 0]$ has total weight 9. What every row has in common is its value.

## Why Root Expressions Need Pairwise Coprime Sizes

When two blocks' sizes share a factor, their fields share more than $\mathbb{F}_2$. In this case:

- **Cross-field independence breaks down.** Proposition 3 requires pairwise coprime field degrees, and sizes that merely differ but share a factor (such as 2 and 4) already break it, because the two fields share a subfield. Two blocks of the same size contribute roots to the same axis in the [coset region grid](../conventions/Coset%20Regions.md), so their products are same-field products (Proposition 2) that add coset weights, not independent axes whose counts multiply.
- **The weight vector no longer has one coordinate per field.** The mesh gives each block its own coordinate and treats different coordinates as different fields, whose counts multiply. With shared sizes, two coordinates denote the same field, and the correct result requires tracking how same-field products interact through Jordan partitions (see [Root Multiplicities — Jordan Decomposition and Products](../theory/Root%20Multiplicities%20and%20Jordan%20Decomposition.md#jordan-decomposition-and-products)).

Monomial profiles need no condition on sizes. Their coordinates are registers rather than fields, and Proposition 3 of [Monomial Profile Theory](../theory/Monomial%20Profile%20Theory.md) multiplies counts across different registers whatever their sizes.

T-functions are excluded from the root-expression mesh by a different condition. Every block of a T-function uses the same polynomial $(1 + x)$ and has size 1, and repeated blocks of size 1 do not violate coprimality, since $\gcd(1, 1) = 1$; it is the one-bit condition of [Eligibility Conditions](#eligibility-conditions) that excludes them. `TFunction` overrides both `root_expressions` and `monomial_profiles` in any case ([T-Functions](Feedback%20Family%20Implementation%20Details/TFunctions.md)).

## The Iterator

The mesh iterator (`_re_mesh_iterator` in `MeshOptimization.py`) enumerates the vectors of $R_m$ for the driven block $m$ by a depth-first traversal of trades, skips some of them, and yields those that pass two further checks. The iterator is numba-jitted for performance. Its inputs are the block sizes, the degree list, and the lock flags. [Correctness of the Iterator](#correctness-of-the-iterator) proves that, when no block is locked, it yields each maximal term of $T_m$ exactly once, and states what it computes when blocks are locked.

### Degree Extraction

Before the mesh runs, the CMPR driver extracts a degree for each block. The pipeline is:

1. `MonomialProfile.from_merged()` computes the monomial profile of the XOR of the block's update functions (each bit's MPR update together with its chaining), recording each variable by the index of its block. The XOR is evaluated on profiles, where it is a union, so a monomial in one bit's update is not cancelled by the same monomial in another's; the profile covers every bit of the block.
2. `.to_BooleanFunction()` converts to an ANF representation where each ANF term becomes an `AND` of `VAR(block_id)` references, with repeated VARs representing same-field powers. Because of the repeats, the ordinary degree of the function cannot be used: ANF degree counts distinct variables, so a repeated VAR would collapse by $x \cdot x = x$.
3. The degree is `max(len(term.args))` across all ANF terms — the maximum number of factors in any single product.

The profile bounds the block's functions as they are written, not their true ANFs. An `AND` of the same variable twice counts 2, although $x \cdot x = x$, and a monomial occurring twice in one XOR is kept, although the two copies cancel. It never misses a monomial of the true ANF, which is what the next claim states.

**Profile coverage.** Let $f$ be a node of a function's DAG, and $\mu$ a monomial of $f$'s true ANF. Then the profile $P(f)$ that `eval_ANF` computes contains a term $t$ whose block set is the set of blocks $\mu$ meets, and whose count at each such block $b$ is at least $|\mu \cap b|$, the number of $\mu$'s variables in block $b$.

*Proof.* By induction over the nodes in the postorder `eval_ANF` uses.

- A `VAR` of a bit in block $b$ gets the profile `TermSet({b: s_b}, {b: 1})` (`from_merged`), which covers its one monomial. `CONST(1)` gets the empty term (`MonomialProfile.logical_one()`), which covers the empty monomial; `CONST(0)` has no monomials.
- **XOR.** A monomial appears in $\mathrm{ANF}(f \oplus g)$ only if it appears in exactly one of $\mathrm{ANF}(f)$ and $\mathrm{ANF}(g)$, so it is covered by a term of $P(f)$ or of $P(g)$. `MonomialProfile.__add__` keeps the maximal terms of both (`maximal_elements` with `is_monomial_subset`), and removes a term only when another term with the same blocks and counts at least as large remains, which covers whatever the removed one covered.
- **AND.** Since $x \cdot x = x$, every monomial of $\mathrm{ANF}(fg)$ is $\mu_f \cup \mu_g$ for monomials $\mu_f$ of $f$ and $\mu_g$ of $g$. If $t_f$ and $t_g$ cover them, `TermSet.__mul__` gives $t_f t_g$ the union of their blocks and, at block $b$, the count $\min(s_b, c_f + c_g) \ge \min(s_b, |\mu_f \cap b| + |\mu_g \cap b|) \ge |(\mu_f \cup \mu_g) \cap b|$; the last step holds because a monomial has at most $s_b$ variables in block $b$. The pruning in `MonomialProfile.__mul__` is the same as for XOR.
- Every other gate is evaluated through these: the `_eval_ANF` methods in `Gates.py` combine `&`, `^`, and `invert`, and `invert` on a profile is `MonomialProfile.__invert__`, an XOR with `logical_one()`. $\square$

Write $D_i$ for the largest count at block $i - 1$ among the profile terms of block $i$ whose only block is $i - 1$. Under (H), take a bit $b$ of block $i$ and a monomial of degree $d_i$ in $\mathcal{C}_{i,b}$. It lies in the true ANF of that bit's update, because the MPR update's monomials are single variables of block $i$. By the claim it is covered by a term of the bit's profile with only block $i - 1$ and count at least $d_i$, and the merged XOR keeps that term or one covering it. So $D_i \ge d_i$. The cap in `TermSet.__mul__` (`MonomialProfile.py`, lines 40–41) gives $D_i \le s_{i-1}$. Equality $D_i = d_i$ holds when each chaining function is written as an XOR of distinct monomials, each an AND of distinct variables. The default path evaluates the same profile, so the two paths use the same exponent $D_i$; and by step 3 of [The Recursion](#the-recursion), every term of $(\text{RE}_{i-1})^{d_i}$ is dominated by a term of $(\text{RE}_{i-1})^{D_i}$, so the larger exponent keeps the result an upper bound. The rates the iterator receives are these $D_i$. The recursion, Lemma 1, and [Correctness of the Iterator](#correctness-of-the-iterator) assume only rates of at least 1, so they apply with $D_i$ in place of $d_i$. The same claim at the empty monomial shows that a constant term of the true ANF appears as `CONST(1)` in the profile, so `constants_possible`, which is read from the profile, over-approximates the constant terms.

The extracted degree also counts factors from the block's own bits. The block's MPR update is linear, so it contributes only terms with a single variable of the block itself, each counting 1. Under (H), when the block chains from upstream ($D_i \ge 1$), the extracted degree is therefore $\max(1, D_i) = D_i$. A block after the source with no upstream input has $d_i = D_i = 0$, but its extracted degree is 1, from its own linear update. Run with that degree, the mesh would trade through the block as if it chained linearly. The drivers detect such a block (`chains_upstream` in `_re_mesh_optimization` and `_mp_mesh_optimization`) and compute the whole register with the default algorithm; `test_mesh_path_on_a_block_without_chaining_input` in `test/Tools/linear_complexity_test.py` covers this case. A term multiplying own bits would also enter `len(term.args)`, but such terms are rejected before any degree is extracted ([The Model Hypothesis](#the-model-hypothesis)). The degree extracted for the source block itself is never used as a rate.

### Tree Traversal

The traversal runs from the driven block, the rightmost array position, toward the source, the leftmost; array position $j$ is block $j$.

1. **Initialization:** `arr[0][-1] = 1` places weight 1 at block $m$, giving the vector $e_m$, which is yielded at once. It is the summand $\langle s_m \cdot 1 \rangle$: the driven block's own roots, which the resolvent contributes whatever the chaining ([Root Expressions — Why the Resolvent Justifies Propagation](../theory/Root%20Expressions%20and%20LC%20Estimation.md#why-the-resolvent-justifies-propagation)). By Lemma 1(c) it stands alone.

2. **Trading:** The row `arr[depth]` handles block $j = m - \texttt{depth}$. Each pass trades `num_traded` units from block $j$ to block $j - 1$ (`arr[depth][block] -= num_traded` and `arr[depth][block-1] += num_traded * degrees[block]`), which is `num_traded` trades of [Weight Vectors and Their Value](#weight-vectors-and-their-value) at rate $d_j$ and preserves the value. After each pass the iterator considers the new vector for yielding, copies the row to the next depth, and descends to block $j - 1$. It returns to the previous depth when block $j$ is empty or $j = 0$; nothing is traded out of the source. The count is `num_traded = max(1, weight - size)`: one unit per pass, except that the first pass out of a block holding more than $s_j + 1$ units moves all but $s_j$ of them at once (the jump rule under [Pruning](#pruning)).

3. **Saturation:** A block with weight $\geq s_j$ is saturated: additional weight does not enlarge its coset class (Proposition 2 caps the weight at $s_j$). The excess beyond $s_j$ is what the jump rule trades on and what the fillable space rule measures.

4. **Accumulated rates:** The `values` array stores $\texttt{values}[i] = d_0 \cdot d_1 \cdots d_i = d_0 P_i$, where $d_0 = $ `degrees[0]` is the degree extracted for the source block ($d_0 \ge 1$, from the source's own linear update). The iterator uses only ratios $\texttt{values}[j] / \texttt{values}[p] = P_j / P_p$, in which $d_0$ cancels. For $j > p$ this ratio is $d_{p+1} \cdots d_j$, the number of units at block $p$ that one unit at block $j$ is worth.

### Pruning

Three rules keep the iterator from yielding dominated or repeated terms. [Correctness of the Iterator](#correctness-of-the-iterator) proves what each one does; this subsection states them.

**Jump rule.** This is the `num_traded` count of step 2. When block $j \ge 1$ holds more than $s_j + 1$ units before its first trade, that trade moves all but $s_j$ of them at once. The vectors in which block $j$ has traded some units yet still holds more than $s_j$ are therefore never visited. The visited vectors are exactly those of $R_m$ in which no block holds more than its size, except possibly the leftmost nonzero block (Lemma 3 below). Every skipped vector is dominated by a visited one, obtained by moving the excess at the overfull block on to the leftmost nonzero block (Lemma 4).

**Leftmost saturation rule.** The iterator yields a vector only if no saturated block other than the leftmost saturated one holds more than its size. Weight above saturation is invisible in the term, so the rule is meant to keep one placement of that excess. Without locked blocks, every vector the traversal visits already satisfies it (Lemma 5), because of the jump rule. With locked blocks, whose sizes are set to 0, it can reject vectors; Lemma 11 says which, and shows that the rejected placement's term is still produced.

**Fillable space rule.** A block $p$ is **partially filled** when $0 < k_p < s_p$. The iterator converts the excess $\max(0, k_j - s_j)$ of every block into units at $p$ at the accumulated rates,

$$\frac{1}{\texttt{values}[p]} \sum_{j} \max(0,\, k_j - s_j)\, \texttt{values}[j],$$

and skips the vector if this is at least 1 for some partially filled $p$. On a visited vector, only the leftmost nonzero block $\ell$ can hold excess, and the quantity is $(k_\ell - s_\ell) P_\ell / P_p$. The rule fires when the excess at $\ell$ is worth at least one unit at $p$. Undoing that much trading moves one unit back to $p$ and leaves $\ell$ saturated, giving a vector of $R_m$ that strictly dominates the skipped one (Lemma 6; Lemma 10 for locked blocks).

### Pruning Example

The example of [Worked Example](#worked-example) triggers neither the jump nor the fillable space rule. This one triggers both. Take sizes $[5, 3, 2, 7]$ and chaining degrees $d_1 = 2$, $d_2 = 3$, $d_3 = 2$ (degree list $[1, 2, 3, 2]$). These degrees are possible: a function of the $s_{i-1}$ bits of block $i - 1$ has degree at most $s_{i-1}$, and $2 \le 5$, $3 \le 3$, $2 \le 2$. Here $P = (1, 2, 6, 12)$, so every vector satisfies $k_0 + 2k_1 + 6k_2 + 12k_3 = 12$. The thirteen vectors of $R_3$, in traversal order, with the trade counts $t_j$ (the number of units traded out of block $j$) and what happens to each:

| Vector $k$ | $(t_1, t_2, t_3)$ | Outcome | Term |
|---|---|---|---|
| $[0,\; 0,\; 0,\; 1]$ | $(0, 0, 0)$ | yielded | $\langle 7 \cdot 1 \rangle$ |
| $[0,\; 0,\; 2,\; 0]$ | $(0, 0, 1)$ | yielded | $\langle 2 \cdot 2 \rangle$ |
| $[0,\; 3,\; 1,\; 0]$ | $(0, 1, 1)$ | yielded | $\langle 3 \cdot 3 \rangle \times \langle 2 \cdot 1 \rangle$ |
| $[2,\; 2,\; 1,\; 0]$ | $(1, 1, 1)$ | yielded | $\langle 5 \cdot 2 \rangle \times \langle 3 \cdot 2 \rangle \times \langle 2 \cdot 1 \rangle$ |
| $[4,\; 1,\; 1,\; 0]$ | $(2, 1, 1)$ | yielded | $\langle 5 \cdot 4 \rangle \times \langle 3 \cdot 1 \rangle \times \langle 2 \cdot 1 \rangle$ |
| $[6,\; 0,\; 1,\; 0]$ | $(3, 1, 1)$ | yielded | $\langle 5 \cdot 5 \rangle \times \langle 2 \cdot 1 \rangle$ |
| $[0,\; 6,\; 0,\; 0]$ | $(0, 2, 1)$ | yielded | $\langle 3 \cdot 3 \rangle$ |
| $[2,\; 5,\; 0,\; 0]$ | $(1, 2, 1)$ | not visited (jump rule) | — |
| $[4,\; 4,\; 0,\; 0]$ | $(2, 2, 1)$ | not visited (jump rule) | — |
| $[6,\; 3,\; 0,\; 0]$ | $(3, 2, 1)$ | yielded | $\langle 5 \cdot 5 \rangle \times \langle 3 \cdot 3 \rangle$ |
| $[8,\; 2,\; 0,\; 0]$ | $(4, 2, 1)$ | rejected (fillable space) | — |
| $[10,\; 1,\; 0,\; 0]$ | $(5, 2, 1)$ | rejected (fillable space) | — |
| $[12,\; 0,\; 0,\; 0]$ | $(6, 2, 1)$ | yielded | $\langle 5 \cdot 5 \rangle$ |

When $t_2 = 2$, block 1 receives 6 units, more than $s_1 + 1 = 4$, so its first trade moves 3 at once and $t_1$ jumps from 0 to 3. The two vectors skipped, $[2, 5, 0, 0]$ and $[4, 4, 0, 0]$, are those in which block 1 has traded but still holds more than 3. Moving their excess at block 1 to block 0 (2 units become 4, and 1 unit becomes 2) gives $[6, 3, 0, 0]$ in both cases, which dominates them: saturated weights $(5, 3)$ against $(2, 3)$ and $(4, 3)$. The vector $[8, 2, 0, 0]$ has excess 3 at block 0, worth $3/2 \ge 1$ unit at block 1, which is partially filled ($2 < 3$), so the fillable space rule rejects it; moving one unit back to block 1 gives $[6, 3, 0, 0]$, which strictly dominates it. The vector $[10, 1, 0, 0]$ has excess 5 at block 0, worth $5/2$ units at block 1; the vector Lemma 6 gives for it is $[8, 2, 0, 0]$, which strictly dominates it (saturated weights $(5, 2)$ against $(5, 1)$). The vector $[6, 0, 1, 0]$ has excess 1 at block 0, but that is worth only $1/6$ of a unit at the partially filled block 2, so it is kept; no other vector of $R_3$ uses exactly blocks 0 and 2. The two vectors over blocks $\{0, 1, 2\}$, $[2, 2, 1, 0]$ and $[4, 1, 1, 0]$, are incomparable. The nine yielded terms are the maximal terms of $T_3$.

With three blocks the jump rule cannot fire: the weight arriving at block 1 is the rate $D_2$, and $D_2 \le s_1$ because `TermSet.__mul__` caps counts at the block size ([Degree Extraction](#degree-extraction)), so it is never more than $s_1 + 1$.

## Correctness of the Iterator

This section proves that, when no block is locked, the vectors `_re_mesh_iterator` yields for the driven block $m$ represent exactly the maximal terms of $T_m$, each once (Theorem 1). By fact (b) of [Terms and Dominance](#terms-and-dominance) these are the maximal terms of $\text{RE}_m$ other than the constant term, which the driver adds separately ([Constants and Blocks Without Chaining Input](#constants-and-blocks-without-chaining-input)). [Locked Blocks](#locked-blocks) then treats locking (Theorem 2).

Throughout, $m \ge 1$ and every rate satisfies $d_j \ge 1$ for $1 \le j \le m$; nothing else about the rates is used, so the results apply to the rates $D_j$ the iterator receives. Until [Locked Blocks](#locked-blocks), every size satisfies $s_j \ge 1$. The sizes cover root expressions (sizes at least 2) and monomial profiles (sizes at least 1, repeats allowed). The rates are what the drivers guarantee by falling back when a block has no upstream input.

The induction over blocks is Lemma 1: passing from $\text{RE}_{j-1}$ to $\text{RE}_j$ multiplies the value by $d_j$ and adds the single unit $e_j$. Each call of the iterator enumerates $R_m$ afresh rather than reusing the previous block's output, so the remaining lemmas analyze a single call, for arbitrary $m$.

**Lemma 2** (trade counts). For $k \in R_m$, let $t_j$ be the number of trades at block $j$ ($1 \le j \le m$) in some sequence of trades reaching $k$ from $e_m$, and set $t_0 = 0$. Then

$$k_m = 1 - t_m, \qquad k_j = d_{j+1}\, t_{j+1} - t_j \quad (0 \le j < m).$$

(i) The counts $t$ are determined by $k$, so they do not depend on the sequence. (ii) Conversely, every integer vector $t$ with $t_0 = 0$, $t_m \in \{0, 1\}$, and $0 \le t_j \le d_{j+1} t_{j+1}$ for $1 \le j < m$ is the trade-count vector of a vector of $R_m$. (iii) For $1 \le j \le m$,

$$t_j P_j = \sum_{i < j} k_i P_i ,$$

so $t_j = 0$ exactly when $k_i = 0$ for every $i < j$.

*Proof.* Each trade at block $j$ removes one unit from block $j$ and adds $d_j$ units to block $j - 1$. Starting from $e_m$, block $j < m$ therefore receives $d_{j+1} t_{j+1}$ units and loses $t_j$, and block $m$ starts with 1 and loses $t_m$; this gives the formulas. (i) Solving them from $j = m$ downward, $t_m = 1 - k_m$ and $t_j = d_{j+1} t_{j+1} - k_j$. (ii) The constraint $t_j \le d_{j+1} t_{j+1}$ is $k_j \ge 0$. Perform $t_m$ trades at block $m$, then $t_{m-1}$ at block $m-1$, and so on down to block 1. When the trades at block $j$ begin, block $j$ holds $d_{j+1} t_{j+1} \ge t_j$ units (1 unit if $j = m$), so each of them is allowed. (iii) Substituting the formulas, and using $t_0 = 0$ and $d_{i+1} P_i = P_{i+1}$,

$$\sum_{i<j} k_i P_i = \sum_{i=0}^{j-1} \bigl( t_{i+1} P_{i+1} - t_i P_i \bigr) = t_j P_j .$$

Every $P_i$ is positive and every $k_i$ is nonnegative, so the sum is 0 exactly when every $k_i$ with $i < j$ is 0. $\square$

Call a block $j \ge 1$ of $k \in R_m$ **overfull** if $t_j \ge 1$ and $k_j > s_j$: it has traded weight away, yet still holds more than its size. By Lemma 2(iii), a block $j$ with $k_j > s_j$ is not overfull exactly when every block to its left is empty, that is, when $j$ is the leftmost nonzero block. Let

$$W = \{\, k \in R_m : k \text{ has no overfull block} \,\},$$

the vectors of $R_m$ in which no block holds more than its size except possibly the leftmost nonzero one.

**Lemma 3** (visited vectors). The traversal visits each vector of $W$ exactly once and visits no other vector.

*Proof.* The row at depth $m - j$ holds a vector whose trade counts are fixed at blocks above $j$ and are 0 at blocks below $j$. Block $j$ starts with everything it received, $I = d_{j+1} t_{j+1}$ units ($I = 1$ when $j = m$), and $t_j = 0$. Each pass trades $\max(1, k_j - s_j)$ units. If $I > s_j$, the first pass trades $I - s_j$, leaving $s_j$, and every later pass trades 1. If $I \le s_j$, every pass trades 1. Passes stop when block $j$ is empty. So $t_j$ takes the value 0 and every value $t$ with $\max(1, I - s_j) \le t \le I$. These are exactly the values with $t_j = 0$ or $k_j = I - t_j \le s_j$, the values for which block $j$ is not overfull.

After each pass with $t_j \ge 1$, the row is copied one level down, where block $j - 1$ is handled in the same way; at block 0 nothing is traded. A vector with $t_j = 0$ has $t_i = 0$ for every $i < j$, because $t_i \le d_{i+1} t_{i+1}$, so it is the vector the level above already holds and is not produced again. Hence each vector of $W$ is produced at exactly one moment: by the trade at the lowest block with a nonzero trade count, or, for $e_m$, at initialization. $\square$

**Lemma 4** (the jump rule loses nothing). Every vector of $R_m$ is dominated by a vector of $W$.

*Proof.* Let block $j$ be overfull in $k \in R_m$. Then $j < m$, since $t_m \ge 1$ forces $k_m = 0$. Because $t_j \ge 1$, Lemma 2(iii) gives a nonzero block to the left of $j$; let $\ell < j$ be the leftmost nonzero block. Then $t_\ell = 0$ if $\ell \ge 1$, by Lemma 2(iii) again, and $t_i \ge 1$ for $\ell < i \le j$. Set $\Delta = k_j - s_j$, and change the trade counts to

$$t'_i = t_i + \Delta \, \frac{P_j}{P_i} \quad (\ell < i \le j), \qquad t'_i = t_i \ \text{otherwise}.$$

Here $P_j / P_i = d_{i+1} \cdots d_j$ is an integer. By the formulas of Lemma 2, using $d_{i+1} P_j / P_{i+1} = P_j / P_i$, the new vector $k'$ has $k'_j = k_j - \Delta = s_j$, has $k'_i = k_i$ for $\ell < i < j$, has $k'_\ell = k_\ell + \Delta P_j / P_\ell$, and agrees with $k$ at every other block. All its entries are nonnegative and $t'_m = t_m$, so $t'$ satisfies the constraints of Lemma 2(ii) and $k' \in R_m$.

$k'$ uses the same blocks as $k$, since $k'_j = s_j \ge 1$ and $k'_\ell > k_\ell > 0$. Its saturated weight at $j$ is $s_j$, as before; at $\ell$ it is no smaller; elsewhere it is unchanged. So $k'$ dominates $k$. Block $j$ is not overfull in $k'$. Each block $i$ with $\ell < i < j$ keeps its weight and still has $t'_i \ge 1$, so it is overfull in $k'$ exactly when it was in $k$. Block $\ell$ is not overfull, because $t'_\ell = t_\ell = 0$ or $\ell = 0$. The other blocks are unchanged. So $k'$ has one overfull block fewer than $k$. Repeating the construction reaches a vector of $W$, and dominance is transitive. $\square$

**Lemma 5** (the leftmost saturation rule). Without locked blocks, the leftmost saturation rule rejects no vector of $W$.

*Proof.* The rule rejects $k$ only if some saturated block $i$ lies to the left of a block $i'$ with $k_{i'} > s_{i'}$. In $W$, a block holding more than its size is the leftmost nonzero block, so $k_i = 0$. But $i$ is saturated, so $k_i \ge s_i \ge 1$, a contradiction. $\square$

**Lemma 6** (the fillable space rule). Let $k \in W$, with no locked block, and let $\ell$ be its leftmost nonzero block. The fillable space rule rejects $k$ exactly when some partially filled block $p$ satisfies

$$(k_\ell - s_\ell)\, P_\ell \;\ge\; P_p ,$$

and in that case the vector $k' = k - (P_p / P_\ell)\, e_\ell + e_p$ lies in $R_m$ and strictly dominates $k$.

*Proof.* In $W$ only block $\ell$ can hold more than its size, so the sum the rule computes has the single term $\max(0, k_\ell - s_\ell)\, \texttt{values}[\ell]$, and the quantity compared with 1 at $p$ is $\max(0, k_\ell - s_\ell)\, P_\ell / P_p$. It is at least 1 exactly when the displayed inequality holds, because $P_p > 0$. Suppose it holds. Then $k_\ell > s_\ell$, so $\ell$ is not partially filled, and $p \ne \ell$; since $p$ is nonzero and $\ell$ is the leftmost nonzero block, $p > \ell$. Also $p < m$: if $k_m \ne 0$, then $k = e_m$ by Lemma 1(c), so $\ell = m$, which contradicts $p > \ell$. Set $D = P_p / P_\ell = d_{\ell+1} \cdots d_p$, and change the trade counts to

$$t'_i = t_i - \frac{P_p}{P_i} \quad (\ell < i \le p), \qquad t'_i = t_i \ \text{otherwise}.$$

These are nonnegative: by Lemma 2(iii), $t_i P_i = \sum_{h < i} k_h P_h \ge k_\ell P_\ell \ge (s_\ell + D) P_\ell \ge D P_\ell = P_p$ for $\ell < i \le p$. By the formulas of Lemma 2, the new vector has $k'_p = k_p + 1$, has $k'_i = k_i$ for $\ell < i < p$, has $k'_\ell = k_\ell - D \ge s_\ell$, and agrees with $k$ elsewhere. All its entries are nonnegative and $t'_m = t_m$, so $k' \in R_m$ by Lemma 2(ii). It uses the same blocks as $k$. Its saturated weight at $\ell$ is $s_\ell$, as before, and at $p$ it is $k_p + 1 \le s_p$, one more than before; elsewhere it is unchanged. So $k'$ strictly dominates $k$. $\square$

**Theorem 1.** Without locked blocks, the vectors `_re_mesh_iterator` yields represent pairwise different terms, and these terms are exactly the maximal terms of $T_m$.

*Proof.* By Lemma 1(a), $T_m = R_m$, so the maximal terms of $T_m$ are those of $R_m$. Every yielded vector lies in $W$, by Lemma 3.

*No term is yielded twice.* Lemma 3 visits each vector once, so two yields are two different vectors of $W$. Suppose $w, w' \in W$ describe the same term: the same support and the same saturated weights. They have the same leftmost nonzero block $\ell$. At every other block neither holds more than its size, so there $w_j = \bar w_j = \bar w'_j = w'_j$. Both have value $P_m$ (Lemma 1(b)), so $w_\ell P_\ell = w'_\ell P_\ell$, and $w = w'$.

*Every yielded term is maximal.* Let $w$ be yielded, and suppose some vector of $R_m$ strictly dominates it. By Lemma 4 that vector is dominated by some $\tilde w \in W$, which therefore strictly dominates $w$. The two have the same support $S$; let $\ell = \min S$. For $j \in S \setminus \{\ell\}$ neither holds more than its size, so $w_j \le \tilde w_j$. Equal values give

$$(w_\ell - \tilde w_\ell)\, P_\ell = \sum_{j \ne \ell} (\tilde w_j - w_j)\, P_j ,$$

where every term on the right is nonnegative. If all of them are 0, then $w$ and $\tilde w$ agree away from $\ell$, hence also at $\ell$, and describe the same term, contradicting strict dominance. Otherwise there is a block $p$ with $\tilde w_p > w_p$. It is partially filled in $w$, because $0 < w_p < \tilde w_p \le s_p$. The right-hand side is at least $P_p > 0$, so $w_\ell > \tilde w_\ell$. Since $\min(w_\ell, s_\ell) \le \min(\tilde w_\ell, s_\ell)$, this forces $\tilde w_\ell \ge s_\ell$, so

$$(w_\ell - s_\ell)\, P_\ell \;\ge\; (w_\ell - \tilde w_\ell)\, P_\ell \;\ge\; P_p ,$$

and by Lemma 6 the fillable space rule rejects $w$, a contradiction.

*Every maximal term is yielded.* Let $k \in R_m$ describe a maximal term. By Lemma 4, some $w \in W$ dominates $k$, and by maximality $w$ describes the same term. The traversal visits $w$ (Lemma 3). The leftmost saturation rule does not reject it (Lemma 5). If the fillable space rule rejected it, Lemma 6 would give a vector strictly dominating $w$, hence strictly dominating $k$, contradicting maximality. So $w$ is yielded. $\square$

### Locked Blocks

`root_expressions(locked_list)` lets a falsy flag **lock** a block: the block's own roots are not added ([Root Expressions — Locked Blocks](../theory/Root%20Expressions%20and%20LC%20Estimation.md#locked-blocks)). Let $L \subseteq \{0, \ldots, m\}$ be the set of locked blocks and $U$ the rest. Monomial profiles take no lock flags, so this subsection concerns root expressions only.

**The locked recursion.** The default path skips `extend` for a locked block, so it does not add the block's own $\langle s_j \cdot 1 \rangle$, and it still composes the block's chaining (`_re_default`). A locked source therefore gets the empty expression $\emptyset$, and $\emptyset \cdot X = \emptyset$ because `RootExpression.__mul__` forms products of pairs of terms and there are none. The recursion becomes

$$\text{RE}^L_0 = \begin{cases} \langle s_0 \cdot 1 \rangle & 0 \in U \\ \emptyset & 0 \in L \end{cases}, \qquad \text{RE}^L_j = (\text{RE}^L_{j-1})^{d_j} + \begin{cases} \langle s_j \cdot 1 \rangle & j \in U \\ \text{nothing} & j \in L \end{cases},$$

with constants handled as before. Its unpruned term set $T^L_j$ satisfies $T^L_0 = \{ e_0 \}$ or $\emptyset$, and $T^L_j$ is the set of sums of $d_j$ vectors of $T^L_{j-1}$, together with $e_j$ when $j \in U$. This is the inductive step of Lemma 1 with the own term left out at a locked block.

**Lemma 7** (locked terms). $T^L_m = \{\, k \in R_m : k_j = 0 \text{ for every } j \in L \,\}$.

*Proof.* The induction proves, for every $j \le m$, that $T^L_j = \{\, k \in R_j : k_i = 0 \text{ for every } i \in L \text{ with } i \le j \,\}$; the lemma is the case $j = m$. For $j = 0$, $R_0 = \{ e_0 \}$, and $e_0$ vanishes on the locked blocks up to 0 exactly when $0 \notin L$, which is when $T^L_0 = \{ e_0 \}$. Let $j \ge 1$. As in the proof of Lemma 1(a), $R_j = \{ e_j \} \cup Q$, where $Q$ is the set of sums of $d_j$ vectors of $R_{j-1}$. The vector $e_j$ vanishes on the locked blocks up to $j$ exactly when $j \notin L$. Every vector of $Q$ has $k_j = 0$, by Lemma 1(c). A sum of nonnegative vectors vanishes on a set of coordinates exactly when every summand does, so the vectors of $Q$ that vanish on the locked blocks below $j$ are exactly the sums of $d_j$ vectors of $T^L_{j-1}$, by the induction hypothesis. These are the two parts of $T^L_j$. $\square$

In terms of trades, a vector of $R_m$ vanishes on the locked blocks exactly when every unit that arrived in a locked block was traded on. Trading weight into a locked block and on through it is the expansion that skips that block's own term while still passing its factors upstream.

The mesh runs on all of $R_m$, not on $T^L_m$. It receives the sizes $s'_j = 0$ for $j \in L$ and $s'_j = s_j$ for $j \in U$ (`sizes *= locked_list`), so `num_traded` $= \max(1, k_j - 0)$ moves a locked block's whole holding in its first trade. On output, `re_compute_single_mesh` drops every field of size 0 and skips a vector that is left with no field (`if not field: continue`). The driver writes the empty expression for a locked source. So the output for a yielded vector $y$ is the term of its **restriction** $\rho(y)$, the vector of its entries at the blocks of $U$, and nothing when $\rho(y) = 0$. Comparisons between restrictions use dominance on the blocks of $U$ with their true sizes $s_j$; call it **$U$-dominance**. A vector of $T^L_m$ is zero on $L$ and nonzero, so its restriction is just the vector itself, and $U$-dominance on $T^L_m$ is ordinary dominance.

**Lemma 8** (locked weight is dominated). Every $k \in R_m$ with $\rho(k) \ne 0$ is $U$-dominated by a vector of $T^L_m$.

*Proof.* By induction on $j$, for $R_j$ and $T^L_j$, with $\rho$ restricted to the blocks of $U$ up to $j$. For $j = 0$, $R_0 = \{ e_0 \}$, and $\rho(e_0) \ne 0$ means $0 \in U$, so $e_0 \in T^L_0$. Let $j \ge 1$ and $k \in R_j$ with $\rho(k) \ne 0$. If $k = e_j$, then $j \in U$ and $e_j \in T^L_j$. Otherwise $k = x_1 + \cdots + x_{d_j}$ with every $x_r \in R_{j-1}$ (proof of Lemma 1(a)). Let $N$ be the set of $r$ with $\rho(x_r) \ne 0$; it is nonempty, because $\rho(k) = \sum_r \rho(x_r)$ and all entries are nonnegative. For $r \in N$, the induction hypothesis gives $x_r^* \in T^L_{j-1}$ that $U$-dominates $x_r$. Fix $r_0 \in N$, set $y_r = x_r^*$ for $r \in N$ and $y_r = x_{r_0}^*$ for $r \notin N$, and let $k^* = y_1 + \cdots + y_{d_j}$, which lies in $T^L_j$. Its blocks are $\bigcup_{r \in N} \mathrm{supp}(x_r^*) = \bigcup_{r \in N} \mathrm{supp}(\rho(x_r)) = \mathrm{supp}(\rho(k))$, because $\rho(x_r) = 0$ for $r \notin N$. At a block $i \in U$, each $\min(s_i, y_{r,i})$ is at least $\min(s_i, x_{r,i})$: for $r \in N$ by the induction hypothesis, and for $r \notin N$ because $x_{r,i} = 0$. By the identity of [Terms and Dominance](#terms-and-dominance), the saturated weight of a sum is $\min(s_i, \sum_r \min(s_i, \cdot))$ of its summands, which is monotone in each, so $k^*$ has saturated weight at least that of $k$ at every block of $U$. $\square$

**Corollary 9.** The maximal terms of $T^L_m$ are exactly the maximal $U$-terms among the vectors $k \in R_m$ with $\rho(k) \ne 0$.

*Proof.* $T^L_m \subseteq R_m$, and its vectors have nonzero restriction. Let $k \in T^L_m$ be maximal in $T^L_m$. If some $x \in R_m$ strictly $U$-dominated $k$, Lemma 8 would give $x^* \in T^L_m$ that $U$-dominates $x$, hence strictly dominates $k$, a contradiction. Conversely, let $x \in R_m$ with $\rho(x) \ne 0$ be $U$-maximal. Lemma 8 gives $x^* \in T^L_m$ that $U$-dominates it, so $x^*$ describes the same term, and $x^*$ is maximal in $T^L_m$ because $T^L_m \subseteq R_m$. $\square$

With sizes $s'$, call a block $j \ge 1$ overfull when $t_j \ge 1$ and $k_j > s'_j$, and let $W'$ be the vectors of $R_m$ with no overfull block. A locked block $j \ge 1$ is overfull exactly when it has traded and still holds any weight.

**Lemma 10** (Lemmas 3, 4, and 6 under locking).

- (a) The traversal visits each vector of $W'$ exactly once and no other vector.
- (b) Every vector of $R_m$ is $U$-dominated by a vector of $W'$.
- (c) Let $k \in W'$ with leftmost nonzero block $\ell$. The fillable space rule rejects $k$ exactly when some block $p$ with $0 < k_p < s'_p$ satisfies $(k_\ell - s'_\ell) P_\ell \ge P_p$, and in that case $k' = k - (P_p / P_\ell)\, e_\ell + e_p$ lies in $R_m$ and strictly $U$-dominates $k$.

*Proof.* (a) The proof of Lemma 3 never uses $s_j \ge 1$. With $s'_j = 0$ it says that $t_j$ takes the values 0 and $I$, which are the values for which block $j$ is not overfull.

(b) Run the construction of Lemma 4 with $s'$. The step "$j < m$" holds as before, and Lemma 2 does not involve sizes. What changes is the check that $k'$ dominates $k$, which now concerns the blocks of $U$ only. At the overfull block $j$: if $j \in U$, then $k'_j = s_j \ge 1$, so $j$ stays in the support with the same saturated weight $s_j$; if $j \in L$, then $k'_j = 0$, and $j$ is not a block of $U$. At $\ell$: if $\ell \in U$, then $k'_\ell > k_\ell > 0$, so its saturated weight does not decrease; if $\ell \in L$, it is not a block of $U$. Every other block keeps its weight. The count of overfull blocks falls by one as before.

(c) The proof of Lemma 6 goes through with $s'$, since its inequality $k'_\ell = k_\ell - D \ge s'_\ell$ and its trade-count bounds do not use $s'_\ell \ge 1$. A block with $s'_p = 0$ is never partially filled, so $p \in U$, and $k'_p = k_p + 1 \le s_p$ raises the saturated weight at $p$ by one. At $\ell$: if $\ell \in U$, then $k'_\ell \ge s_\ell \ge 1$, so its saturated weight stays $s_\ell$; if $\ell \in L$, it is not a block of $U$. The blocks between are unchanged. $\square$

**Lemma 11** (the leftmost saturation rule under locking). Let $k \in W'$ with $k \ne e_m$, and let $\ell$ be its leftmost nonzero block. The rule rejects $k$ exactly when $k_\ell > s'_\ell$ and some locked block lies to the left of $\ell$. In that case, let $i^*$ be the leftmost locked block and $\Delta = k_\ell - s'_\ell$, and change the trade counts to $t''_h = t_h + \Delta P_\ell / P_h$ for $i^* < h \le \ell$. The resulting vector $k''$ lies in $W'$, describes the same $U$-term as $k$, and is not rejected by the rule.

*Proof.* Since $k \ne e_m$, Lemma 1(c) gives $k_m = 0$, so $\ell < m$. In $W'$, a block holding more than its size is the leftmost nonzero block (as for $W$, by Lemma 2(iii)), so the rule can reject $k$ only because $k_\ell > s'_\ell$ and some block to the left of $\ell$ is saturated. Those blocks are empty, and an empty block is saturated exactly when its size is 0, that is, when it is locked. This proves the first statement.

For the second, $i^* < \ell$ because a locked block lies left of $\ell$. By the formulas of Lemma 2, as in Lemma 4, $k''_\ell = k_\ell - \Delta = s'_\ell$, $k''_h = k_h = 0$ for $i^* < h < \ell$, $k''_{i^*} = \Delta P_\ell / P_{i^*}$, and every other entry is unchanged. All entries are nonnegative and $t''_m = t_m$, so $k'' \in R_m$. No block of $k''$ is overfull: blocks $h$ with $i^* < h \le \ell$ hold at most their size; $t''_{i^*} = t_{i^*} = 0$ if $i^* \ge 1$, because the blocks left of $i^*$ are empty; and the blocks right of $\ell$ are unchanged. So $k'' \in W'$. The $U$-term is unchanged: if $\ell \in U$, its saturated weight was $s_\ell$ (since $k_\ell > s_\ell$) and is $s_\ell \ge 1$; if $\ell \in L$, it is not a block of $U$; the blocks between stay empty; and $i^*$ is locked. Finally, the only block of $k''$ holding more than its size is $i^*$, and every block left of $i^*$ is empty and unlocked, hence of size at least 1 and not saturated. So the rule does not reject $k''$. $\square$

**Theorem 2.** For any set $L$ of locked blocks:

- (i) every output of `re_compute_single_mesh` (the term of $\rho(y)$ for a yielded $y$ with $\rho(y) \ne 0$) is dominated by a term of $T^L_m$;
- (ii) every maximal term of $T^L_m$ is an output.

Consequently the union of the outputs' root sets equals the union of the root sets of the terms of $T^L_m$.

*Proof.* (i) is Lemma 8, since every yielded vector lies in $R_m$.

(ii) Let $x \in T^L_m$ describe a maximal term. By Corollary 9, it is a maximal $U$-term among the vectors of $R_m$ with nonzero restriction. Lemma 10(b) gives $w \in W'$ that $U$-dominates $x$, and by maximality $w$ describes the same $U$-term. If $w = e_m$, then $\rho(e_m) \ne 0$, so $m \in U$, and $e_m$ is yielded at initialization. Otherwise, if the leftmost saturation rule rejects $w$, replace $w$ by the vector $w''$ of Lemma 11, which lies in $W'$, describes the same $U$-term, and is not rejected by that rule. The traversal visits it (Lemma 10(a)). If the fillable space rule rejected it, Lemma 10(c) would give a vector of $R_m$ strictly $U$-dominating it, contradicting maximality (Corollary 9). So it is yielded, its restriction is nonzero, and its output is the term of $x$.

For the consequence: a dominated term's roots are among the dominating term's, so by (i) every output's roots lie in the union for $T^L_m$; every term of $T^L_m$ is dominated by a maximal one, which by (ii) is an output. $\square$

Theorem 2 does not show that every output is maximal, or that no term is output twice, when blocks are locked. Neither affects the union of root sets, which is what `upper()` measures when field degrees are pairwise coprime (see the note after [Proposition 5](../theory/Root%20Expressions%20and%20LC%20Estimation.md#the-five-root-expression-propositions)). The same holds for `lower()` (`RootExpression.py`, lines 220–250). It files each term under its tuple of fields, without splitting full classes, and measures, for each tuple, the union of origin-anchored boxes whose sides are the term's largest Jordan length and `pessimistic_expected_value(b, c)` for each field $b$ with weight $c$. That side does not decrease as $c$ grows (checked numerically for every $b < 40$), so the box of a dominated term lies inside the box of the term dominating it, under the same tuple. By Theorem 2 (i) and (ii), each tuple's union of boxes is then the same on both paths, and so is `lower()`, up to the order in which floating-point volumes are summed; that is the difference the 576-case probe below observed. Theorem 1 gives the same conclusion without locking.

**Example.** Take the [Pruning Example](#pruning-example) with the source locked, so the iterator receives sizes $s' = [0, 3, 2, 7]$ (the values are still $1, 2, 6, 12$). By Lemma 7, $T^L_3 = \{ (0,0,0,1), (0,0,2,0), (0,3,1,0), (0,6,0,0) \}$, the vectors of the earlier table with $k_0 = 0$. Their terms $\langle 7 \cdot 1 \rangle$, $\langle 2 \cdot 2 \rangle$, $\langle 3 \cdot 3 \rangle \times \langle 2 \cdot 1 \rangle$ and $\langle 3 \cdot 3 \rangle$ use four different sets of fields, so all four are maximal. With $s'_0 = 0$, all of block 0's weight is excess, and the run goes as follows.

- The iterator yields, in order, $(0,0,0,1)$, $(0,0,2,0)$, $(0,3,1,0)$, $(6,3,0,0)$ and $(12,0,0,0)$.
- The fillable space rule rejects $(2,2,1,0)$, whose excess 2 at block 0 is worth $2/P_1 = 1$ unit at the partially filled block 1, and $(4,1,1,0)$, worth $4/2 = 2$ units there. It also rejects $(6,0,1,0)$: its excess 6 is worth exactly $6/6 = 1$ unit at the partially filled block 2, and the code compares with $\ge 1$. Without locking this vector was kept, because its excess was $6 - 5 = 1$, worth only $1/6$. It rejects $(8,2,0,0)$ and $(10,1,0,0)$ as before.
- The leftmost saturation rule rejects $(0,6,0,0)$: block 0 is saturated ($0 \ge 0$) and block 1 holds $6 > 3$. Lemma 11's replacement, with $i^* = 0$, $\Delta = 3$ and $t''_1 = 3$, is $(6,3,0,0)$, which is yielded and has the same restriction term $\langle 3 \cdot 3 \rangle$.
- $(12,0,0,0)$ is yielded but produces no output: its only weight is in the locked block 0, so `re_compute_single_mesh` skips it (`if not field: continue`).

The five yields give four outputs, $\langle 7 \cdot 1 \rangle$, $\langle 2 \cdot 2 \rangle$, $\langle 3 \cdot 3 \rangle \times \langle 2 \cdot 1 \rangle$ and $\langle 3 \cdot 3 \rangle$, exactly the maximal terms of $T^L_3$. In this instance the outputs equal the maximal terms exactly; Theorem 2 claims less in general. This run was confirmed by executing `_re_mesh_iterator` during review.

### What the Tests Check

The theorems concern terms. The tests in `test/Tools/linear_complexity_test.py` compare counts.

- `test_mesh_path_agrees_with_the_default_path` checks that the mesh and default paths give the same `upper()` and `lower()` on every bit for registers of sizes $(5, 3)$, $(7, 4, 3)$, and $(7, 5, 3, 2)$.
- `test_mesh_path_agrees_with_the_default_path_under_locking` does the same for every pattern of locked blocks on sizes $(7, 4, 3)$, with chaining from `arman_template`. That template gives each bit's function a `CONST(1)` with probability 1/2 (`arman_template.py`, line 20), which would hide a wrongly added constant term.
- `test_mesh_path_locking_without_chaining_constants` therefore repeats the comparison for all 8 lock patterns on sizes $(7, 4, 3)$ with hand-written chaining that has no constant.
- `test_mesh_path_requires_pairwise_coprime_sizes` checks that sizes $(6, 4)$ take the default path and sizes $(5, 4)$ the mesh.
- `test_chaining_nonlinear_in_own_bits_is_rejected` checks that `root_expressions()`, `monomial_profiles()`, and `estimate_LC` raise `ValueError` on a $(7, 4, 3)$ register whose middle block's chaining multiplies an own bit by a source bit.
- `test_monomial_profile_mesh_path_agrees_with_the_default_path` and `test_monomial_profile_mesh_path_with_repeated_and_one_bit_blocks` check equal monomial-profile `upper()` values, the second with repeated sizes and a one-bit block.

A larger probe, run during review and not recorded as a test, compared the two paths over 576 combinations: registers of sizes $(7, 5, 3, 2)$, $(9, 7, 5, 4)$, $(11, 9, 7, 5, 4)$ and $(7, 5, 4, 3, 2)$, chaining from `arman_template` with `max_and` 2, 3, and 4, two seeds, and every lock pattern. The upper bounds agreed in all 576. The lower values agreed up to floating-point rounding: in 18 cases, with values of order $10^6$ to $10^9$, they differed only beyond the displayed digits, consistent with the rectangle volumes being summed in a different order.

## Output Construction

Each yielded weight vector is converted to a `JordanSet` or `TermSet`:

- **Root expressions** (`re_compute_single_mesh`): For each nonzero entry $k_j$ at a block of nonzero size, the field degree $s_j$ and coset weight $\min(s_j, k_j)$ are recorded; a vector left with no field is skipped. The Jordan multiplicity set is $\{1\}$. The default path agrees under pairwise coprime sizes: `RootExpression.extend` records a Jordan length above 1 for the own term $\langle s_i \cdot 1 \rangle$ only when `is_embedded_subset(⟨s_i·1⟩, term)` holds for some existing term, which requires the field $s_i$ among that term's fields, and no upstream term has it when the sizes are pairwise coprime. The upstream terms themselves keep Jordan set $\{1\}$, because a product of two Jordan sets takes its lengths from `JP_solve` (`JordanSet.py`, lines 28–31), and `JP_solve(1, 1, 2)` returns `[(1, 1)]`; for example, `JordanSet({3: 1}, {1}) * JordanSet({5: 1}, {1})` has Jordan set $\{1\}$ (checked during review).
- **Monomial profiles** (`mp_compute_single_mesh`): The analogous conversion to `TermSet` objects, recording block IDs and saturated weights.

After iteration, the CMPR driver handles two things the mesh does not:

1. **Constants:** If block $i$'s chaining can produce a constant term, or block $i - 1$'s expression contains the constant term (the term with no fields, not the root 1, which a full class $\langle s \cdot s \rangle$ also contains), the driver adds `RootExpression.logical_one()` (or `MonomialProfile.logical_one()`); these are the two ways a constant term enters, derived in [Constants and Blocks Without Chaining Input](#constants-and-blocks-without-chaining-input). A block's own constant terms are read from `CONST(1)` in its as-written profile during degree extraction, which over-approximates the true ANF's constant term ([Degree Extraction](#degree-extraction)), and stored in the `constants_possible` array; `constants_possible[block_id]` is also set whenever `constants_possible[block_id - 1]` is.
2. **Locked blocks:** The `locked_list` parameter zeros out the sizes of locked blocks (falsy flags; `sizes *= locked_list`), so the first trade out of a locked block moves its whole holding. On output a locked block contributes no field, and a vector whose weight lies only in locked blocks is skipped (`if not field: continue`). The driver writes the empty expression for a locked source. [Locked Blocks](#locked-blocks) proves what this computes.

## Eligibility Conditions

Before either path runs, `CMPR.root_expressions()` and `CMPR.monomial_profiles()` raise `ValueError` if any block's update multiplies its own bits (`_require_chaining_linear_in_own_bits()`; see [The Model Hypothesis](#the-model-hypothesis)). That is a precondition of both paths, not a condition for choosing the mesh. Updates that read a block's own bits linearly beyond its MPR update are not detected.

The mesh optimization is used when every condition below that applies to the method being called holds. The first three are checked in `CMPR.root_expressions()` and `CMPR.monomial_profiles()`, the last in the mesh drivers `_re_mesh_optimization` and `_mp_mesh_optimization`.

| Condition | Required for | Reason |
|---|---|---|
| No 1-bit MPR blocks | Root expressions | Size-1 blocks produce roots in $\mathbb{F}_2$ ($\alpha = 1$ only). The default path records a 1-bit block's own root as the term with no fields (`JordanSet({}, {1})` in `_re_default`), the same term a constant produces, not as a field of its own. See the open question below the table. |
| Pairwise coprime block sizes | Root expressions | The counts of different fields multiply only when the fields meet in $\mathbb{F}_2$, which is the case exactly for pairwise coprime sizes (Proposition 3); sizes such as 2 and 4 share the subfield $\mathbb{F}_4$. See [Why Root Expressions Need Pairwise Coprime Sizes](#why-root-expressions-need-pairwise-coprime-sizes). |
| Simple chaining | Both | Each block's update functions may read only the block's own bits and those of the immediately preceding block, and the source only its own (`used_bits <= allowed_bits`). Together with the rejection above, this is still weaker than (H), since linear own-bit terms pass both ([The Model Hypothesis](#the-model-hypothesis)). Under (H), each block's expression is a power of the preceding block's alone ([The Recursion](#the-recursion)), with the single exponent $D_i$; chaining from several upstream blocks would need an exponent per upstream block. |
| Every block after the source takes input from upstream | Both | Such a block has $D_i \ge 1$. A block with $D_i = 0$ has extracted degree 1, which would make the mesh trade through it as if it chained linearly ([Degree Extraction](#degree-extraction)). |

> **OPEN:** Is the one-bit exclusion about the repeated root 1, about the representation, or both? — needs author input.
> The default path adds a 1-bit block's own root through `RootExpression.extend` with `JordanSet({}, {1})`. `is_embedded_subset` holds for that set against any existing term all of whose classes are full (every field outside the empty key set must have weight equal to its degree), and against the constant term, and `extend` then records Jordan length 2. The mesh's Jordan sets are always $\{1\}$. Which of these differences the exclusion guards against is not recorded.

Monomial profiles need only the last two conditions. Monomials are indexed by register, not by field, so block sizes place no condition on them (the comment in `CMPR.monomial_profiles()`, and [Why Root Expressions Need Pairwise Coprime Sizes](#why-root-expressions-need-pairwise-coprime-sizes)). `test_monomial_profile_mesh_path_with_repeated_and_one_bit_blocks` runs repeated sizes and a one-bit block through both paths.

When any condition fails, `CMPR` falls back to the default gate-by-gate composition algorithm (`_re_default` / `_mp_default`), which handles Jordan interactions and arbitrary chaining topologies at higher computational cost. When block sizes share a factor, its root count is still an upper bound, but it can count a root of a common subfield more than once: for sizes 4 and 2 it lists both $\langle 4 \cdot 2 \rangle$ and $\langle 2 \cdot 1 \rangle$, and the two elements of $\mathbb{F}_4$ other than 0 and 1 lie in both.

## Performance

The default composition algorithm has cost $O(nk^2 N^4)$, where $n$ is the size of the register in bits, $k$ is the number of blocks ($m + 1$ in the numbering above, where the blocks are $0, \ldots, m$; not the weights $k_j$), and $N \leq \prod s_i / \max(s_i)$ is the maximum number of terms in any root expression (see [Root Expressions — Runtime Complexity](../theory/Root%20Expressions%20and%20LC%20Estimation.md#runtime-complexity-of-root-counting)). The gate-by-gate composition builds many intermediate expression objects and performs expensive dominance-pruning at each step.

The mesh avoids intermediate objects entirely — it generates the final set of maximal terms directly via the numba-jitted iterator. The pruning is integrated into the tree traversal (checked before yielding) rather than applied as a post-processing pass over $O(N^2)$ term pairs.

> **OPEN:** How much faster is the mesh than the default path, and on which registers? — needs author input.
> An earlier version of this paragraph said "orders of magnitude faster" for CMPRs with many blocks of pairwise coprime sizes and high chaining degrees. No timing of the two paths is recorded in the repository.

## Connections

- [Root Expressions and LC Estimation](../theory/Root%20Expressions%20and%20LC%20Estimation.md) — the five propositions that the mesh relies on (especially Props 2 and 3), and the propagation rule from which the recursion is derived
- [Monomial Profile Theory](../theory/Monomial%20Profile%20Theory.md) — the propositions the monomial-profile mesh relies on, under which block sizes place no condition
- [Coset Regions](../conventions/Coset%20Regions.md) — the geometric picture of weight distributions as grid positions
- [Root Multiplicities and Jordan Decomposition](../theory/Root%20Multiplicities%20and%20Jordan%20Decomposition.md) — why same-field products (repeated sizes) require Jordan partitions that the mesh cannot provide
