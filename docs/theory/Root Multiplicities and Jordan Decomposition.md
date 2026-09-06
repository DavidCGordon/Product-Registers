# Root Multiplicities and Jordan Decomposition

Root expressions describe possible roots of a register sequence, but a root alone is not enough when the update operator is not semisimple. A repeated eigenvalue can occur in a Jordan block, and the associated generalized eigenspace contributes a polynomial factor in time as well as an exponential factor. This document explains the multiplicity data carried by `JordanSet` and why multiplication uses Jordan partitions.

## From Eigenvalues to Sequences

Let $U$ be the linear update map of a register, extended from $\mathbb{F}_2$ to the algebraic closure $\overline{\mathbb{F}_2}$. On a Jordan block of size $m$ with eigenvalue $\alpha$, write

$$
J_m(\alpha)=\alpha I+N,
$$

where $N$ is nilpotent and $N^m=0$. Powers of the block have the form

$$
J_m(\alpha)^t=\alpha^t(I+\alpha^{-1}N)^t.
$$

Consequently, an observed coordinate is a linear combination of terms

$$
\binom{t}{r}\alpha^t, \qquad 0\le r<m,
$$

over the field of characteristic $2$. The eigenvalue $\alpha$ identifies the exponential root, while the Jordan size controls the generalized-root multiplicity. In characteristic $2$, the binomial coefficients can vanish in patterns that depend on the binary expansion of $t$; therefore the multiplicity bookkeeping is not an ordinary count of repeated copies.

For a Frobenius orbit of $\alpha$, all conjugate roots occur together in a sequence over $\mathbb{F}_2$. A root expression therefore tracks field-degree/root coverage and the possible Jordan lengths attached to that coverage.

## What `JordanSet` Stores

A `JordanSet` has two pieces of data:

- `roots` maps a field degree $b$ to a coverage value $c$. The pair $(b,c)$ corresponds to the coset class $\langle b \cdot c \rangle$ — roots in $\mathbb{F}_{2^b}$ whose cyclotomic coset has weight at most $c$ (see [Root Expressions — Coset Classes](Root%20Expressions%20and%20LC%20Estimation.md#coset-classes)).
- `mults` is a set of possible Jordan block lengths. It is a set of lengths, not a scalar multiplicity and not the number of roots in the coset class.

A `RootExpression` is implemented as table of terms, grouped by the field degrees appearing in each term. It is an upper-bound object: addition and multiplication preserve possible contributions, even when a particular register or initial state causes cancellation.

## Jordan Decomposition and Products

Suppose two generalized eigenspaces have eigenvalues $\alpha$ and $\beta$, with Jordan lengths $s$ and $t$. On their tensor product, the eigenvalue is $\alpha\beta$. The nilpotent part has a decomposition into Jordan blocks whose lengths are given by the Jordan-partition calculation. In the implementation this is `JP_solve(s, t, 2)`.

Thus, for two `JordanSet` objects $A$ and $B$,

$$
A\cdot B:
\quad
c_b=\min(b,c_b^A+c_b^B),
\qquad
M_{A\cdot B}=\bigcup_{s\in M_A,\,t\in M_B}J_2(s,t),
$$

where $J_2(s,t)$ is the set of Jordan lengths returned by the characteristic-two Jordan partition. The truncation at $b$ records that coset weight cannot exceed the field degree (see [Coset Regions](../conventions/Coset%20Regions.md) for the geometric picture).

This explains why same-field multiplication cannot be represented by simply adding scalar multiplicities. For distinct component fields, the eigenvalue choices are independent and the root counts multiply; for equal fields, generalized eigenspaces interact through tensor products and require the Jordan partition.

## Counting and the Upper Bound

The `upper()` method expands full coset classes (weight $= e$) into their embedded subfield classes, forms the corresponding counting regions (see [Coset Regions](../conventions/Coset%20Regions.md)), and uses the largest possible Jordan length in each term. This is intentionally optimistic: it assumes every allowed root/coset and every allowed generalized-root layer is present. The result is an upper bound on trajectory linear complexity, not a claim that all listed contributions survive in a concrete sequence.

The `lower()` method uses the same structure with pessimistic expected values, capping each coset class's weight at $e - 1$ to exclude the weight-$e$ roots that cancel with probability $1/2$ (see [Root Expressions — When the Upper Bound is Not Tight](Root%20Expressions%20and%20LC%20Estimation.md#when-the-upper-bound-is-not-tight)). It is statistical rather than an algebraic lower bound for every initial state.

## Relation to Root Expressions

The five combination rules in [Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md) are the counting-level form of this decomposition:

1. a single coset class is counted by its available Frobenius cosets;
2. same-field products add coverage and combine Jordan lengths through Jordan partitions;
3. independent field components multiply their available choices;
4. XOR combines possible contributions and corrects overlaps;
5. intersections retain only common field coverage and common compatible Jordan lengths.

These rules are exact for the represented upper-bound sets when their stated hypotheses hold. They do not remove cancellations caused by a particular chaining function or initial state.

## Multiplicity Brackets and the Period Staircase

### The Frobenius Identity and Power-of-2 Periodicity

In characteristic 2, the Frobenius endomorphism gives

$$(1 + D)^{2^k} = 1 + D^{2^k}.$$

This identity controls the period of sequences annihilated by powers of $(1 + D)$. If a sequence satisfies $(1 + D)^m s[t] = 0$, then $(1 + D)^{2^k} s[t] = (1 + D^{2^k}) s[t] = 0$ for the smallest $k$ with $2^k \geq m$. The minimal polynomial therefore divides $(1 + D^{2^k})$, and the sequence has period dividing $2^k$. In particular:

$$\text{period} = 2^{\lceil \log_2 m \rceil}$$

where $m$ is the maximum Jordan block length for root $\alpha = 1$. Jordan lengths do not produce arbitrary periods — they produce periods that jump in powers of 2. We call the range $[2^{k-1}+1,\; 2^k]$ the **$k$-th multiplicity bracket**: any Jordan length in this range yields the same period $2^k$.

### The Power-of-2 Ceiling on Products

Tensor products of Jordan blocks in characteristic 2 cannot escape the current bracket. Specifically:

> **Theorem (Product Ceiling):** If $s \leq 2^k$ and $t \leq 2^k$, then every Jordan block in the decomposition of $J_s \otimes J_t$ (over a field of characteristic 2) has length $\leq 2^k$.

This is verified computationally via `JP_solve(s, t, 2)`:

| Product | Decomposition | Max length | Bracket |
|---|---|---|---|
| $J_2 \otimes J_2$ | $2 \times J_2$ | 2 | $2^1$ |
| $J_4 \otimes J_4$ | $4 \times J_4$ | 4 | $2^2$ |
| $J_8 \otimes J_8$ | $8 \times J_8$ | 8 | $2^3$ |

But once either input exceeds the bracket ceiling:

| Product | Decomposition | Max length | Bracket |
|---|---|---|---|
| $J_2 \otimes J_3$ | $J_4 + J_2$ | 4 | $2^2$ |
| $J_4 \otimes J_5$ | $J_8 + 3 \times J_4$ | 8 | $2^3$ |

The product of two inputs within bracket $k$ stays in bracket $k$, but the product of an input at the ceiling $2^k$ with an input exceeding the ceiling can reach the next bracket $2^{k+1}$.

The underlying reason is the Frobenius identity: the sequence $\binom{t}{r} \bmod 2$ for $r < 2^k$ satisfies $(1+D)^{2^k} = 1 + D^{2^k}$, and convolution of two such sequences (the tensor product) cannot generate components annihilated by a higher power than $(1+D)^{2^k}$.

### The Resolvent Pushes Past the Ceiling

The resolvent adds exactly 1 to the maximum Jordan length (see [Resolvent Analysis](Resolvent%20Analysis.md) for the derivation). In the D-transform framework:

$$B(D) = \frac{1}{\chi_U(D)} \operatorname{Adj}(I \oplus UD)\bigl(D\,\mathcal{C}(D) \oplus B[0]\bigr)$$

When the block has characteristic polynomial $\chi_U(D) = (1+D)^s$ and the chaining input has maximum Jordan length $m$, the driven output's maximum Jordan length is $m + 1$. Unlike the product (which stays within the bracket), the resolvent increment can push a Jordan length from $2^k$ to $2^k + 1$, crossing into bracket $k+1$ and doubling the period.

This asymmetry — products preserve the ceiling, the resolvent breaks it — is the engine of period growth in chained registers.

### T-Functions and Binomial Sequences

T-functions provide the clearest illustration because they eliminate all root structure except multiplicities. Every block uses the polynomial $1 + x$, whose only root is $\alpha = 1$. The entire output is therefore a linear combination of **binomial sequences**:

$$b_r[t] = \binom{t}{r} \bmod 2, \qquad r = 0, 1, 2, \ldots$$

Each $b_r$ is generated by a Jordan block $J_{r+1}(1)$: the sequence $(1+D)^{r+1} b_r[t] = 0$ but $(1+D)^r b_r[t] \neq 0$. The period of $b_r$ is $2^{\lceil \log_2(r+1) \rceil}$, matching the bracket structure. All complexity in a T-function's output comes from which binomial layers $b_r$ are present — which is exactly the Jordan multiplicity set tracked by `JordanSet.mults`.

### Period Doubling Through Chaining

Consider an $n$-bit T-function where bit $k$ has update $b_k[t+1] = b_k[t] \oplus f_k(b_0[t], \ldots, b_{k-1}[t])$. The initial bit (bit 0) has no chaining inputs, but the T-function construction adds a constant driving term $\oplus\, 1$ — the "+1" of the binary counter. Since this constant shares root $\alpha = 1$ with the block's polynomial $(1+x)$, the resolvent interaction pushes the Jordan length from 1 to 2: bit 0 alternates with period 2, not period 1. Its Jordan set is $\{1, 2\}$.

Each successive bit's Jordan set is then determined by:

1. **Chaining:** The function $f_k$ combines earlier bits' Jordan sets via products (AND) and unions (XOR). By the ceiling theorem, these products cannot exceed the current bracket ceiling.
2. **Resolvent:** The block's own $(1+D)$ factor adds 1 to the maximum, potentially pushing past the ceiling.

The period doubles at bit $k$ if and only if the chaining products reach the current bracket ceiling $2^j$, so that the resolvent increment $2^j \to 2^j + 1$ crosses into the next bracket. With full nonlinear mixing, the Jordan set evolution is:

| Bit | Chaining max | Resolvent max | Bracket | Period |
|---|---|---|---|---|
| 0 | 1 (constant) | 2 | $2^1$ | 2 |
| 1 | 2 | 3 | $2^2$ | 4 |
| 2 | 4 | 5 | $2^3$ | 8 |
| 3 | 8 | 9 | $2^4$ | 16 |
| $k$ | $2^k$ | $2^k+1$ | $2^{k+1}$ | $2^{k+1}$ |

At each step, the chaining products saturate the current bracket (achieving max $= 2^k$), and the resolvent pushes to $2^k+1$, entering the next bracket with double the period. The general formula is: bit $k$ of a full-period T-function has maximum Jordan length $2^k + 1$ and period $2^{k+1}$.

If the chaining at bit $k$ fails to reach the ceiling, the resolvent stays within the current bracket and the period does not increase.

### Why Every Prior Bit Is Required

In the full-mixing T-function, the chaining at bit $k$ is the product of all bits $0, \ldots, k-1$, and every one of them is essential. Dropping any single bit $j$ from the product reduces the maximum Jordan length by exactly $2^j$:

| Dropped bit $j$ | Max Jordan of bit $j$ | Product max without $j$ | Shortfall |
|---|---|---|---|
| 0 | 2 | $2^k - 1$ | $2^0 = 1$ |
| 1 | 3 | $2^k - 2$ | $2^1 = 2$ |
| 2 | 5 | $2^k - 4$ | $2^2 = 4$ |
| $j$ | $2^j + 1$ | $2^k - 2^j$ | $2^j$ |

In every case, $2^k - 2^j < 2^k$, so the product without bit $j$ fails to reach the ceiling and the period does not double. The reason this is tight is the binary sum identity:

$$\sum_{j=0}^{k-1} 2^j = 2^k - 1$$

Each bit $j$ contributes $2^j$ to the product maximum. The contributions sum to $2^k - 1$, and the base (the identity element of the product, with max 1) brings the total to $2^k$. Since the sum is tight — there is no slack — every term is load-bearing, and omitting any one leaves the product strictly below the ceiling.

This gives a Jordan-theoretic proof of the T-function full periodicity condition: **the chaining at bit $k$ must include all prior bits** (with sufficient nonlinear mixing) **for the period to double**. In the standard T-function construction, this is equivalent to requiring each $f_k$ to depend nontrivially on $b_{k-1}$, since $b_{k-1}$ is the unique bit whose Jordan set exceeds the running product's bracket ceiling.

### Connection to $(1+x)^n$ Recurrences

A T-function bit $k$ with maximum Jordan length $m = 2^k + 1$ satisfies the recurrence $(1+x)^m \cdot b_k[t] = 0$. By the Frobenius identity, $(1+x)^{2^{k+1}} = 1 + x^{2^{k+1}}$, and since $m \leq 2^{k+1}$, the minimal polynomial divides $(1+x)^{2^{k+1}}$, giving period $2^{k+1}$. For the last bit of an $n$-bit T-function ($k = n-1$), the maximum Jordan length is $2^{n-1} + 1$, the annihilating polynomial is $(1+x)^{2^{n-1}+1}$, and the period is $2^n$.

## Practical Consequences

Repeated component sizes are the main reason the [mesh optimization](../architecture/Mesh%20Optimization.md) is disabled in some CMPR configurations: equal-sized components share field structure, so their products need multiplicity-aware Jordan handling. T-Functions are especially exposed because every component uses the same polynomial $(1+x)$, making the bracket/staircase structure (above) the dominant factor in period and linear complexity analysis. For those cases, use the default root-expression path and treat the resulting bounds as implementation-backed estimates whose tightness still requires empirical checking.
