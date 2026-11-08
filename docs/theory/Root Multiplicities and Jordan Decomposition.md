# Root Multiplicities and Jordan Decomposition

[Root expressions](Root%20Expressions%20and%20LC%20Estimation.md) describe possible roots of a register sequence, but a root alone is not enough when the update operator is not semisimple. A repeated eigenvalue can occur in a Jordan block, and the associated generalized eigenspace contributes a polynomial factor in time as well as an exponential factor. This document explains the multiplicity data carried by `JordanSet` (`src/PyPR/Tools/RootCounting/JordanSet.py`) and why multiplication uses Jordan partitions.

It assumes the root and Frobenius-orbit background of [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md) and the coset classes $\langle e \cdot w \rangle$ of [Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md#coset-classes).

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

- `roots` maps a field degree $b$ to a coverage value $c$. The pair $(b,c)$ corresponds to the coset class $\langle b \cdot c \rangle$ — roots in $\mathbb{F}_{2^b}$ whose cyclotomic coset has weight from 1 to $c$ (see [Root Expressions — Coset Classes](Root%20Expressions%20and%20LC%20Estimation.md#coset-classes)).
- `mults` is a set of possible Jordan block lengths. It is a set of lengths, not a scalar multiplicity and not the number of roots in the coset class.

A `RootExpression` is implemented as a table of terms, grouped by the field degrees appearing in each term. It is an upper-bound object: addition and multiplication preserve possible contributions, even when a particular register or initial state causes cancellation.

## Jordan Decomposition and Products

Suppose two generalized eigenspaces have eigenvalues $\alpha$ and $\beta$, with Jordan lengths $s$ and $t$. On their tensor product, the eigenvalue is $\alpha\beta$. The nilpotent part has a decomposition into Jordan blocks whose lengths are given by the Jordan-partition calculation. In the implementation this is `JP_solve(s, t, 2)` in `src/PyPR/Tools/RootCounting/JordanPartition.py`.

Thus, for two `JordanSet` objects $A$ and $B$,

$$
A\cdot B:
\quad
c_b=\min(b,c_b^A+c_b^B),
\qquad
M_{A\cdot B}=\bigcup_{s\in M_A,\,t\in M_B}J_2(s,t),
$$

where $J_2(s,t)$ is the set of Jordan lengths returned by the characteristic-two Jordan partition. The truncation at $b$ records that coset weight cannot exceed the field degree (see [Coset Regions](../conventions/Coset%20Regions.md) for the geometric picture).

This explains why same-field multiplication cannot be represented by simply adding scalar multiplicities. For component fields of pairwise coprime degree, the eigenvalue choices are independent and the root counts multiply ([Proposition 3](Root%20Expressions%20and%20LC%20Estimation.md#the-five-root-expression-propositions)); for equal fields, generalized eigenspaces interact through tensor products and require the Jordan partition.

## Counting and the Upper Bound

The `upper()` method expands full coset classes (weight $= e$) into their embedded subfield classes, forms the corresponding counting regions (see [Coset Regions](../conventions/Coset%20Regions.md)), and uses the largest possible Jordan length in each term. This is intentionally optimistic: it assumes every allowed root/coset and every allowed generalized-root layer is present. The result is an upper bound on trajectory linear complexity, not a claim that all listed contributions survive in a concrete sequence.

The `lower()` method counts the same regions with a pessimistic weight for each coset class. It caps the class's weight at $e - 1$, which excludes the weight-$e$ root that the statistical model treats as absent with probability $1/2$ (see [Root Expressions — When the Upper Bound is Not Tight](Root%20Expressions%20and%20LC%20Estimation.md#when-the-upper-bound-is-not-tight)), and it multiplies the count by the factor $\left((2^e - 1)/2^e\right)^2$, so its result is in general not an integer. What that factor models is an open question recorded in [Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md) under the `estimate_LC` algorithm. The result is statistical rather than an algebraic lower bound for every initial state.

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

Throughout this section $E$ denotes the one-step shift on sequences, $(Es)[t] = s[t+1]$, so $(1 + E)s$ is the sequence of sums of consecutive terms, $s[t] \oplus s[t+1]$, and $(1 + E^j)s$ is the sequence $s[t] \oplus s[t+j]$. The letter $D$ is reserved for the variable of the D-transform, which enters in [The Resolvent Pushes Past the Ceiling](#the-resolvent-pushes-past-the-ceiling). In characteristic 2, the Frobenius endomorphism gives

$$(1 + E)^{2^k} = 1 + E^{2^k}.$$

This identity controls the period of sequences annihilated by powers of $(1 + E)$. If a sequence satisfies $(1 + E)^m s = 0$ with $m \geq 1$, let $k = \lceil \log_2 m \rceil$ be the smallest $k$ with $2^k \geq m$. Then $(1 + E)^{2^k} s = (1 + E^{2^k}) s = 0$, which says $s[t + 2^k] = s[t]$: the period divides $2^k$.

Equality holds when $m$ is exact, meaning $(1 + E)^{m-1} s \neq 0$. For $m = 1$ the sequence is a nonzero constant, with period $1 = 2^0$. For $m \geq 2$ we have $k \geq 1$ and $m > 2^{k-1}$, so $(1 + E^{2^{k-1}}) s = (1 + E)^{2^{k-1}} s \neq 0$ and the period does not divide $2^{k-1}$; a divisor of $2^k$ that does not divide $2^{k-1}$ is $2^k$ itself. So

$$\text{period} = 2^{\lceil \log_2 m \rceil}$$

where $m$ is the exact maximum Jordan block length for root $\alpha = 1$. When $m$ is only an upper bound on that length, as the Jordan sets of a root expression are, the conclusion is that the period divides $2^{\lceil \log_2 m \rceil}$.

Jordan lengths therefore do not produce arbitrary periods — they produce periods that jump in powers of 2. We call the range $[2^{k-1}+1,\; 2^k]$ the **$k$-th multiplicity bracket**: any exact Jordan length in this range yields the same period $2^k$.

### The Power-of-2 Ceiling on Products

Tensor products of Jordan blocks in characteristic 2 do not exceed the top of the bracket their factors lie under:

> **Theorem (Product Ceiling):** If $s \leq 2^k$ and $t \leq 2^k$, then every Jordan block in the decomposition of $J_s \otimes J_t$ (over a field of characteristic 2) has length $\leq 2^k$.

*Proof.* Let $A = J_s(\alpha) = \alpha I + N$ and $B = J_t(\beta) = \beta I + N'$, with $N^s = 0$ and $N'^t = 0$. Since $\alpha I$ and $N$ commute and $2^k \geq s$, the Frobenius identity gives

$$A^{2^k} = \alpha^{2^k} I + N^{2^k} = \alpha^{2^k} I,$$

and likewise $B^{2^k} = \beta^{2^k} I$. For $X = A \otimes B$ this gives $X^{2^k} = A^{2^k} \otimes B^{2^k} = (\alpha\beta)^{2^k} I$, and therefore

$$(X + \alpha\beta I)^{2^k} = X^{2^k} + (\alpha\beta)^{2^k} I = 0.$$

The only eigenvalue of $X$ is $\alpha\beta$, and $X + \alpha\beta I$ is nilpotent of index at most $2^k$, so every Jordan block of $X$ has length at most $2^k$. $\square$

The following instances are computed with `JP_solve(s, t, 2)`:

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

When both inputs are at most $2^k$, the product does not exceed the ceiling $2^k$. Once one input exceeds $2^k$, the product can reach $2^{k+1}$.

### The Resolvent Pushes Past the Ceiling

For a one-bit block whose chaining input shares its root, the resolvent adds exactly 1 to the maximum Jordan length. The argument is carried out on D-transforms, so the first step is to restate "annihilated by $(1+E)^m$" in transform terms.

Let $S(D) = \sum_{t \geq 0} s[t] D^t$ and $m \geq 1$. For $t \geq m$, the coefficient of $D^t$ in $(1+D)^m S(D)$ is

$$\sum_{j=0}^{m} \binom{m}{j} s[t - j] = \bigl((1+E)^m s\bigr)[t - m].$$

So $(1+E)^m s = 0$ if and only if $(1+D)^m S(D)$ is a polynomial $P$ of degree less than $m$, that is, $S(D) = P(D)/(1+D)^m$ with $\deg P < m$. Moreover $m$ is exact if and only if $P(1) \neq 0$. If $P(1) = 0$, then $P = (1+D)P'$ with $\deg P' < m - 1$, so $S = P'/(1+D)^{m-1}$ and $(1+E)^{m-1} s = 0$. Conversely, if $(1+E)^{m-1} s = 0$, then $S = P'/(1+D)^{m-1}$ with $\deg P' < m - 1$, so $P = (1+D)P'$ and $P(1) = 0$.

In the D-transform framework of [Resolvent Analysis](Resolvent%20Analysis.md), a block with update matrix $U$ and chaining input $\mathcal{C}$ has transform

$$B(D) = \frac{1}{\chi_U(D)} \operatorname{Adj}(I \oplus UD)\bigl(D\,\mathcal{C}(D) \oplus B[0]\bigr),$$

a one-bit block with polynomial $1 + x$ has $U = (1)$ and $\chi_U(D) = 1 + D$, so

$$B(D) = \frac{D\,\mathcal{C}(D) + B[0]}{1 + D}.$$

Suppose the chaining input has only the root 1, with exact maximum Jordan length $m \geq 1$. By the correspondence above, $\mathcal{C}(D) = P(D)/(1+D)^m$ for a polynomial $P$ with $\deg P < m$ and $P(1) \neq 0$. Then

$$B(D) = \frac{D\,P(D) + B[0]\,(1+D)^m}{(1+D)^{m+1}}.$$

The numerator has degree at most $m$, which is less than $m + 1$, and at $D = 1$ it equals $P(1) \neq 0$. By the correspondence again, the driven output is annihilated by $(1+E)^{m+1}$ exactly: its maximum Jordan length is $m + 1$.

Two limits on this statement. It applies when the chaining input shares the block's root; that is the hypothesis on $\mathcal{C}$ above. And it is specific to one-bit blocks: for a block with $\chi_U(D) = (1+D)^s$ and $s > 1$ the increment can be larger. With $U = J_2(1)$,

$$(I \oplus UD)^{-1} = \frac{1}{(1+D)^2}\begin{pmatrix} 1+D & D \\ 0 & 1+D \end{pmatrix},$$

which has the entry $D/(1+D)^2$ and so can add 2.

Unlike the product, which stays under the ceiling, the increment can take a Jordan length from $2^k$ to $2^k + 1$, crossing into bracket $k+1$ and doubling the period. The T-function analysis below uses exactly this contrast: products stay under the ceiling, and the one-bit resolvent step crosses it.

### T-Functions and Binomial Sequences

T-functions isolate the multiplicity structure, because they have no root structure other than multiplicities. Every block uses the polynomial $1 + x$, whose only root is $\alpha = 1$. The entire output is therefore a linear combination of **binomial sequences**:

$$u_r[t] = \binom{t}{r} \bmod 2, \qquad r = 0, 1, 2, \ldots$$

Each $u_r$ is generated by a Jordan block $J_{r+1}(1)$: $(1+E)^{r+1} u_r = 0$ but $(1+E)^r u_r \neq 0$. This follows from Pascal's rule, $\binom{t+1}{r} + \binom{t}{r} \equiv \binom{t}{r-1} \pmod 2$, which says $(1+E)u_r = u_{r-1}$ for $r \geq 1$, together with $(1+E)u_0 = 0$ for the constant sequence $u_0 = 1$. By the first subsection, the period of $u_r$ is $2^{\lceil \log_2(r+1) \rceil}$. Which binomial layers $u_r$ can be present in a T-function's output is what the Jordan length set `JordanSet.mults` records.

### Period Doubling Through Chaining

Consider an $n$-bit T-function where bit $i$ has update $b_i[t+1] = b_i[t] \oplus f_i(b_0[t], \ldots, b_{i-1}[t])$. The initial bit (bit 0) has no chaining inputs, but the T-function construction adds a constant driving term $\oplus\, 1$ — the "+1" of the binary counter. This constant has root $\alpha = 1$ with Jordan length 1, the root of the block's polynomial $1 + x$, so the resolvent step raises the Jordan length from 1 to 2: bit 0 alternates with period 2, not period 1. Its Jordan set is $\{1, 2\}$.

Each successive bit's Jordan set is then determined by:

1. **Chaining:** The function $f_i$ combines earlier bits' Jordan sets via products (AND) and unions (XOR). By the Product Ceiling theorem, these products cannot exceed the ceiling their inputs lie under.
2. **Resolvent:** The one-bit block adds 1 to the maximum, which can push it past the ceiling.

In terms of these sets, the period bound doubles at bit $i$ exactly when the chaining maximum reaches the ceiling $2^i$, so that the resolvent increment $2^i \to 2^i + 1$ crosses into the next bracket. With full nonlinear mixing, the Jordan set evolution is:

| Bit | Chaining max | Resolvent max | Bracket | Period |
|---|---|---|---|---|
| 0 | 1 (constant) | 2 | $2^1$ | 2 |
| 1 | 2 | 3 | $2^2$ | 4 |
| 2 | 4 | 5 | $2^3$ | 8 |
| 3 | 8 | 9 | $2^4$ | 16 |
| $i$ | $2^i$ | $2^i+1$ | $2^{i+1}$ | $2^{i+1}$ |

If the chaining at bit $i$ fails to reach the ceiling, the resolvent stays within the current bracket and the period bound does not increase.

The Jordan length and period in the last row are exact for a single-cycle T-function — one whose state sequence has the full period $2^n$ — and can be derived without the Jordan sets. In a single-cycle T-function the low $j$ bits have period $2^j$, for every $j \leq n$. The reason is that the update is triangular: bit $l$ depends only on bits $0, \ldots, l$, so the low $j$ bits evolve by their own map on $2^j$ states. A single cycle on all $n$ bits visits every state, so its projection visits all $2^j$ states of the low bits, and an orbit of a map that visits all $2^j$ of its states is a single cycle of length $2^j$. Apply this with $j = i$ and $j = i + 1$. After $2^i$ steps the low $i$ bits have returned to their values and the low $i + 1$ bits have not, so bit $i$ has flipped:

$$b_i[t + 2^i] = b_i[t] \oplus 1 \quad \text{for all } t.$$

Hence $(1+E)^{2^i} b_i = (1 + E^{2^i}) b_i$ is the constant sequence 1, which is nonzero, and one more factor gives $(1+E)^{2^i + 1} b_i = 0$. So bit $i$ has exact Jordan length $2^i + 1$, which lies in bracket $i + 1$, and period $2^{i+1}$. For the last bit of an $n$-bit T-function ($i = n - 1$), the Jordan length is $2^{n-1} + 1$, the sequence is annihilated by $(1+E)^{2^{n-1}+1}$, and the period is $2^n$.

### Why Every Prior Bit Is Required

This subsection concerns the binary counter $x \mapsto x + 1$, and its figures are exact linear complexities of the counter's sequences, not Jordan-set bounds. All the sequences involved have only the root 1, so the linear complexity of each is its exact Jordan length $m$, read from its D-transform by the correspondence of [The Resolvent Pushes Past the Ceiling](#the-resolvent-pushes-past-the-ceiling).

For the counter, the chaining at bit $i$ is the product of all bits $0, \ldots, i-1$. That product is 1 exactly when the low $i$ bits are all 1, so it is the indicator sequence of a single residue class $t \equiv c \pmod{2^i}$, with $0 \leq c < 2^i$. Its transform is

$$\sum_{l \geq 0} D^{c + l 2^i} = \frac{D^c}{1 + D^{2^i}} = \frac{D^c}{(1+D)^{2^i}}.$$

The numerator is nonzero at $D = 1$ and has degree $c < 2^i$, so the linear complexity is $2^i$, the ceiling.

Without bit $j$, the product is 1 when the other low bits are all 1 and bit $j$ is arbitrary. That is the indicator of two residue classes, $c$ and $c + 2^j$, with $0 \leq c < 2^i - 2^j$. Its transform is

$$\frac{D^c\,(1 + D^{2^j})}{(1+D)^{2^i}} = \frac{D^c\,(1+D)^{2^j}}{(1+D)^{2^i}} = \frac{D^c}{(1+D)^{2^i - 2^j}},$$

with numerator nonzero at $D = 1$ and of degree $c < 2^i - 2^j$. So the linear complexity is $2^i - 2^j$:

| Dropped bit $j$ | LC of bit $j$ | LC of the product without bit $j$ | Shortfall from $2^i$ |
|---|---|---|---|
| 0 | 2 | $2^i - 1$ | $2^0 = 1$ |
| 1 | 3 | $2^i - 2$ | $2^1 = 2$ |
| 2 | 5 | $2^i - 4$ | $2^2 = 4$ |
| $j$ | $2^j + 1$ | $2^i - 2^j$ | $2^j$ |

In every case $2^i - 2^j < 2^i$. The resolvent step then gives bit $i$ exact Jordan length $2^i - 2^j + 1 \leq 2^i$, so its period divides $2^i$ and does not double.

This recovers, for the counter, the necessity direction of the T-function full-periodicity condition: **if the chaining at bit $i$ omits a prior bit, the period does not double.** The argument is about the counter with one bit removed from its carry product. It does not address sufficiency, and it is not established for general T-functions.

> **OPEN:** Which T-function construction makes full periodicity equivalent to each $f_i$ depending nontrivially on $b_{i-1}$? — needs author input.
> An earlier version of this section stated that equivalence for "the standard T-function construction". It contradicts the table above: the chaining $f_2 = b_1$ depends on $b_1$ but omits $b_0$, and the linear complexity of the product is then $2^2 - 2^0 = 3 < 4$.

## Practical Consequences

Repeated component sizes are the main reason the [mesh optimization](../architecture/Mesh%20Optimization.md) is disabled in some CMPR configurations: equal-sized components share field structure, so their products need multiplicity-aware Jordan handling. T-Functions are especially exposed because every component uses the same polynomial $1 + x$, making the bracket/staircase structure of [Multiplicity Brackets and the Period Staircase](#multiplicity-brackets-and-the-period-staircase) the dominant factor in period and linear complexity analysis. For those cases, use the default root-expression path. The tightness of the resulting bounds in these configurations has not been established, and they should be checked against measured values.
