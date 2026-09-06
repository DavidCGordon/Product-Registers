# Möbius Inversion on the Bit-Subset Lattice

This document develops the Möbius transform on the Boolean lattice — the lattice-theoretic tool that converts between a Boolean function's truth table and its [Algebraic Normal Form](Algebraic%20Normal%20Form.md). Matrix rows, columns, and vector orientation follow [Matrix Indexing](../conventions/Matrix%20Indexing.md). The approach here proves the transform works via interval counting on the bit-subset poset; for the equivalent matrix perspective (Kronecker products, Pascal's triangle mod 2, butterfly decomposition), see the [ANF doc](Algebraic%20Normal%20Form.md#the-transform-matrix). The connection to binomial coefficients mod 2, via Lucas' theorem, makes this the natural framework for decomposing expressions like $\binom{dn}{k} \pmod{2}$ (see [Binomial Decimation Identity](Binomial%20Decimation%20Identity.md)).

## The Bit-Subset Partial Order

**Definition.** For non-negative integers $j$ and $i$, write $j \preccurlyeq i$ (read "$j$ is a bit-subset of $i$") when every bit that is set in $j$ is also set in $i$. Equivalently, $j \mathop{\&} i = j$, where $\&$ is bitwise AND.

**Examples.**
- $5 = 101_2 \preccurlyeq 7 = 111_2$, because bits 0 and 2 of 5 are both set in 7.
- $5 = 101_2 \not\preccurlyeq 6 = 110_2$, because bit 0 of 5 is not set in 6.
- $0 \preccurlyeq i$ for every $i$ (the empty bit-set is a subset of any bit-set).
- $i \preccurlyeq i$ for every $i$ (every bit-set is a subset of itself).

This is a partial order: it is reflexive, antisymmetric, and transitive. When restricted to integers with at most $m$ bits (i.e., $\{0, \ldots, 2^m - 1\}$), it is isomorphic to the inclusion order on subsets of $\{0, \ldots, m-1\}$, via the identification $i \leftrightarrow \text{bits}(i)$.

## Connection to Binomial Coefficients

Lucas' theorem (for the prime $p = 2$) states that a binomial coefficient mod 2 factors into the product of its binary digits:

$$\binom{n}{i} \;\equiv\; \prod_{b \geq 0} \binom{n_b}{i_b} \pmod{2}$$

where $n_b$ and $i_b$ are the $b$-th binary digits of $n$ and $i$. Each factor $\binom{n_b}{i_b}$ is a binomial coefficient of single bits (0 or 1), so the only way a factor can be zero is $\binom{0}{1} = 0$ — that is, when $i_b = 1$ but $n_b = 0$. If no such bit exists, every factor is 1 and the product is 1.

This gives the equivalence:

$$\binom{n}{i} \equiv 1 \pmod{2} \quad\Longleftrightarrow\quad i \preccurlyeq n$$

In other words, $\binom{n}{i} \pmod{2}$ is precisely the indicator function of the bit-subset relation. We will use binomial coefficients and the $\preccurlyeq$ relation interchangeably throughout: $\binom{n}{i} \pmod{2}$ is the analytic form, $i \preccurlyeq n$ is the combinatorial form, but they compute the same thing.

## Algebraic Normal Form

A **Boolean function** on $m$ bits is a function $f: \{0,1\}^m \to \{0,1\}$. We encode inputs by identifying the bit-vector $(n_0, \ldots, n_{m-1})$ with the integer $n = \sum_b n_b \, 2^b$, so $f$ becomes a function from $\{0, \ldots, 2^m - 1\}$ to $\{0,1\}$.

**Theorem (ANF existence and uniqueness).** Every Boolean function $f$ on $m$ bits has a unique representation

$$f(n) \;=\; \bigoplus_{i=0}^{2^m - 1}\; a_i \,\binom{n}{i} \pmod{2}$$

with coefficients $a_i \in \{0,1\}$. This is the **Algebraic Normal Form** (ANF) of $f$.

Each term $\binom{n}{i} \pmod{2}$ acts as a monomial: by Lucas' theorem and the bit-product identity above, $\binom{n}{i} \equiv \prod_{b \in \text{bits}(i)} n_b \pmod{2}$. When $i = 0$, the empty product is 1 (the constant monomial). When $i = 3 = 11_2$, the monomial is $n_0 \cdot n_1$. And so on — each $i$ selects a distinct subset of input bits to AND together, and the ANF is the XOR of the selected monomials.

**Proof of uniqueness.** There are $2^{2^m}$ Boolean functions and $2^{2^m}$ possible coefficient vectors $(a_0, \ldots, a_{2^m - 1})$, so it suffices to show the representation map is injective — equivalently, that the $2^m$ functions $\{n \mapsto \binom{n}{i}\}_{i=0}^{2^m - 1}$ are linearly independent over $\mathbb{F}_2$.

Consider the evaluation matrix $M$ with rows indexed by $n$ and columns by $i$, where $M_{n,i} = \binom{n}{i} \pmod{2}$. Order both indices so that $\preccurlyeq$ refines the ordering (any linear extension of $\preccurlyeq$ works — for instance, ordinary $\leq$). Then $M_{n,i} = 1$ requires $i \preccurlyeq n$, so $M$ is lower-triangular with 1s on the diagonal ($\binom{i}{i} = 1$). A triangular matrix with 1s on the diagonal is invertible over any field, so the columns are linearly independent. $\square$

**Remark.** The lower-triangular structure is the reason Möbius inversion works: the matrix $M$ and its inverse have the same sparsity pattern (both supported on $\{(n,i) : i \preccurlyeq n\}$), and over $\mathbb{F}_2$ the inverse turns out to be $M$ itself (see below).

## The Möbius Transform

The ANF representation says $f(n) = \bigoplus_{i \preccurlyeq n} a_i$: to evaluate $f$ at $n$, XOR together all $a_i$ for which $i \preccurlyeq n$. This is the **zeta transform** (or **downward sum**) of the coefficient sequence $a$, evaluated at $n$.

To go the other direction — recovering the ANF coefficients $a_i$ from the truth table $f$ — apply the same operation:

$$\boxed{a_i \;=\; \bigoplus_{j \,\preccurlyeq\, i}\; f(j)}$$

This is the **Möbius transform** on the bit-subset lattice.

### Proof of Correctness

We verify that the formula recovers $a_i$ by substituting the definition of $f$ and showing everything cancels except the $a_i$ term.

Start from the right-hand side and expand each $f(j)$ using the ANF:

$$\bigoplus_{j \,\preccurlyeq\, i}\; f(j) \;=\; \bigoplus_{j \,\preccurlyeq\, i}\; \bigoplus_{l \,\preccurlyeq\, j}\; a_l$$

Each coefficient $a_l$ appears in this double sum once for every $j$ satisfying $l \preccurlyeq j \preccurlyeq i$ — that is, once for every element of the **interval** between $l$ and $i$ in the bit-subset order. Over $\mathbb{F}_2$, repeated terms cancel in pairs, so $a_l$ survives if and only if the interval has an odd number of elements.

**Counting the interval.** If $l \not\preccurlyeq i$, no $j$ can simultaneously satisfy $l \preccurlyeq j$ and $j \preccurlyeq i$, so the interval is empty (size 0, which is even — $a_l$ cancels).

If $l \preccurlyeq i$, then any $j$ in the interval must have all bits of $l$ set (because $l \preccurlyeq j$) and no bits outside $i$ set (because $j \preccurlyeq i$). The remaining bits — those set in $i$ but not in $l$ — may be independently set or cleared in $j$. There are $\text{popcount}(i) - \text{popcount}(l)$ such free bits, giving

$$|\{j : l \preccurlyeq j \preccurlyeq i\}| \;=\; 2^{\,\text{popcount}(i) \,-\, \text{popcount}(l)}$$

This count is odd only when the exponent is 0, which happens only when $l = i$. For every other $l \preccurlyeq i$ (where $l \neq i$), the exponent is at least 1, the count is even, and $a_l$ cancels.

Therefore:

$$\bigoplus_{j \,\preccurlyeq\, i}\; f(j) \;=\; a_i \qquad\square$$

### The Transform is Self-Inverse

Define the operator $\mathcal{M}$ by $(\mathcal{M}g)(i) = \bigoplus_{j \preccurlyeq i} g(j)$. The ANF representation says $f = \mathcal{M}a$ (the truth table is the Möbius transform of the coefficients), and the inversion formula says $a = \mathcal{M}f$ (the coefficients are the Möbius transform of the truth table). Together: $\mathcal{M}(\mathcal{M}g) = g$ for any function $g$. The Möbius transform is its own inverse — an **involution**.

This self-inverse property is special to $\mathbb{F}_2$. Over $\mathbb{Z}$, the zeta transform (sum over predecessors) and Möbius transform (alternating sum) are distinct operations that require different sign patterns to invert each other. Over $\mathbb{F}_2$, addition and subtraction are the same, so the transform and its inverse coincide.

## Example

Consider Boolean functions on 2 bits, so $n \in \{0, 1, 2, 3\}$.

Define $f$ by its truth table: $f(0) = 0,\; f(1) = 1,\; f(2) = 0,\; f(3) = 0$. This function outputs 1 only when $n = 1$ — bit 0 is set and bit 1 is not.

Apply the Möbius transform to recover the ANF coefficients:

| $i$ | bit-subsets $j \preccurlyeq i$ | $\bigoplus_j f(j)$ | $a_i$ |
|-----|-------------------------------|---------------------|-------|
| $0 = 00$ | $\{0\}$ | $f(0) = 0$ | 0 |
| $1 = 01$ | $\{0, 1\}$ | $f(0) \oplus f(1) = 0 \oplus 1$ | 1 |
| $2 = 10$ | $\{0, 2\}$ | $f(0) \oplus f(2) = 0 \oplus 0$ | 0 |
| $3 = 11$ | $\{0, 1, 2, 3\}$ | $f(0) \oplus f(1) \oplus f(2) \oplus f(3) = 0 \oplus 1 \oplus 0 \oplus 0$ | 1 |

The ANF is $f(n) = \binom{n}{1} + \binom{n}{3} \pmod{2}$, or equivalently $f(n) = n_0 + n_0 n_1 \pmod{2}$.

Verification:

| $n$ | $\binom{n}{1}$ | $\binom{n}{3}$ | sum | $f(n)$ |
|-----|----------------|----------------|-----|--------|
| 0 | 0 | 0 | 0 | 0 ✓ |
| 1 | 1 | 0 | 1 | 1 ✓ |
| 2 | 0 | 0 | 0 | 0 ✓ |
| 3 | 1 | 1 | 0 | 0 ✓ |

Note that $\binom{3}{3} = 1$ contributes at $n = 3$, which is exactly the cancellation that makes $f(3) = 0$: the monomial $n_0$ alone would give 1 at $n = 3$, but the correction term $n_0 n_1$ flips it back.

## Nomenclature and Relation to Other Transforms

### The Classical Möbius Function (Number Theory)

In number theory, **Möbius inversion** is a technique for recovering a function from its cumulative sums over divisors. It uses two ingredients:

The **divisor-sum relation.** Many arithmetic functions arise as sums over divisors. If $g$ is some function on positive integers, define $f(n) = \sum_{d \mid n} g(d)$ — the sum of $g$ over all divisors of $n$. For example, the sum-of-divisors function $\sigma(n) = \sum_{d \mid n} d$ has this form with $g(d) = d$.

The **Möbius function** $\mu(n)$ inverts divisor sums. It is defined by:

$$\mu(n) = \begin{cases} 1 & \text{if } n = 1 \\ (-1)^k & \text{if } n = p_1 p_2 \cdots p_k \text{ for distinct primes } p_i \\ 0 & \text{if } n \text{ has any squared prime factor} \end{cases}$$

The **inversion formula** recovers $g$ from $f$: if $f(n) = \sum_{d \mid n} g(d)$, then $g(n) = \sum_{d \mid n} \mu(n/d)\, f(d)$.

This works because $\mu$ satisfies the fundamental identity $\sum_{d \mid n} \mu(d) = [n = 1]$ — the Möbius function cancels all terms except at $n = 1$, which is the same cancellation mechanism behind the Boolean lattice proof in this document.

### Rota's Generalization: Möbius Inversion on Posets

Gian-Carlo Rota observed that the divisor-sum / Möbius-inversion pattern is not special to divisibility — it works on any **partially ordered set** (poset): a set $P$ with a relation $\leq$ that is reflexive ($x \leq x$), antisymmetric ($x \leq y$ and $y \leq x$ implies $x = y$), and transitive ($x \leq y$ and $y \leq z$ implies $x \leq z$). The poset is **locally finite** if every interval $[x, y] = \{z : x \leq z \leq y\}$ is a finite set.

On any locally finite poset, define:

- The **zeta function**: $\zeta(x, y) = 1$ if $x \leq y$, and $0$ otherwise. This is the indicator of the partial order — the analog of "$d$ divides $n$."

- The **Möbius function**: $\mu(x, y)$ is defined recursively for $x \leq y$ by $\mu(x, x) = 1$ and

$$\mu(x, y) = -\sum_{x \leq z < y} \mu(x, z) \qquad\text{for } x < y$$

This ensures $\sum_{x \leq z \leq y} \mu(x, z) = [x = y]$ — the Möbius function sums to zero over every non-trivial interval, providing the cancellation needed for inversion.

The **generalized inversion formula**: if $f(x) = \sum_{y \leq x} g(y)$ (sum over all predecessors of $x$), then $g(x) = \sum_{y \leq x} \mu(y, x)\, f(y)$.

The number-theoretic case is recovered by taking $P = \mathbb{Z}_{>0}$ with $\leq$ being divisibility ($d \mid n$). For this poset, $\mu(1, n)$ gives exactly the classical $\mu(n)$.

### Specialization to the Boolean Lattice

Our setting is the **Boolean lattice**: the set $\{0, \ldots, 2^m - 1\}$ ordered by the bit-subset relation $\preccurlyeq$. This is a locally finite poset (every interval is finite), so Rota's theory applies.

The interval $[l, i] = \{j : l \preccurlyeq j \preccurlyeq i\}$ is isomorphic to a Boolean lattice on $\text{popcount}(i) - \text{popcount}(l)$ bits (the "free" bits from the [proof of correctness](#proof-of-correctness)). The Möbius function on a Boolean lattice of $k$ bits evaluates to $(-1)^k$ at the endpoints. So:

$$\mu(l, i) = (-1)^{\text{popcount}(i) - \text{popcount}(l)}$$

Over $\mathbb{F}_2$, every $(-1)^k$ reduces to 1, so $\mu \equiv 1$ on all intervals. This is why the inversion formula over $\mathbb{F}_2$ uses the same unsigned sum as the zeta transform:

$$a_i = \sum_{j \preccurlyeq i} \mu(j, i)\, f(j) = \sum_{j \preccurlyeq i} 1 \cdot f(j) = \bigoplus_{j \preccurlyeq i} f(j)$$

The zeta and Möbius transforms coincide — the transform is self-inverse — because the Möbius function is constantly 1 mod 2. Over $\mathbb{Z}$, the signs would alternate and the two transforms would differ: the zeta direction would sum, and the Möbius direction would alternate (inclusion-exclusion).

### The Walsh-Hadamard Transform: A Different Decomposition

The **Walsh-Hadamard transform** (WHT) is a separate transform that also operates on Boolean-indexed vectors but uses a different base matrix, works over $\mathbb{R}$ instead of $\mathbb{F}_2$, and measures something fundamentally different.

**The $\{0,1\} \to \{-1,+1\}$ lifting.** The WHT works by reinterpreting a Boolean function $f: \{0,1\}^m \to \{0,1\}$ as a real-valued function $\hat{f}: \{0,1\}^m \to \{-1,+1\}$ via the substitution $b \mapsto (-1)^b$: output 0 becomes $+1$ and output 1 becomes $-1$. This lifts $f$ from the additive world of $\mathbb{F}_2$ (where XOR is addition) to the multiplicative world of $\mathbb{R}$ (where XOR becomes multiplication of $\pm 1$ signs).

**Linear functions and characters.** For each subset $S \subseteq \{0, \ldots, m-1\}$, define the **linear function** $\ell_S(x) = \bigoplus_{b \in S} x_b$ — the XOR (parity) of the input bits in $S$. After lifting, this becomes the **character** $\chi_S(x) = (-1)^{\ell_S(x)} = \prod_{b \in S} (-1)^{x_b}$, which outputs $+1$ when the parity is even and $-1$ when odd.

**Walsh coefficients as correlations.** The Walsh-Hadamard transform computes, for each subset $S$, how correlated $f$ is with the linear function $\ell_S$:

$$W_f(S) = \sum_{x \in \{0,1\}^m} (-1)^{f(x) \oplus \ell_S(x)}$$

Each term is $+1$ when $f(x)$ agrees with $\ell_S(x)$ and $-1$ when they disagree, so $W_f(S)$ counts agreements minus disagreements. A large positive $W_f(S)$ means $f$ closely resembles $\ell_S$; a large negative value means $f$ closely resembles its complement; $W_f(S) = 0$ means $f$ and $\ell_S$ agree on exactly half the inputs.

**The base matrix.** For a single bit ($m = 1$), the WHT matrix is the **Hadamard matrix**

$$H = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$$

Applied to the vector $(f(0), f(1))$, it computes:
- $f(0) + f(1)$: the sum — measures the "constant" component (correlation with $\ell_\emptyset$, which is always 0).
- $f(0) - f(1)$: the difference — measures correlation with the input bit $x_0$.

The $m$-bit WHT is the Kronecker power $H^{\otimes m}$, just as the Möbius transform is $T_1^{\otimes m}$.

**Comparison.** Both transforms decompose Boolean functions into a basis, but the bases are fundamentally different:

| | Möbius / ANF | Walsh-Hadamard |
|---|---|---|
| **Base matrix** | $T_1 = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$ | $H = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$ |
| **Operates over** | $\mathbb{F}_2$ | $\mathbb{R}$ |
| **Basis functions** | AND-monomials $\prod_{b \in S} x_b$ | parity characters $(-1)^{\bigoplus_{b \in S} x_b}$ |
| **Coefficients tell you** | which monomials appear in the polynomial | how correlated $f$ is with each linear function |
| **Used for** | algebraic degree, ANF structure | nonlinearity, correlation immunity, Walsh spectrum |
| **Self-inverse?** | yes (over $\mathbb{F}_2$) | yes (up to scaling by $2^m$) |

Over $\mathbb{F}_2$, the Hadamard matrix collapses to $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ (since $-1 \equiv 1$), which is singular — the WHT does not exist over $\mathbb{F}_2$. The Möbius transform replaces it entirely for the $\mathbb{F}_2$ setting.

The two transforms answer different questions about the same function: the ANF/Möbius decomposition reveals the algebraic structure (which products of variables appear), while the Walsh-Hadamard decomposition reveals the statistical structure (how far $f$ is from each linear function). Both are computed via butterfly operations from their Kronecker factorizations, which is the main source of the naming confusion — some sources call any such Kronecker-butterfly transform "Walsh-Hadamard." In this project, we use **Möbius transform** exclusively for the $T_1^{\otimes m}$ operation over $\mathbb{F}_2$.

## Connections

- **[Algebraic Normal Form](Algebraic%20Normal%20Form.md):** Full treatment of ANF — the standard polynomial definition, the binomial basis, and the matrix/Kronecker perspective on the same transform proved here.
- **[Binomial Decimation Identity](Binomial%20Decimation%20Identity.md):** Uses this transform to decompose $\binom{dn}{k} \pmod{2}$ into the binomial basis.
- **[Algebraic Attacks](Algebraic%20Attacks.md):** The ANF representation is the foundation of algebraic cryptanalysis — annihilators, low-degree relations, and equation generation all operate in the ANF ring.
