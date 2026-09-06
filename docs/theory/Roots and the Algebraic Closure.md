# Roots and the Algebraic Closure

This document develops the algebraic setting in which roots of feedback register sequences live: the algebraic closure $\overline{\mathbb{F}_2}$, its subfield lattice, the exponential representation of roots, and the key polynomial types (minimal, primitive, characteristic) that connect field elements to polynomials. These ideas underpin the [primal/dual framework](../conventions/Polynomial%20Conventions.md), [root expression estimation](Root%20Expressions%20and%20LC%20Estimation.md), [root multiplicities and Jordan decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md), and [resolvent analysis](Resolvent%20Analysis.md).

For notation conventions (coefficient ordering, block indexing, vocabulary), see [Notation and Terminology](../conventions/Notation%20and%20Terminology.md).

## The Algebraic Closure of $\mathbb{F}_2$

The base field $\mathbb{F}_2 = \{0, 1\}$ is too small to contain the roots of most polynomials over it. For instance, $x^2 + x + 1$ is irreducible over $\mathbb{F}_2$ — it has no roots in $\{0, 1\}$. To factor polynomials completely and express sequences in terms of their spectral components, we need a field large enough to contain every root of every polynomial over $\mathbb{F}_2$.

The **algebraic closure** $\overline{\mathbb{F}_2}$ is that field. It is the smallest field extension of $\mathbb{F}_2$ in which every polynomial $P \in \mathbb{F}_2[x]$ splits into linear factors. It is unique up to isomorphism, infinite, and has characteristic 2 (so $1 + 1 = 0$, and addition is XOR).

Concretely, $\overline{\mathbb{F}_2}$ is the union of all finite fields of characteristic 2:

$$\overline{\mathbb{F}_2} = \bigcup_{n \geq 1} \mathbb{F}_{2^n}$$

This union is well-defined because the finite fields nest in a lattice governed by divisibility: $\mathbb{F}_{2^m} \subseteq \mathbb{F}_{2^n}$ if and only if $m \mid n$. The lattice means that, for example, $\mathbb{F}_{2^3}$ and $\mathbb{F}_{2^6}$ share a common subfield $\mathbb{F}_{2^3}$, but $\mathbb{F}_{2^3}$ and $\mathbb{F}_{2^5}$ meet only at $\mathbb{F}_2$ itself (since $\gcd(3, 5) = 1$).

Every element of $\overline{\mathbb{F}_2}$ lives in some specific finite subfield $\mathbb{F}_{2^n}$, and the smallest such $n$ determines the element's algebraic degree over $\mathbb{F}_2$. This is why the subfield structure matters: it tells us how "complex" each root is and which roots can interact algebraically.

### Why the Algebraic Closure Matters for Registers

A feedback register of width $n$ bits has its state space naturally identified with $\mathbb{F}_{2^n}$ (for linear registers) or $\mathbb{F}_2^n$ (as a vector space). The output sequence of a linear register satisfies a linear recurrence, and the roots of that recurrence's characteristic polynomial live in $\overline{\mathbb{F}_2}$ — typically spread across several subfields. The [root expression](Root%20Expressions%20and%20LC%20Estimation.md) framework tracks which subfields contain active roots, and the [primal/dual distinction](../conventions/Polynomial%20Conventions.md) determines which polynomial's roots we mean.

## The Exponential Representation

Elements of $\overline{\mathbb{F}_2}$ can be named using a universal exponential notation that makes their period and field membership immediately visible.

### The Notation $\alpha^{p/q}$

Fix a compatible system of primitive elements (described in more rigor below): for each $n$, let $\alpha_n$ be a primitive element of $\mathbb{F}_{2^n}$ (a generator of $\mathbb{F}_{2^n}^\times$), chosen so that the $\alpha_n$ are compatible across subfield embeddings. Then every nonzero element of $\overline{\mathbb{F}_2}$ can be written as $\alpha^{p/q}$ where:

- $q$ is a positive odd integer (the **period** of the element under repeated squaring in an appropriate sense),
- $p$ is an integer with $\gcd(p, q) = 1$ and $0 < p < q$,
- the fraction $p/q$ is interpreted modulo 1 in the group $\mathbb{Q}/\mathbb{Z}$ restricted to 2-adic rationals.

More concretely, $\alpha^{p/q}$ denotes the element of multiplicative order $q$ in the unique subfield that contains it. Since $\mathbb{F}_{2^n}^\times$ is cyclic of order $2^n - 1$, an element of multiplicative order $q$ exists in $\mathbb{F}_{2^n}$ precisely when $q \mid 2^n - 1$.

The exponential notation makes a wealth of structural properties immediately visible — period, field membership, Frobenius conjugates, multiplicative arithmetic — developed in full in the [Properties of Roots](#properties-of-roots) section below. First, we establish the infrastructure for working with these elements concretely: the compatible hierarchy that makes the notation well-defined, the analogy with $\mathbb{C}$ that frames the interplay between representations, and the concrete polynomial/matrix/vector forms the library uses.

## The Compatible Primitive Element Hierarchy

The exponential notation $\alpha^{p/q}$ presupposes a "compatible system of primitive elements" — one primitive element $\alpha_n$ per field $\mathbb{F}_{2^n}$, chosen so that they agree across subfield embeddings. This section makes that precise, proves it exists, and explains why it matters.

### The Compatibility Condition

When $m \mid n$, the subfield $\mathbb{F}_{2^m}$ embeds into $\mathbb{F}_{2^n}$. A primitive element $\alpha_n$ of $\mathbb{F}_{2^n}$ does not automatically restrict to a primitive element of $\mathbb{F}_{2^m}$ — the norm map $\alpha_n \mapsto \alpha_n^{(2^n - 1)/(2^m - 1)}$ sends $\alpha_n$ to an element of $\mathbb{F}_{2^m}^\times$, but that element might not be a generator.

A system $\{\alpha_n\}_{n \geq 1}$ is **compatible** if for every pair $m \mid n$:

$$\alpha_n^{(2^n - 1)/(2^m - 1)} = \alpha_m$$

This ensures that the exponential label of an element does not depend on which ambient field we view it in: the element $\alpha_m^j$ in $\mathbb{F}_{2^m}$ has the same identity whether we name it using $\alpha_m$ or by embedding into $\mathbb{F}_{2^n}$ and using $\alpha_n$.

### Existence

**Claim.** A compatible system of primitive elements exists.

*Proof sketch.* The system is built by an inverse limit construction over the divisibility poset of positive integers. For each $n$, we need $\alpha_n$ primitive in $\mathbb{F}_{2^n}$ with the compatibility constraint. The key observation is that the norm map $N_{n/m}: \mathbb{F}_{2^n}^\times \to \mathbb{F}_{2^m}^\times$ given by $x \mapsto x^{(2^n - 1)/(2^m - 1)}$ is a surjective group homomorphism (surjectivity follows because $\mathbb{F}_{2^m}^\times$ is cyclic and the image contains an element of order $2^m - 1$). Since each $N_{n/m}$ is surjective, the inverse system $(\mathbb{F}_{2^n}^\times, N_{n/m})$ satisfies the Mittag-Leffler condition (all transition maps are surjective between finite groups), so the inverse limit $\varprojlim \mathbb{F}_{2^n}^\times$ is nonempty. Any element of the inverse limit gives a compatible system.

More concretely: start with any primitive $\alpha_1 = 1$ in $\mathbb{F}_2^\times$. For each successive $n$ (in an order compatible with divisibility), choose $\alpha_n$ primitive in $\mathbb{F}_{2^n}$ such that $\alpha_n^{(2^n-1)/(2^m-1)} = \alpha_m$ for all $m \mid n$ with $m < n$. The surjectivity of $N_{n/m}$ guarantees that among the $\phi(2^n - 1)$ primitive elements, at least one satisfies all the finitely many constraints from smaller fields.

The system is far from unique — different choices of compatible hierarchy give different "coordinate systems" on $\overline{\mathbb{F}_2}$, but the algebraic structure is the same in each. $\square$

### Why the Hierarchy Matters

The compatible hierarchy provides a canonical bridge between the exponential and polynomial representations.

**Naming across fields.** Without compatibility, saying "$\alpha^{3/7}$" is ambiguous: it depends on which primitive element we chose for $\mathbb{F}_{2^3}$. With a compatible hierarchy, the name is universal — the element $\alpha^{3/7}$ is the same whether we work in $\mathbb{F}_{2^3}$, $\mathbb{F}_{2^6}$, $\mathbb{F}_{2^{12}}$, or any other field containing $\mathbb{F}_{2^3}$.

**Enabling addition.** The exponential representation has no direct formula for addition: $\alpha^{p_1/q_1} + \alpha^{p_2/q_2}$ cannot in general be simplified to $\alpha^{p_3/q_3}$ without computing in a concrete representation. To add two elements, one must:

1. Identify the smallest common field $\mathbb{F}_{2^N}$ containing both (where $N = \text{lcm}(\text{Ord}(q_1), \text{Ord}(q_2))$; see the LCM rule in Properties below).
2. Choose a representation — a primitive polynomial $P_N$ of degree $N$ — which fixes $\alpha_N$ as the root of $P_N$.
3. Express both elements as polynomials of degree $< N$ in $\alpha_N$ (their vector representations).
4. Add the vectors (XOR the coefficient arrays).
5. If needed, convert the result back to exponential form (a discrete logarithm).

Steps 2–4 are polynomial arithmetic modulo $P_N$. The compatible hierarchy ensures that the embedding in step 3 is consistent: $\alpha^{p_1/q_1}$ is $\alpha_N^{p_1(2^N - 1)/q_1}$, which can be computed as a polynomial in $\alpha_N$ by repeated squaring modulo $P_N$.

This asymmetry — multiplication is trivial in exponential form, addition requires a detour through polynomial form — is a fundamental feature of the algebraic closure, not a limitation of the notation. It mirrors the situation in other domains (see the analogy with $\mathbb{C}$ below).

## The Two Faces of $\overline{\mathbb{F}_2}$: Analogy with $\mathbb{C}$

The interplay between exponential and polynomial representations in $\overline{\mathbb{F}_2}$ closely parallels the interplay between polar and rectangular coordinates in $\mathbb{C}$. The analogy is deeper than it first appears, and making it precise clarifies both what the compatible hierarchy buys us and when each representation is the right tool.

### The Analogy

| | $\mathbb{C}$ | $\overline{\mathbb{F}_2}$ |
|---|---|---|
| **"Polar" form** | $r e^{i\theta}$ | $\alpha^{p/q}$ |
| **"Rectangular" form** | $a + bi$ | $b_0 + b_1 \alpha + \cdots + b_{n-1}\alpha^{n-1}$ |
| **Multiplication** | Easy in polar: multiply magnitudes, add angles | Easy in exponential: add exponents |
| **Addition** | Easy in rectangular: add components | Easy in polynomial: XOR coefficient vectors |
| **Converting polar → rectangular** | $r\cos\theta + i \cdot r\sin\theta$ | Repeated squaring of $\alpha_N$ modulo $P_N$ |
| **Converting rectangular → polar** | $r = \|z\|$, $\theta = \arg(z)$ | Discrete logarithm (expensive) |
| **"Basis" choice** | $i$ vs $-i$ (two options) | The compatible hierarchy $\{\alpha_n\}$ (infinitely many) |

### The Choice of Roots

The last row deserves elaboration. Over $\mathbb{C}$, there is a genuine choice: either of the two roots of $x^2 + 1$ could serve as "$i$." Calling one $i$ and the other $-i$ is a convention — the algebraic structure of $\mathbb{C}$ is invariant under the swap $i \leftrightarrow -i$ (this swap is the nontrivial element of $\text{Gal}(\mathbb{C}/\mathbb{R}) \cong \mathbb{Z}/2\mathbb{Z}$). But the choice does not change the multiplicative structure: $|z_1 z_2| = |z_1||z_2|$ and $\arg(z_1 z_2) = \arg(z_1) + \arg(z_2)$ regardless of which root we call $i$. The choice only matters when we want to write $z = a + bi$ — that is, when we want to do **addition** in coordinates.

The situation in $\overline{\mathbb{F}_2}$ is exactly analogous, just larger. The compatible hierarchy $\{\alpha_n\}$ is a choice of "which root to call $\alpha$" in every finite subfield simultaneously, subject to cross-field consistency. The Galois group $\text{Gal}(\overline{\mathbb{F}_2}/\mathbb{F}_2)$ is the profinite completion $\hat{\mathbb{Z}}$ (topologically generated by the Frobenius $x \mapsto x^2$), so there are infinitely many automorphisms permuting the roots — compared to just 2 over $\mathbb{C}$. But the principle is the same: the **multiplicative structure is intrinsic** and does not depend on the hierarchy at all. Periods, field membership, Frobenius orbits, the LCM rule, coset sizes — all of these are well-defined properties of the abstract elements, visible from any labeling. The hierarchy is needed only to give a concrete form of **addition**: to write $\alpha^{p_1/q_1} + \alpha^{p_2/q_2}$ as a polynomial in some chosen generator.

In both settings, the practical pattern is the same: work in polar/exponential form when reasoning about multiplicative structure, switch to rectangular/polynomial form when you need addition, and use the chosen basis ($i$, or the hierarchy) to convert between them. The conversion from polynomial to exponential form (the discrete logarithm) is computationally expensive over $\overline{\mathbb{F}_2}$ — unlike $\mathbb{C}$, where $r = |z|$ and $\theta = \arg(z)$ have closed-form expressions. This asymmetry in conversion cost is in fact the basis of some cryptographic systems.

## Representations of Field Elements

With the analogy in mind, we can now describe the "rectangular" representations concretely. These are the forms the library actually computes with — register states, matrix updates, and polynomial arithmetic all live here.

Any root $\alpha \in \overline{\mathbb{F}_2}$ with $\text{Ord}(q) = n$ (so $\alpha \in \mathbb{F}_{2^n}$) can be concretely represented in several equivalent ways:

### Polynomial Quotient Representation

The standard algebraic construction identifies $\mathbb{F}_{2^n}$ with the quotient ring:

$$\mathbb{F}_{2^n} \cong \mathbb{F}_2[x] / \langle P(x) \rangle$$

where $P$ is an irreducible polynomial of degree $n$ over $\mathbb{F}_2$. Under this identification, elements of $\mathbb{F}_{2^n}$ are polynomials of degree less than $n$ with coefficients in $\mathbb{F}_2$, and multiplication is polynomial multiplication modulo $P$. The element $x$ (the coset of the indeterminate) is a root of $P$ in this field.

When $P$ is the minimal polynomial of $\alpha$, there is a natural isomorphism sending $x \mapsto \alpha$. Under this isomorphism, $\alpha$ (the abstract root) and $x$ (the indeterminate modulo $P$) are technically distinct objects but become practically interchangeable — evaluating any polynomial expression at $x$ gives the same result as evaluating it at $\alpha$.

The choice of $P$ matters: different irreducible polynomials of degree $n$ give isomorphic copies of $\mathbb{F}_{2^n}$, but the isomorphism is not canonical. Choosing a different $P$ changes which element plays the role of "$x$" — it's a different primitive element of the same abstract field. The library's register constructors fix this choice: the constructor polynomial $P$ determines the representation.

### Matrix Representation

Under the polynomial quotient representation, multiplication by $\alpha$ is an $\mathbb{F}_2$-linear map on the $n$-dimensional vector space $\mathbb{F}_{2^n} \cong \mathbb{F}_2^n$. This gives a faithful representation:

$$\alpha \mapsto M_\alpha \in \text{GL}(n, \mathbb{F}_2)$$

where $M_\alpha$ is the **companion matrix** of the minimal polynomial of $\alpha$. More generally, any element $\beta \in \mathbb{F}_{2^n}$ maps to a matrix $M_\beta$ whose action on $\mathbb{F}_2^n$ is multiplication by $\beta$.

This is exactly the representation that Galois LFSRs use: the register state is a column vector in $\mathbb{F}_2^n$, with coordinate 0 at the top under the [Matrix Indexing](../conventions/Matrix%20Indexing.md) convention. Clocking the register multiplies (or divides) by $\alpha$ — depending on the shift direction convention (see [Polynomial Conventions — Hardware Conventions](../conventions/Polynomial%20Conventions.md#hardware-conventions)).

### Vector Representation

Elements can also be represented directly as coefficient vectors in $\mathbb{F}_2^n$:

$$\beta = b_0 + b_1 x + b_2 x^2 + \cdots + b_{n-1} x^{n-1} \longleftrightarrow (b_0, b_1, \ldots, b_{n-1})$$

This is the representation the library uses for register states. The vector $(1, 0, 0, \ldots, 0)$ corresponds to the field identity $1$, and $(0, 1, 0, \ldots, 0)$ corresponds to $\alpha$ (the root of the representation polynomial). Clocking the register applies $M_\alpha$ to this vector.

The rest of this document develops the properties that each representation reveals, starting with what is visible purely from the exponential form.

## Properties of Roots

### Period and Field Membership

The period $q$ of a root $\alpha^{p/q}$ determines which subfield it inhabits:

$$\alpha^{p/q} \in \mathbb{F}_{2^n} \iff q \mid 2^n - 1 \iff \text{Ord}_q(2) \mid n$$

where $\text{Ord}_q(2)$ is the multiplicative order of 2 modulo $q$ — the smallest positive $k$ such that $2^k \equiv 1 \pmod{q}$. The **smallest** field containing $\alpha^{p/q}$ is $\mathbb{F}_{2^{\text{Ord}_q(2)}}$.

We write $\text{Ord}(q)$ as shorthand for $\text{Ord}_q(2)$ throughout the library, since the order of 2 is almost always what is meant.

**Example:** Consider $q = 7$. Since $2^3 = 8 \equiv 1 \pmod{7}$, we have $\text{Ord}(7) = 3$. So any root $\alpha^{p/7}$ (with $\gcd(p, 7) = 1$) lives in $\mathbb{F}_{2^3}$. The six such roots ($p = 1, 2, 3, 4, 5, 6$) are exactly the primitive elements of $\mathbb{F}_{2^3}$ — the elements of multiplicative order 7 in the cyclic group $\mathbb{F}_{2^3}^\times \cong \mathbb{Z}/7\mathbb{Z}$.

**Example:** Consider $q = 21 = 3 \cdot 7$. We have $\text{Ord}_{21}(2) = \text{lcm}(\text{Ord}_3(2), \text{Ord}_7(2)) = \text{lcm}(2, 3) = 6$. So roots of period 21 live in $\mathbb{F}_{2^6}$. These roots are neither primitive in $\mathbb{F}_{2^6}$ (since $21 \neq 2^6 - 1 = 63$) nor do they live in any proper subfield (since no proper divisor of 6 gives $21 \mid 2^k - 1$).

### Primitive Roots

A root is **primitive** if it has the form $\alpha^{p/q}$ where $q = 2^k - 1$ for some $k$, placing it as a primitive element of $\mathbb{F}_{2^k}$ (a generator of the full multiplicative group $\mathbb{F}_{2^k}^\times$). Equivalently, $\alpha^{p/q}$ is primitive iff $\text{Ord}(q) = k$ and $q = 2^k - 1$ — the period exhausts the entire multiplicative group of the field.

Primitive roots are special because their minimal polynomials are primitive polynomials (see below), and LFSRs with primitive polynomials achieve the maximum possible period $2^n - 1$. Most of the library's register types are built from primitive polynomials, making primitive roots the common case.

Non-primitive roots arise when a register's characteristic polynomial is reducible or when the polynomial is irreducible but not primitive. They also appear naturally when analyzing CMPRs, where products of roots from different component registers produce roots of composite period.

### The Frobenius Automorphism and Conjugates

The map $\phi: x \mapsto x^2$ is the **Frobenius automorphism** of $\overline{\mathbb{F}_2}$. It fixes $\mathbb{F}_2$ pointwise and permutes the elements of each finite subfield $\mathbb{F}_{2^n}$. The Frobenius is the fundamental symmetry of characteristic-2 algebra — it plays the role that complex conjugation plays over $\mathbb{R}$ (and in fact generates the full Galois group $\text{Gal}(\overline{\mathbb{F}_2}/\mathbb{F}_2) \cong \hat{\mathbb{Z}}$ discussed in the analogy above).

**Frobenius orbits (cyclotomic cosets).** The orbit of an element $\alpha$ under repeated application of $\phi$ is the set:

$$\{\alpha, \alpha^2, \alpha^{2^2}, \alpha^{2^3}, \ldots\}$$

This orbit is finite: if $\alpha \in \mathbb{F}_{2^n}$, then $\alpha^{2^n} = \alpha$ (by Fermat's little theorem for finite fields), so the orbit has length dividing $n$. The orbit length equals the degree of $\alpha$ over $\mathbb{F}_2$ — the degree of its minimal polynomial.

The elements of a Frobenius orbit are called the **conjugates** of $\alpha$ over $\mathbb{F}_2$. In exponential notation, the conjugates of $\alpha^{p/q}$ are:

$$\alpha^{p/q},\ \alpha^{2p/q},\ \alpha^{4p/q},\ \ldots,\ \alpha^{2^{d-1}p/q}$$

where all exponent numerators are taken modulo $q$, and $d = \text{Ord}(q)$ is the orbit length. The set $\{p, 2p, 4p, \ldots, 2^{d-1}p\} \pmod{q}$ is a **cyclotomic coset** of 2 modulo $q$.

**Example:** Over $\mathbb{F}_{2^4}$ with $q = 15$, the element $\alpha^{1/15}$ (a primitive 15th root of unity, i.e. a primitive element of $\mathbb{F}_{2^4}$) has conjugates $\alpha^{1/15}, \alpha^{2/15}, \alpha^{4/15}, \alpha^{8/15}$, corresponding to the cyclotomic coset $\{1, 2, 4, 8\} \pmod{15}$.

**Why conjugates matter.** Conjugates appear together or not at all in GF(2)-linear combinations. Generalized eigenspaces add Jordan-length information to this Frobenius-orbit picture; see [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md). This is because:

1. **Minimal polynomials have GF(2) coefficients.** The minimal polynomial of $\alpha$ has exactly the Frobenius orbit of $\alpha$ as its roots. If one conjugate is a root of a polynomial over $\mathbb{F}_2$, all conjugates are.

2. **The trace function sums over conjugates.** The field trace $\text{Tr}: \mathbb{F}_{2^n} \to \mathbb{F}_2$ is $\text{Tr}(\alpha) = \alpha + \alpha^2 + \alpha^{2^2} + \cdots + \alpha^{2^{n-1}}$, which sums over the full Frobenius orbit (with repetitions if $\alpha$ lives in a proper subfield). The trace is the bridge between field elements and GF(2)-valued sequences: a sequence $s[t] = \text{Tr}(c \cdot \alpha^t)$ is a GF(2) sequence parameterized by a root $\alpha$ and a coefficient $c$.

3. **Root expressions count by cosets.** The [root expression](Root%20Expressions%20and%20LC%20Estimation.md) framework counts roots in units of Frobenius cosets, not individual roots. A coset of size $d$ in $\mathbb{F}_{2^e}$ contributes $d$ to the linear complexity — all $d$ conjugates are either all present or all absent. This is why random cancellation of a full coset has probability only $1/2^d$ (see [Monomial Profile Theory](Monomial%20Profile%20Theory.md) for the statistical model).

### Arithmetic in Exponential Form

The exponential representation makes multiplication and several structural properties immediate, at the cost of making addition opaque (requiring the detour through polynomial form described above).

**Multiplication.** Products of roots correspond to sums of exponents. If two elements have exponential labels $\alpha^{p_1/q_1}$ and $\alpha^{p_2/q_2}$, their product is $\alpha^{p_1/q_1 + p_2/q_2}$, where the sum is taken in $\mathbb{Q}/\mathbb{Z}$ (reduce the combined fraction and extract its period). The resulting element has multiplicative order dividing $\text{lcm}(q_1, q_2)$.

**The LCM rule for field membership.** The smallest field containing two elements $\alpha^{p_1/q_1}$ and $\alpha^{p_2/q_2}$ — and therefore any expression formed from them — has degree $\text{lcm}(\text{Ord}(q_1), \text{Ord}(q_2))$ over $\mathbb{F}_2$. This follows from a more fundamental fact:

$$\text{Ord}(\text{lcm}(q_1, q_2)) = \text{lcm}(\text{Ord}(q_1), \text{Ord}(q_2))$$

*Proof.* Factor each $q_i$ into prime powers $q_i = \prod_j p_j^{a_{ij}}$. Then $\text{lcm}(q_1, q_2) = \prod_j p_j^{\max(a_{1j}, a_{2j})}$. By CRT, $\text{Ord}_q(2) = \text{lcm}_j(\text{Ord}_{p_j^{a_j}}(2))$ for any $q = \prod_j p_j^{a_j}$. Since each prime-power component $p_j^{\max(a_{1j}, a_{2j})}$ of $\text{lcm}(q_1, q_2)$ comes from whichever $q_i$ contributes the larger power of $p_j$, the lcm of the resulting orders equals the lcm of the orders from $q_1$ and $q_2$ individually. $\square$

This rule is pervasive: it governs the field degree needed for [root expression](Root%20Expressions%20and%20LC%20Estimation.md) computations, determines where products of roots from different CMPR blocks land, and underlies the period analysis of composite registers.

**Counting elements of a given period.** The number of elements of multiplicative order exactly $q$ in $\overline{\mathbb{F}_2}$ is $\phi(q)$, where $\phi$ is Euler's totient function. These $\phi(q)$ elements partition into $\phi(q) / \text{Ord}(q)$ Frobenius orbits, each of size $\text{Ord}(q)$ — so there are exactly $\phi(q)/\text{Ord}(q)$ irreducible polynomials over $\mathbb{F}_2$ whose roots have period $q$. For the primitive case $q = 2^n - 1$, this recovers the count $\phi(2^n - 1)/n$ of primitive polynomials of degree $n$.

**The Carmichael function and the multiplicative group.** The multiplicative group $\mathbb{F}_{2^n}^\times \cong \mathbb{Z}/(2^n - 1)\mathbb{Z}$ is cyclic, so every element's order divides $2^n - 1$. The Carmichael function $\lambda(m)$ gives the exponent of the group $(\mathbb{Z}/m\mathbb{Z})^\times$ — the largest order of any element. For cyclic groups, $\lambda(m) = \phi(m)$ when $m$ is 1, 2, 4, an odd prime power, or twice an odd prime power; otherwise $\lambda(m) < \phi(m)$.

In the register context, $\lambda(2^n - 1)$ is less directly relevant (since $\mathbb{F}_{2^n}^\times$ is cyclic, every divisor of $2^n - 1$ occurs as an element order). The Carmichael function becomes important when working with $(\mathbb{Z}/q\mathbb{Z})^\times$ for composite $q$ — for instance, the decimation-based primitive polynomial enumeration in `NumberTheory.py` uses $\lambda(p^k)$ to identify primitive elements modulo prime powers and then lifts them via CRT. The key identity is:

$$\lambda(p^k) = \begin{cases} p^{k}/2 & \text{if } p^k \in \{2, 4\} \\ p^{k-2} & \text{if } p = 2, k \geq 3 \\ p^{k-1}(p-1) & \text{if } p \text{ is odd} \end{cases}$$

**Exponentiation (decimation).** Raising a root to a power $d$ maps $\alpha^{p/q} \mapsto \alpha^{dp/q}$ (reducing $dp \bmod q$). In sequence terms, this corresponds to **decimation**: sampling every $d$-th element of the output sequence. Decimation by $d$ permutes the roots within $\mathbb{F}_{2^n}$ (when $\gcd(d, 2^n - 1) = 1$) and sends one primitive polynomial to another. This is the mechanism behind `fast_decimate` and `generate_primitive_polynomials` in the library. For the full treatment — including the inverse operation (expansion), effects on linear complexity, root fans, and Jordan structure — see [Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md).

### What Each Form Reveals

**Exponential form** (hierarchy-free — these properties are intrinsic):
- **Period and field membership:** Immediately visible from the denominator $q$ and $\text{Ord}(q)$.
- **Multiplication:** Exponents add in $\mathbb{Q}/\mathbb{Z}$.
- **Frobenius orbits:** Conjugates are $\alpha^{2^k p/q}$; the orbit structure is determined by the cyclotomic coset of $p$ mod $q$.
- **Minimal polynomial degree:** Equals the orbit size $\text{Ord}(q)$.
- **Subfield containment:** $\alpha^{p/q} \in \mathbb{F}_{2^m}$ iff $\text{Ord}(q) \mid m$.
- **Root counting and LC bounds:** The [root expression](Root%20Expressions%20and%20LC%20Estimation.md) framework operates entirely in this world, tracking coset sizes and periods without materializing any polynomial representatives.

**Polynomial form** (requires a chosen hierarchy and representation polynomial):
- **Addition:** XOR of coefficient vectors.
- **Concrete state vectors:** Register states are vectors in $\mathbb{F}_2^n$, directly corresponding to polynomials in $\alpha$.
- **Matrix action:** Clocking a register applies $M_\alpha$ (the companion matrix) to the state vector.
- **Equation solving:** Linear algebra over $\mathbb{F}_2$ — Gaussian elimination, LU decomposition, Gröbner bases — all operate on the polynomial/vector form.
- **Evaluation and interpolation:** Computing outputs, building equation systems for [algebraic attacks](Algebraic%20Attacks.md).

### Working with Both Forms: A Concrete Example

Consider $\mathbb{F}_{2^3}$ with primitive polynomial $P(x) = x^3 + x + 1$ (coefficient list `[1,1,0,1]`). Choosing $\alpha = \alpha^{1/7}$ as the root of $P$ fixes the hierarchy at level 3. Every nonzero element has an exponential label $\alpha^{j/7}$ and a polynomial representative in $\mathbb{F}_2[\alpha]/\langle P(\alpha) \rangle$:

| Exponential | Power | Polynomial | Vector |
|---|---|---|---|
| $\alpha^{1/7}$ | $\alpha^1$ | $\alpha$ | $(0,1,0)$ |
| $\alpha^{2/7}$ | $\alpha^2$ | $\alpha^2$ | $(0,0,1)$ |
| $\alpha^{3/7}$ | $\alpha^3$ | $\alpha + 1$ | $(1,1,0)$ |
| $\alpha^{4/7}$ | $\alpha^4$ | $\alpha^2 + \alpha$ | $(0,1,1)$ |
| $\alpha^{5/7}$ | $\alpha^5$ | $\alpha^2 + \alpha + 1$ | $(1,1,1)$ |
| $\alpha^{6/7}$ | $\alpha^6$ | $\alpha^2 + 1$ | $(1,0,1)$ |

(The polynomial column is computed by reducing powers of $\alpha$ modulo $P$: $\alpha^3 = \alpha + 1$, then $\alpha^4 = \alpha \cdot \alpha^3 = \alpha(\alpha+1) = \alpha^2 + \alpha$, and so on.)

**Multiplication (easy in exponential form).** What is $\alpha^{1/7} \cdot \alpha^{3/7}$? In exponential form: $1/7 + 3/7 = 4/7$, so the answer is $\alpha^{4/7}$. No computation needed — the representation makes it trivial. In polynomial form, the same result requires multiplying $\alpha \cdot (\alpha + 1) = \alpha^2 + \alpha$ and checking the table to verify this is indeed $\alpha^4$.

**Addition (requires polynomial form).** What is $\alpha^{1/7} + \alpha^{2/7}$? The exponential labels give no hint. Converting to vectors: $(0,1,0) \oplus (0,0,1) = (0,1,1)$. Reading off the table: $(0,1,1) = \alpha^2 + \alpha = \alpha^{4/7}$. The fact that $\alpha^{1/7} + \alpha^{2/7} = \alpha^{4/7}$ is not predictable from the exponential labels — it depends entirely on the choice of representation polynomial $P$.

**Minimal polynomial (using both forms together).** The minimal polynomial of $\alpha^{3/7}$ has roots at its Frobenius conjugates. From the exponential form, the cyclotomic coset of 3 mod 7 is $\{3, 6, 5\}$ (since $3 \cdot 2 = 6$, $6 \cdot 2 = 12 \equiv 5$, $5 \cdot 2 = 10 \equiv 3 \pmod{7}$). So the minimal polynomial is:

$$(x + \alpha^{3/7})(x + \alpha^{6/7})(x + \alpha^{5/7})$$

The exponential form identified the roots; to expand this product into coefficients over $\mathbb{F}_2$, we switch to polynomial form: $(\alpha+1)$, $(\alpha^2+1)$, $(\alpha^2+\alpha+1)$. Expanding and reducing modulo 2 gives $x^3 + x^2 + 1$ — the other primitive polynomial of degree 3.

### Where the Library's Frameworks Live

The [root expression](Root%20Expressions%20and%20LC%20Estimation.md) framework lives almost entirely in the exponential world — it tracks periods and coset sizes without ever materializing polynomial representatives. The [resolvent analysis](Resolvent%20Analysis.md), by contrast, works in the polynomial/power-series world. The primal/dual framework in [Polynomial Conventions](../conventions/Polynomial%20Conventions.md) connects the two: the primal perspective (roots of a polynomial) is exponential thinking, while the dual perspective (convolution annihilation) is polynomial thinking.

The `fast_decimate` function in `NumberTheory.py` illustrates the alternation pattern in practice: it starts with a primitive polynomial (polynomial form), generates a sequence by running a register (vector form), decimates the sequence (an operation defined by exponential-form reasoning — decimation by $d$ maps $\alpha^{p/q} \mapsto \alpha^{dp/q}$), then applies Berlekamp-Massey to recover the result as a polynomial. Each step uses whichever representation makes that particular operation natural.

## Polynomial Types

Minimal, characteristic, and primitive polynomials each play a distinct role in register analysis. See [Polynomial Types](Polynomial%20Types.md) for their definitions, relationships, and how they connect to the concepts developed above.

## Connections to Other Documents

- **[Polynomial Conventions](../conventions/Polynomial%20Conventions.md):** Defines the primal/dual framework — how a polynomial relates to a sequence through its roots (primal) or through convolution (dual). The algebraic closure is where those roots live.
- **[Polynomial Types](Polynomial%20Types.md):** Minimal, characteristic, and primitive polynomials — their definitions, relationships, and role in register analysis.
- **[Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md):** Uses the subfield structure of $\overline{\mathbb{F}_2}$ to bound linear complexity. Root expressions track how many roots from each subfield $\mathbb{F}_{2^e}$ can appear in a signal.
- **[Monomial Profile Theory](Monomial%20Profile%20Theory.md):** Parallels root expressions for ANF degree. The full-coset degeneracy (§1 of that document) is a direct consequence of $\alpha^{2^e - 1} = 1$ in $\mathbb{F}_{2^e}$.
- **[Resolvent Analysis](Resolvent%20Analysis.md):** The resolvent matrix has entries in $\text{GF}(2)[[D]]$, and its eigenstructure feeds the root expression computation.
- **[Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md):** Develops the root-level theory of decimation ($\alpha^{p/q} \mapsto \alpha^{dp/q}$) and expansion (root fans), their effects on LC and Jordan structure, and the duality between them.
- **[Algebraic Closure](Algebraic%20Closure.md):** Explores the algebraic closure from the perspective of completions and localizations — relating the field-theoretic and string-theoretic (sequence) viewpoints.
