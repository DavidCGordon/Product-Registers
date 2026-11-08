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

The nonzero elements of $\overline{\mathbb{F}_2}$ are named in this library by a fractional exponent, $\alpha^{p/q}$. The denominator $q$ is the multiplicative order of the element and involves no choice. The numerator $p$ is defined relative to a choice of one generator per finite field. This section treats the denominator; the next section makes the choice of generators and then defines the notation.

### Period and Field Membership

Let $\beta$ be a nonzero element of $\overline{\mathbb{F}_2}$. It lies in some $\mathbb{F}_{2^n}$, whose multiplicative group $\mathbb{F}_{2^n}^\times$ is cyclic of order $2^n - 1$, so $\beta$ has a finite multiplicative order $q$ dividing $2^n - 1$. In particular $q$ is odd. We call $q$ the **period** of $\beta$: it is the period of the sequence $t \mapsto \beta^t$. It is not the period of $\beta$ under repeated squaring; that is the length of its Frobenius orbit, treated under [The Frobenius Automorphism and Conjugates](#the-frobenius-automorphism-and-conjugates).

The period determines which subfields contain $\beta$:

$$\beta \in \mathbb{F}_{2^n} \iff q \mid 2^n - 1 \iff \text{Ord}_q(2) \mid n$$

where $\text{Ord}_q(2)$ is the multiplicative order of 2 modulo $q$ — the smallest positive $k$ such that $2^k \equiv 1 \pmod{q}$. For the first equivalence, the nonzero elements of $\mathbb{F}_{2^n}$ are exactly the solutions of $x^{2^n - 1} = 1$ in $\overline{\mathbb{F}_2}$ (all $2^n - 1$ of them satisfy it, and the equation has at most $2^n - 1$ solutions), and $\beta^{2^n - 1} = 1$ exactly when $q \mid 2^n - 1$. For the second, $q \mid 2^n - 1$ says $2^n \equiv 1 \pmod q$, which holds exactly when the order of 2 modulo $q$ divides $n$. The **smallest** field containing $\beta$ is therefore $\mathbb{F}_{2^{\text{Ord}_q(2)}}$. In the other direction, since $\mathbb{F}_{2^n}^\times$ is cyclic of order $2^n - 1$, an element of multiplicative order $q$ exists in $\mathbb{F}_{2^n}$ precisely when $q \mid 2^n - 1$.

We write $\text{Ord}(q)$ as shorthand for $\text{Ord}_q(2)$ throughout the library, since the order of 2 is almost always what is meant.

**Example:** Consider $q = 7$. Since $2^3 = 8 \equiv 1 \pmod{7}$, we have $\text{Ord}(7) = 3$. So any element of period 7 lives in $\mathbb{F}_{2^3}$. There are six of them, and they are exactly the primitive elements of $\mathbb{F}_{2^3}$ — the elements of multiplicative order 7 in the cyclic group $\mathbb{F}_{2^3}^\times \cong \mathbb{Z}/7\mathbb{Z}$. In the notation defined below they are $\alpha^{p/7}$ for $p = 1, 2, 3, 4, 5, 6$.

**Example:** Consider $q = 21 = 3 \cdot 7$. We have $\text{Ord}_{21}(2) = \text{lcm}(\text{Ord}_3(2), \text{Ord}_7(2)) = \text{lcm}(2, 3) = 6$. So roots of period 21 live in $\mathbb{F}_{2^6}$. These roots are neither primitive in $\mathbb{F}_{2^6}$ (since $21 \neq 2^6 - 1 = 63$) nor do they live in any proper subfield (since no proper divisor of 6 gives $21 \mid 2^k - 1$).

## The Compatible Primitive Element Hierarchy

The numerator of the exponential notation is defined relative to a "compatible system of primitive elements" — one primitive element $\alpha_n$ per field $\mathbb{F}_{2^n}$, chosen so that they agree across subfield embeddings. A **primitive element** of $\mathbb{F}_{2^n}$ is a generator of the cyclic group $\mathbb{F}_{2^n}^\times$. This section states the compatibility condition, defines the notation $\alpha^{p/q}$ from it, proves that a compatible system exists, and says what depends on the choice of system.

### The Compatibility Condition

When $m \mid n$, the subfield $\mathbb{F}_{2^m}$ embeds into $\mathbb{F}_{2^n}$. The **norm map** $N_{n/m}: \mathbb{F}_{2^n}^\times \to \mathbb{F}_{2^m}^\times$ is the product of the conjugates of $x$ over $\mathbb{F}_{2^m}$, and it is a power map:

$$N_{n/m}(x) = \prod_{j=0}^{n/m - 1} x^{2^{mj}} = x^{(2^n - 1)/(2^m - 1)}, \qquad \text{since } \sum_{j=0}^{n/m - 1} 2^{mj} = \frac{2^n - 1}{2^m - 1}.$$

A primitive element $\alpha_n$ of $\mathbb{F}_{2^n}$ is not an element of $\mathbb{F}_{2^m}$ when $m < n$. The norm map sends it to a primitive element of $\mathbb{F}_{2^m}$ (shown under [Existence](#existence)), but not necessarily to the one that was chosen for $\mathbb{F}_{2^m}$.

A system $\{\alpha_n\}_{n \geq 1}$ of primitive elements is **compatible** if for every pair $m \mid n$:

$$\alpha_n^{(2^n - 1)/(2^m - 1)} = \alpha_m$$

This ensures that the exponential label of an element does not depend on which ambient field we view it in: the element $\alpha_m^j$ in $\mathbb{F}_{2^m}$ has the same identity whether we name it using $\alpha_m$ or by embedding into $\mathbb{F}_{2^n}$ and using $\alpha_n$. The definition below makes this precise.

### Fractional Exponent Notation

Fix a compatible system $\{\alpha_n\}$; one exists by the [next subsection](#existence). Let $q > 1$ be an odd integer and $p$ an integer with $\gcd(p, q) = 1$ and $0 < p < q$. Define

$$\alpha^{p/q} := \alpha_n^{\,p(2^n - 1)/q} \qquad \text{for any } n \text{ with } q \mid 2^n - 1.$$

Such an $n$ exists, for instance $n = \text{Ord}(q)$, and the exponent is an integer because $q \mid 2^n - 1$.

*The definition does not depend on $n$.* If $m \mid n$ and $q \mid 2^m - 1$, then by compatibility

$$\alpha_n^{\,p(2^n - 1)/q} = \left(\alpha_n^{(2^n - 1)/(2^m - 1)}\right)^{p(2^m - 1)/q} = \alpha_m^{\,p(2^m - 1)/q}.$$

Any $n$ with $q \mid 2^n - 1$ is a multiple of $m = \text{Ord}(q)$ by the equivalence above, so every admissible $n$ gives the same element as $\text{Ord}(q)$ does.

*The element has period $q$.* Since $\alpha_n$ has order $2^n - 1$, the power $\alpha_n^{(2^n - 1)/q}$ has order $q$, and raising an element of order $q$ to a power $p$ coprime to $q$ leaves its order unchanged.

*Every nonzero element other than 1 has exactly one such name.* An element $\beta$ of period $q > 1$ lies in $\mathbb{F}_{2^n}$ for $n = \text{Ord}(q)$, so $\beta = \alpha_n^j$ for a unique $j$ with $0 < j < 2^n - 1$. The order of $\alpha_n^j$ is $(2^n - 1)/\gcd(j, 2^n - 1)$, and this equals $q$ exactly when $j = p(2^n - 1)/q$ with $\gcd(p, q) = 1$; the range of $j$ gives $0 < p < q$. The element 1 is the case $q = 1$, for which no such $p$ exists; it is written $1$.

*Exponents are fractions modulo 1.* Write $\mathbb{Z}_{(2)}$ for the ring of fractions with odd denominator. The formula $a/q \mapsto \alpha_n^{\,a(2^n - 1)/q}$ makes sense for any integer $a$ and odd $q$ with $q \mid 2^n - 1$. It is unchanged when the fraction is reduced to lowest terms, and unchanged when an integer is added to $a/q$, because $\alpha_n^{2^n - 1} = 1$. So it is a map on the group $\mathbb{Z}_{(2)}/\mathbb{Z}$ of odd-denominator fractions modulo 1, and it sends sums to products, because for a common $n$ the exponents of $\alpha_n$ add. This is the sense in which exponents add under [Arithmetic in Exponential Form](#arithmetic-in-exponential-form).

Period, field membership, Frobenius conjugates, and multiplicative arithmetic are read from this notation in the [Properties of Roots](#properties-of-roots) section below. The sections before it supply what is needed to work with these elements concretely: the existence of a compatible hierarchy, a comparison with $\mathbb{C}$ that sets the exponential and polynomial representations side by side, and the polynomial, matrix, and vector forms the library computes with.

### Existence

**Claim.** A compatible system of primitive elements exists.

*Proof.* The construction chooses generators along the chain of fields $\mathbb{F}_{2^{k!}}$, $k = 1, 2, 3, \ldots$, and obtains every other $\alpha_m$ as a norm.

1. *The norm is a surjection of cyclic groups.* For $m \mid n$, the map $N_{n/m}(x) = x^{(2^n - 1)/(2^m - 1)}$ is a homomorphism on the cyclic group $\mathbb{F}_{2^n}^\times$ of order $2^n - 1$. It sends a generator to an element of order $2^m - 1$, so its image is the subgroup of order $2^m - 1$. A cyclic group has only one subgroup of each order, and $\mathbb{F}_{2^m}^\times$ is a subgroup of that order, so the image is $\mathbb{F}_{2^m}^\times$.

2. *A surjection of cyclic groups sends generators onto generators.* Fix a generator $g$ of $\mathbb{F}_{2^n}^\times$ and let $h = N_{n/m}(g)$, which generates the image. Writing $N = 2^n - 1$ and $M = 2^m - 1$, the map sends $g^u \mapsto h^u$, so on exponents it is the reduction $\mathbb{Z}/N \to \mathbb{Z}/M$. The generators are the powers $g^u$ with $u$ a unit modulo $N$, respectively $h^u$ with $u$ a unit modulo $M$, and for $M \mid N$ reduction maps the units $(\mathbb{Z}/N)^\times$ onto the units $(\mathbb{Z}/M)^\times$. To see this, take $u$ prime to $M$ and put $u' = u + tM$, where $t$ is the product of the primes dividing $N$ that divide neither $u$ nor $M$. For each prime $p \mid N$: if $p \mid M$ then $p \nmid u$, so $p \nmid u'$; if $p \nmid M$ and $p \mid u$ then $p \nmid tM$, so $p \nmid u'$; and if $p \nmid M$ and $p \nmid u$ then $p \mid t$, so $u' \equiv u \not\equiv 0 \pmod p$. Thus $u'$ is a unit modulo $N$ that reduces to $u$ modulo $M$. So the image of a generator of $\mathbb{F}_{2^n}^\times$ under $N_{n/m}$ is a generator of $\mathbb{F}_{2^m}^\times$, and every generator of $\mathbb{F}_{2^m}^\times$ is the norm of some generator of $\mathbb{F}_{2^n}^\times$.

3. *Choose generators along the factorials.* Set $\alpha_{1!} = 1$, the only element of $\mathbb{F}_2^\times$. Given a generator $\alpha_{k!}$ of $\mathbb{F}_{2^{k!}}^\times$, step 2 applied to $k! \mid (k+1)!$ gives a generator $\alpha_{(k+1)!}$ of $\mathbb{F}_{2^{(k+1)!}}^\times$ with $N_{(k+1)!/k!}(\alpha_{(k+1)!}) = \alpha_{k!}$.

4. *Define the rest by norms.* Every $m$ divides $k!$ for $k \geq m$. Define $\alpha_m := N_{k!/m}(\alpha_{k!})$ for any such $k$.

5. *This is well defined and compatible.* Norms compose, because the exponents multiply: for $m \mid n \mid l$,
   $$N_{n/m} \circ N_{l/n} = N_{l/m}, \qquad \frac{2^l - 1}{2^n - 1} \cdot \frac{2^n - 1}{2^m - 1} = \frac{2^l - 1}{2^m - 1}.$$
   By step 3 and composition, $N_{k'!/k!}(\alpha_{k'!}) = \alpha_{k!}$ for $k \leq k'$, so $N_{k'!/m}(\alpha_{k'!}) = N_{k!/m}(\alpha_{k!})$ and the definition in step 4 does not depend on $k$. For $m \mid n$, choose $k$ with $n \mid k!$; then $N_{n/m}(\alpha_n) = N_{n/m}(N_{k!/n}(\alpha_{k!})) = N_{k!/m}(\alpha_{k!}) = \alpha_m$, which is the compatibility condition. Each $\alpha_m$ is a generator by step 2. $\square$

The system is not unique, because step 3 can offer more than one admissible generator. Already at $k = 2$, both generators of $\mathbb{F}_{2^2}^\times$ have norm $x^3 = 1 = \alpha_1$. More generally, any generator $\gamma$ of $\mathbb{F}_{2^m}^\times$ occurs as $\alpha_m$ in some compatible system: by step 2 choose a generator $\alpha_{m!}$ with $N_{m!/m}(\alpha_{m!}) = \gamma$, take $\alpha_{k!}$ for $k < m$ to be its norms, and continue step 3 from $k = m$. What does and does not depend on the choice is taken up under [The Choice of Roots](#the-choice-of-roots).

### Why the Hierarchy Matters

The compatible hierarchy is what connects the exponential and polynomial representations.

**Naming across fields.** Without compatibility, saying "$\alpha^{3/7}$" is ambiguous: it depends on which primitive element we chose for $\mathbb{F}_{2^3}$. With a compatible hierarchy, the name is universal — the element $\alpha^{3/7}$ is the same whether we work in $\mathbb{F}_{2^3}$, $\mathbb{F}_{2^6}$, $\mathbb{F}_{2^{12}}$, or any other field containing $\mathbb{F}_{2^3}$.

**Enabling addition.** The exponential representation has no direct formula for addition: $\alpha^{p_1/q_1} + \alpha^{p_2/q_2}$ cannot in general be simplified to $\alpha^{p_3/q_3}$ without computing in a concrete representation. To add two elements, one must:

1. Identify the smallest common field $\mathbb{F}_{2^N}$ containing both. By [Period and Field Membership](#period-and-field-membership), $\mathbb{F}_{2^N}$ contains both exactly when $\text{Ord}(q_1) \mid N$ and $\text{Ord}(q_2) \mid N$, so the smallest is $N = \text{lcm}(\text{Ord}(q_1), \text{Ord}(q_2))$.
2. Choose a representation — a primitive polynomial $P_N$ of degree $N$ — which fixes $\alpha_N$ as the root of $P_N$.
3. Express both elements as polynomials of degree $< N$ in $\alpha_N$ (their vector representations).
4. Add the vectors (XOR the coefficient arrays).
5. If needed, convert the result back to exponential form (a discrete logarithm).

Steps 2–4 are polynomial arithmetic modulo $P_N$. The compatible hierarchy ensures that the embedding in step 3 is consistent: by the definition of the notation, $\alpha^{p_1/q_1}$ is $\alpha_N^{p_1(2^N - 1)/q_1}$, which can be computed as a polynomial in $\alpha_N$ by repeated squaring modulo $P_N$.

Multiplication is immediate in exponential form, while addition requires this detour through polynomial form. Polar and rectangular coordinates on $\mathbb{C}$ show the same asymmetry, and the next section sets the two cases side by side.

## The Two Faces of $\overline{\mathbb{F}_2}$: Analogy with $\mathbb{C}$

The interplay between exponential and polynomial representations in $\overline{\mathbb{F}_2}$ parallels the interplay between polar and rectangular coordinates in $\mathbb{C}$. The table below states the correspondence row by row; the discussion after it uses the last row to say what the compatible hierarchy provides and when each representation is the right tool.

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

In $\overline{\mathbb{F}_2}$, the compatible hierarchy $\{\alpha_n\}$ is the corresponding choice: a choice of "which root to call $\alpha$" in every finite subfield simultaneously, subject to cross-field consistency. The Galois group $\text{Gal}(\overline{\mathbb{F}_2}/\mathbb{F}_2)$ is the profinite completion $\hat{\mathbb{Z}}$ (topologically generated by the Frobenius $x \mapsto x^2$), so there are infinitely many automorphisms permuting the roots — compared to just 2 over $\mathbb{C}$.

Some properties of an element do not depend on the hierarchy: its period $q$, the degree $\text{Ord}(q)$, the size of its Frobenius orbit, and which subfields contain it. These are properties of the abstract element. The numerator $p$ is not: the same element can have different numerators in different hierarchies.

Here the analogy with $\mathbb{C}$ reaches its limit. Over $\mathbb{C}$ the two choices of $i$ are exchanged by a field automorphism, so additive identities do not depend on the choice. In $\overline{\mathbb{F}_2}$, additive identities between labels do depend on the hierarchy. Take two hierarchies that differ at level 3. In the first, $\alpha = \alpha^{1/7}$ is a root of $x^3 + x + 1$; then $\alpha^3 = \alpha + 1$ and $\alpha^4 = \alpha^2 + \alpha$, so $\alpha + \alpha^2 = \alpha^4$. In the second, $\alpha' = \alpha^{1/7}$ is a root of $x^3 + x^2 + 1$; then $\alpha'^3 = \alpha'^2 + 1$, $\alpha'^4 = \alpha'^2 + \alpha' + 1$, $\alpha'^5 = \alpha' + 1$, and $\alpha'^6 = \alpha'^2 + \alpha'$, so $\alpha' + \alpha'^2 = \alpha'^6$. The labelled sum $\alpha^{1/7} + \alpha^{2/7}$ is $\alpha^{4/7}$ in the first hierarchy and $\alpha^{6/7}$ in the second.

In both settings, the practical pattern is the same: work in polar/exponential form when reasoning about multiplicative structure, switch to rectangular/polynomial form when you need addition, and use the chosen basis ($i$, or the hierarchy) to convert between them. The conversion from polynomial to exponential form (the discrete logarithm) is computationally expensive over $\overline{\mathbb{F}_2}$ — unlike $\mathbb{C}$, where $r = |z|$ and $\theta = \arg(z)$ have closed-form expressions. This asymmetry in conversion cost is in fact the basis of some cryptographic systems.

## Representations of Field Elements

This section describes the "rectangular" representations concretely. These are the forms the library computes with: register states, matrix updates, and polynomial arithmetic are all expressed in them.

Any root $\alpha \in \overline{\mathbb{F}_2}$ of multiplicative order $q$ with $\text{Ord}(q) = n$ (so $\alpha \in \mathbb{F}_{2^n}$) can be concretely represented in several equivalent ways.

### Polynomial Quotient Representation

The standard algebraic construction identifies $\mathbb{F}_{2^n}$ with the quotient ring:

$$\mathbb{F}_{2^n} \cong \mathbb{F}_2[x] / \langle P(x) \rangle$$

where $P$ is an irreducible polynomial of degree $n$ over $\mathbb{F}_2$. Under this identification, elements of $\mathbb{F}_{2^n}$ are polynomials of degree less than $n$ with coefficients in $\mathbb{F}_2$, and multiplication is polynomial multiplication modulo $P$. The element $x$ (the coset of the indeterminate) is a root of $P$ in this field.

When $P$ is the minimal polynomial of $\alpha$, there is a natural isomorphism sending $x \mapsto \alpha$. Under this isomorphism, $\alpha$ (the abstract root) and $x$ (the indeterminate modulo $P$) are technically distinct objects but become practically interchangeable — evaluating any polynomial expression at $x$ gives the same result as evaluating it at $\alpha$.

The choice of $P$ matters: different irreducible polynomials of degree $n$ give isomorphic copies of $\mathbb{F}_{2^n}$, but the isomorphism is not canonical. Choosing a different $P$ changes which element plays the role of "$x$" — it is a different generator of the field over $\mathbb{F}_2$. (It need not be a primitive element: modulo the irreducible polynomial $x^4 + x^3 + x^2 + x + 1$, the element $x$ has multiplicative order 5.) The library's register constructors fix this choice: the constructor polynomial $P$ determines the representation.

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

The period $q$ of a root and the degree $\text{Ord}(q)$ of the smallest field containing it were established under [Period and Field Membership](#period-and-field-membership). The subsections below build on them.

### Primitive Roots

A root is **primitive** if it has the form $\alpha^{p/q}$ where $q = 2^k - 1$ for some $k$, placing it as a primitive element of $\mathbb{F}_{2^k}$ (a generator of the full multiplicative group $\mathbb{F}_{2^k}^\times$). Equivalently, $\alpha^{p/q}$ is primitive iff $\text{Ord}(q) = k$ and $q = 2^k - 1$ — the period exhausts the entire multiplicative group of the field.

Primitive roots are special because their minimal polynomials are primitive polynomials (see [Polynomial Types](Polynomial%20Types.md)), and LFSRs with primitive polynomials achieve the maximum possible period $2^n - 1$. Most of the library's register types are built from primitive polynomials, making primitive roots the common case.

Non-primitive roots arise when a register's characteristic polynomial is reducible or when the polynomial is irreducible but not primitive. They also appear naturally when analyzing CMPRs, where products of roots from different component registers produce roots of composite period.

### The Frobenius Automorphism and Conjugates

The map $\phi: x \mapsto x^2$ is the **Frobenius automorphism** of $\overline{\mathbb{F}_2}$. It fixes $\mathbb{F}_2$ pointwise and permutes the elements of each finite subfield $\mathbb{F}_{2^n}$. It topologically generates the Galois group $\text{Gal}(\overline{\mathbb{F}_2}/\mathbb{F}_2) \cong \hat{\mathbb{Z}}$ of [The Choice of Roots](#the-choice-of-roots), as complex conjugation generates $\text{Gal}(\mathbb{C}/\mathbb{R})$.

**Frobenius orbits (cyclotomic cosets).** The orbit of an element $\alpha$ under repeated application of $\phi$ is the set:

$$\{\alpha, \alpha^2, \alpha^{2^2}, \alpha^{2^3}, \ldots\}$$

This orbit is finite: if $\alpha \in \mathbb{F}_{2^n}$, then $\alpha^{2^n} = \alpha$ (by Fermat's little theorem for finite fields), so the orbit has length dividing $n$. The orbit length equals the degree of $\alpha$ over $\mathbb{F}_2$ — the degree of its minimal polynomial.

The elements of a Frobenius orbit are called the **conjugates** of $\alpha$ over $\mathbb{F}_2$. In exponential notation, the conjugates of $\alpha^{p/q}$ are:

$$\alpha^{p/q},\ \alpha^{2p/q},\ \alpha^{4p/q},\ \ldots,\ \alpha^{2^{d-1}p/q}$$

where all exponent numerators are taken modulo $q$, and $d = \text{Ord}(q)$ is the orbit length. The set $\{p, 2p, 4p, \ldots, 2^{d-1}p\} \pmod{q}$ is a **cyclotomic coset** of 2 modulo $q$.

**Example:** Over $\mathbb{F}_{2^4}$ with $q = 15$, the element $\alpha^{1/15}$ (a primitive 15th root of unity, i.e. a primitive element of $\mathbb{F}_{2^4}$) has conjugates $\alpha^{1/15}, \alpha^{2/15}, \alpha^{4/15}, \alpha^{8/15}$, corresponding to the cyclotomic coset $\{1, 2, 4, 8\} \pmod{15}$.

**Why conjugates matter.** In the root representation of an $\mathbb{F}_2$-valued sequence, conjugates appear together or not at all. Precisely, let

$$s[t] = \sum_\alpha A_\alpha\, \alpha^t$$

be a sequence with values in $\mathbb{F}_2$, where the sum runs over finitely many distinct nonzero $\alpha \in \overline{\mathbb{F}_2}$ and the coefficients $A_\alpha$ lie in $\overline{\mathbb{F}_2}$. Since $s[t] \in \mathbb{F}_2$, we have $s[t]^2 = s[t]$, and squaring is additive in characteristic 2, so

$$\sum_\alpha A_\alpha^2\, (\alpha^2)^t = \sum_\alpha A_\alpha\, \alpha^t.$$

The sequences $t \mapsto \alpha^t$ for distinct $\alpha$ are linearly independent (a Vandermonde determinant), so the representation is unique, and comparing the coefficient of $(\alpha^2)^t$ on both sides gives $A_{\alpha^2} = A_\alpha^2$. Hence $A_\alpha = 0$ if and only if $A_{\alpha^2} = 0$: a root is present exactly when every conjugate is.

Three places where this shows up:

1. **Minimal polynomials have $\mathbb{F}_2$ coefficients.** The minimal polynomial of $\alpha$ has exactly the Frobenius orbit of $\alpha$ as its roots. If one conjugate is a root of a polynomial over $\mathbb{F}_2$, all conjugates are.

2. **The trace function sums over conjugates.** The field trace $\text{Tr}: \mathbb{F}_{2^n} \to \mathbb{F}_2$ is $\text{Tr}(\alpha) = \alpha + \alpha^2 + \alpha^{2^2} + \cdots + \alpha^{2^{n-1}}$, which sums over the full Frobenius orbit (with repetitions if $\alpha$ lives in a proper subfield). The trace is the bridge between field elements and $\mathbb{F}_2$-valued sequences: a sequence $s[t] = \text{Tr}(c \cdot \alpha^t)$ is an $\mathbb{F}_2$-valued sequence parameterized by a root $\alpha$ and a coefficient $c$.

3. **Root expressions count by cosets.** The [root expression](Root%20Expressions%20and%20LC%20Estimation.md) framework counts roots in units of Frobenius cosets, not individual roots. A coset of size $d$ in $\mathbb{F}_{2^e}$ contributes $d$ to the linear complexity — all $d$ conjugates are either all present or all absent. This is why random cancellation of an entire Frobenius coset has probability only $1/2^d$ (see [Monomial Profile Theory](Monomial%20Profile%20Theory.md) for the statistical model).

Generalized eigenspaces add Jordan-length information to this Frobenius-orbit picture; see [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md).

### Arithmetic in Exponential Form

The exponential representation makes multiplication and several structural properties immediate, at the cost of making addition opaque (requiring the detour through polynomial form described under [Why the Hierarchy Matters](#why-the-hierarchy-matters)).

**Multiplication.** Products of roots correspond to sums of exponents. If two elements have exponential labels $\alpha^{p_1/q_1}$ and $\alpha^{p_2/q_2}$, their product is $\alpha^{p_1/q_1 + p_2/q_2}$, where the sum is taken in $\mathbb{Z}_{(2)}/\mathbb{Z}$, the odd-denominator fractions modulo 1 of [Fractional Exponent Notation](#fractional-exponent-notation) (reduce the combined fraction and extract its period). The resulting element has multiplicative order dividing $\text{lcm}(q_1, q_2)$.

**The LCM rule for field membership.** The smallest field containing two elements $\alpha^{p_1/q_1}$ and $\alpha^{p_2/q_2}$ — and therefore any expression formed from them — has degree $\text{lcm}(\text{Ord}(q_1), \text{Ord}(q_2))$ over $\mathbb{F}_2$, as shown in step 1 of [Why the Hierarchy Matters](#why-the-hierarchy-matters). The same degree is the order of 2 modulo $\text{lcm}(q_1, q_2)$:

$$\text{Ord}(\text{lcm}(q_1, q_2)) = \text{lcm}(\text{Ord}(q_1), \text{Ord}(q_2))$$

*Proof.* Factor each $q_i$ into prime powers $q_i = \prod_j p_j^{a_{ij}}$. Then $\text{lcm}(q_1, q_2) = \prod_j p_j^{\max(a_{1j}, a_{2j})}$. By CRT, $\text{Ord}_q(2) = \text{lcm}_j(\text{Ord}_{p_j^{a_j}}(2))$ for any $q = \prod_j p_j^{a_j}$. Since each prime-power component $p_j^{\max(a_{1j}, a_{2j})}$ of $\text{lcm}(q_1, q_2)$ comes from whichever $q_i$ contributes the larger power of $p_j$, the lcm of the resulting orders equals the lcm of the orders from $q_1$ and $q_2$ individually. $\square$

This rule is pervasive: it governs the field degree needed for [root expression](Root%20Expressions%20and%20LC%20Estimation.md) computations, determines where products of roots from different CMPR blocks land, and underlies the period analysis of composite registers.

**Counting elements of a given period.** The number of elements of multiplicative order exactly $q$ in $\overline{\mathbb{F}_2}$ is $\phi(q)$, where $\phi$ is Euler's totient function. These $\phi(q)$ elements partition into $\phi(q) / \text{Ord}(q)$ Frobenius orbits, each of size $\text{Ord}(q)$ — so there are exactly $\phi(q)/\text{Ord}(q)$ irreducible polynomials over $\mathbb{F}_2$ whose roots have period $q$. For the primitive case $q = 2^n - 1$, this recovers the count $\phi(2^n - 1)/n$ of primitive polynomials of degree $n$.

**The Carmichael function and the multiplicative group.** The multiplicative group $\mathbb{F}_{2^n}^\times \cong \mathbb{Z}/(2^n - 1)\mathbb{Z}$ is cyclic, so every element's order divides $2^n - 1$. The Carmichael function $\lambda(m)$ gives the exponent of the group $(\mathbb{Z}/m\mathbb{Z})^\times$ — the largest order of any element. The group $(\mathbb{Z}/m\mathbb{Z})^\times$ is cyclic exactly when $m$ is 1, 2, 4, an odd prime power $p^k$, or twice an odd prime power $2p^k$, and in those cases its exponent $\lambda(m)$ equals its order $\phi(m)$; otherwise $\lambda(m) < \phi(m)$.

In the register context, $\lambda(2^n - 1)$ is less directly relevant (since $\mathbb{F}_{2^n}^\times$ is cyclic, every divisor of $2^n - 1$ occurs as an element order). The Carmichael function becomes important when working with $(\mathbb{Z}/q\mathbb{Z})^\times$ for composite $q$ — for instance, the decimation-based primitive polynomial enumeration in `NumberTheory.py` uses $\lambda(p^k)$ to identify primitive elements modulo prime powers and then lifts them via CRT. The key identity is:

$$\lambda(p^k) = \begin{cases} p^{k}/2 & \text{if } p^k \in \{2, 4\} \\ p^{k-2} & \text{if } p = 2, k \geq 3 \\ p^{k-1}(p-1) & \text{if } p \text{ is odd} \end{cases}$$

**Exponentiation (decimation).** Raising a root to a power $d$ maps $\alpha^{p/q} \mapsto \alpha^{dp/q}$ (reducing $dp \bmod q$). In sequence terms, this corresponds to **decimation**: sampling every $d$-th element of the output sequence. Decimation by $d$ permutes the roots within $\mathbb{F}_{2^n}$ (when $\gcd(d, 2^n - 1) = 1$) and sends one primitive polynomial to another. This is the mechanism behind `fast_decimate` and `generate_primitive_polynomials` in the library. For the full treatment — including the inverse operation (expansion), effects on linear complexity, root fans, and Jordan structure — see [Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md).

### What Each Form Reveals

**Exponential form** (the period $q$, $\text{Ord}(q)$, orbit sizes, and field membership do not depend on the hierarchy; the numerator $p$ of a given element does):
- **Period and field membership:** Immediately visible from the denominator $q$ and $\text{Ord}(q)$.
- **Multiplication:** Exponents add in $\mathbb{Z}_{(2)}/\mathbb{Z}$.
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

The [root expression](Root%20Expressions%20and%20LC%20Estimation.md) framework lives almost entirely in the exponential world — it tracks periods and coset sizes without ever materializing polynomial representatives. The [resolvent analysis](Resolvent%20Analysis.md), by contrast, works in the polynomial/power-series world. The primal/dual framework in [Polynomial Conventions](../conventions/Polynomial%20Conventions.md) uses both: the primal perspective describes a polynomial by its roots, which are named in exponential form, while the dual perspective describes it by the sequences it annihilates under convolution, which is a computation on coefficients.

The `fast_decimate` function in `NumberTheory.py` illustrates the alternation pattern in practice: it starts with a primitive polynomial (polynomial form), generates a sequence by running a register (vector form), decimates the sequence (an operation defined by exponential-form reasoning — decimation by $d$ maps $\alpha^{p/q} \mapsto \alpha^{dp/q}$), then applies Berlekamp-Massey to recover the result as a polynomial. Each step uses whichever representation makes that particular operation natural.

## Polynomial Types

Minimal, characteristic, and primitive polynomials each play a distinct role in register analysis. See [Polynomial Types](Polynomial%20Types.md) for their definitions, relationships, and how they connect to the concepts developed above.

## Connections to Other Documents

- **[Polynomial Conventions](../conventions/Polynomial%20Conventions.md):** Defines the primal/dual framework — how a polynomial relates to a sequence through its roots (primal) or through convolution (dual). The algebraic closure is where those roots live.
- **[Polynomial Types](Polynomial%20Types.md):** Minimal, characteristic, and primitive polynomials — their definitions, relationships, and role in register analysis.
- **[Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md):** Uses the subfield structure of $\overline{\mathbb{F}_2}$ to bound linear complexity. Root expressions track how many roots from each subfield $\mathbb{F}_{2^e}$ can appear in a signal.
- **[Monomial Profile Theory](Monomial%20Profile%20Theory.md):** Parallels root expressions for ANF degree. The full-coset degeneracy (§1 of that document) is a direct consequence of $\alpha^{2^e - 1} = 1$ in $\mathbb{F}_{2^e}$.
- **[Resolvent Analysis](Resolvent%20Analysis.md):** The resolvent matrix has entries in $\mathbb{F}_2[[D]]$, and its eigenstructure feeds the root expression computation.
- **[Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md):** Develops the root-level theory of decimation ($\alpha^{p/q} \mapsto \alpha^{dp/q}$) and expansion (root fans), their effects on LC and Jordan structure, and the duality between them.
- **[Algebraic Closure](Algebraic%20Closure.md):** A diagram of three quotients of $\mathbb{Z}[x]$ with their localizations and completions, read as sets of binary strings (eventually constant, eventually periodic, arbitrary).
