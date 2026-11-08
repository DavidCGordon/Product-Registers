# Decimation, Expansion, and Linear Complexity

This document describes what two operations on a binary sequence do to its [roots](Roots%20and%20the%20Algebraic%20Closure.md): **decimation**, which keeps every $d$-th sample, and **expansion**, which inserts $d-1$ zeros between consecutive samples. For each operation it gives the effect on roots, on [Frobenius orbits](Roots%20and%20the%20Algebraic%20Closure.md), and on [Jordan structure](Root%20Multiplicities%20and%20Jordan%20Decomposition.md), and it shows how the effects compose. The two operations are the root-level mechanics of [interleaving](Interleaving%20and%20Clock-Controlled%20Registers.md) and [shrinking](Interleaving%20and%20Clock-Controlled%20Registers.md#shrinking): interleaving introduces new roots through expansion (root fans and larger Jordan blocks), and extracting one phase of an interleaved stream is a decimation, which multiplies each root's exponent by $d$.

The document assumes the exponential notation $\alpha^{p/q}$ for roots and the Frobenius orbit of a root (its orbit under repeated squaring), both from [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md). It also assumes Jordan blocks as developed in [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md): a root $\alpha$ carrying a Jordan block of size $m$ contributes the terms $\binom{t}{r}\alpha^t$, $0 \le r < m$, to the sequence. $L(s)$ denotes the linear complexity of a sequence $s$, abbreviated LC.

The letter $q$ plays two roles. In a root label $\alpha^{p/q}$, as in the root rule and the fan formula, $q$ is the denominator of that root. In statements about the period of a sequence, and in the LC and Period rows of the table in [The Duality](#the-duality), $q$ is the least period of the sequence $s$.

## Decimation

**Decimation by $d$** extracts every $d$-th sample from a sequence:

$$s'[t] = s[dt].$$

This is the operation that extracts a single phase from a $d$-phase interleaved stream. It is also how the phase sequences of the shrinking generator sample the data: at stride $q_A$, the period of the control sequence.

### Root Rule

A root $\alpha^{p/q}$ of the original sequence maps to $\alpha^{(dp \bmod q)/q}$ in the decimated sequence: the exponent is multiplied by $d$.

$$\boxed{\alpha^{p/q} \;\longmapsto\; \alpha^{(dp \bmod q)/q}}$$

Whether this map is injective governs what can happen to the LC:

- **$\gcd(d, q) = 1$:** The map $p \mapsto dp$ is a bijection on $\mathbb{Z}/q\mathbb{Z}$. The roots are permuted and none collide. If in addition every Jordan size is 1, **LC and period are preserved.** The same conclusion holds, without the Jordan condition, when $q$ is read as the least period of $s$, which is the form used in the table in [The Duality](#the-duality). With $q$ a root denominator the Jordan condition is needed: $q$ is odd, so $d = 2$ satisfies $\gcd(d, q) = 1$, yet decimation by 2 compresses Jordan blocks (the sequence $t \bmod 2$ has LC 2 and decimates to the zero sequence; see [Jordan Blocks Under Decimation](#jordan-blocks-under-decimation)). This is the case in the shrinking generator when $\gcd(q_A, q_B) = 1$, where $q_B$ is the period of the data sequence: the data roots $\rho \mapsto \rho^{q_A}$ stay distinct across phases.

- **$\gcd(d, q) = g > 1$:** Roots can collide. Distinct roots $\alpha^{p_1/q}$ and $\alpha^{p_2/q}$ map to the same value when $p_1 \equiv p_2 \pmod{q/g}$. When roots collide, their coefficients, which are elements of $\overline{\mathbb{F}_2}$, add, and the merged root is present in the decimated sequence if and only if the sum is nonzero. **LC can drop.** The worst case is $d = q$, where every root maps to $\alpha^0 = 1$.

Two period-3 sequences decimated by 3 show both outcomes of a collision. Each has the two roots $\alpha^{1/3}$ and $\alpha^{2/3}$, which both map to 1. The sequence $1,1,0,1,1,0,\ldots$ decimates to the all-ones sequence, of LC 1: the coefficients sum to a nonzero value and the merged root survives. The sequence $0,1,1,0,1,1,\ldots$ decimates to the zero sequence: the coefficients cancel.

When $\gcd(d, q) > 1$ the LC may drop, but it need not. The degree-4 m-sequence of period 15, decimated by 3, has period 5 and still has LC 4: its roots have numerators $\{1, 2, 4, 8\}$ over the denominator 15, and these are distinct modulo $15/3 = 5$, so no two roots collide.

The drop is bounded as follows. If $s$ has least period $T$, then

$$L(s') \leq \min\bigl(L(s),\; T/\gcd(d, T)\bigr),$$

The second bound holds because the period of $s'$ divides $T/\gcd(d,T)$ and a sequence of period $N$ has LC at most $N$. For the first, let $m_\rho$ be the Jordan size at each root $\rho$ of $s$, so that $L(s) = \sum_\rho m_\rho$. Decimation turns each term $\binom{t}{r}\rho^t$ of $s$ into $\binom{dt}{r}(\rho^d)^t$, and by the [Binomial Decimation Identity](Binomial%20Decimation%20Identity.md) $\binom{dt}{r}$ is, modulo 2, a combination of the $\binom{t}{i}$ with $i \leq r$. So $s'$ lies in the span of the sequences $\binom{t}{i}(\rho^d)^t$ with $i < m_\rho$, a space of dimension at most $\sum_\rho m_\rho = L(s)$.

### Frobenius Orbits Under Decimation

Frobenius (squaring) commutes with decimation: the orbit $\{p, 2p, 4p, \ldots\}$ maps to $\{dp, 2dp, 4dp, \ldots\}$. When $\gcd(d,q)=1$, each orbit maps bijectively to another of the same size.

### Coset Weight

In $\mathbb{F}_{2^e}$, write a root as the $k$-th power of the primitive element, with the exponent $k$ taken modulo $2^e - 1$; for a root $\alpha^{p/q}$ this exponent is $k = p(2^e-1)/q$. Squaring the root doubles $k$ modulo $2^e - 1$, which is a cyclic shift of the $e$-bit binary string of $k$. The exponents in one Frobenius orbit therefore share a Hamming weight, called the **coset weight** (see [coset classes](Root%20Expressions%20and%20LC%20Estimation.md#coset-classes)). The statement is about the exponents $k$, not the numerators $p$: for $q = 5$ and $e = 4$ the orbit of numerators $\{1, 2, 4, 3\}$ has binary weights $1, 1, 1, 2$.

Decimation does not preserve coset weight. For $e = 4$ and $d = 7$, the exponent $k = 1$ (weight 1) maps to $7$ (weight 3). Decimation therefore has no clean description in terms of coset classes.

### Jordan Blocks Under Decimation

**Decimation by $2^k$.** A Jordan term $\binom{t}{r}\alpha^t$ becomes $\binom{2^k t}{r}(\alpha^{2^k})^t$. Lucas' theorem writes $\binom{2^k t}{r} \bmod 2$ as a product, over bit positions, of the binomial coefficients of the corresponding bits of $2^k t$ and $r$. The lowest $k$ bits of $2^k t$ are zero, so the product is zero unless the lowest $k$ bits of $r$ are also zero, that is, unless $2^k \mid r$. When $r = 2^k r'$, the remaining bit positions are those of $t$ and $r'$, so

$$\binom{2^k t}{r} \equiv \begin{cases} \binom{t}{r'} & r = 2^k r' \\ 0 & 2^k \nmid r \end{cases} \pmod{2}.$$

Write the part of the sequence at root $\alpha$ as $\sum_{r=0}^{m-1} a_r \binom{t}{r}\alpha^t$. After decimation only the terms with $r = 2^k r'$ remain, and they become $a_{2^k r'}\binom{t}{r'}(\alpha^{2^k})^t$. The surviving indices are $0 \leq r' \leq \lfloor (m-1)/2^k \rfloor$, of which there are $\lceil m/2^k \rceil$. The Jordan size after decimation is therefore

$$m' = 1 + \max\{\, r' : a_{2^k r'} \neq 0 \,\} \;\leq\; \lceil m / 2^k \rceil,$$

with equality if and only if $a_{2^k \lfloor (m-1)/2^k \rfloor} \neq 0$. If every $a_{2^k r'}$ is zero, the root disappears from the decimated sequence.

The table gives the bound $\lceil m/2^k \rceil$ for several block sizes.

| Original $m$ | Bound for $\div 2$ | Bound for $\div 4$ | Bound for $\div 8$ |
|---|---|---|---|
| 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 |
| 4 | 2 | 1 | 1 |
| 5 | 3 | 2 | 1 |
| 8 | 4 | 2 | 1 |

The bound is not always attained. The sequence $s[t] = t \bmod 2 = \binom{t}{1}$ has Jordan size 2 at the root 1, with $a_0 = 0$ and $a_1 = 1$. The table entry for $m = 2$ decimated by 2 is 1, but $s[2t] = 0$: the root disappears.

Decimation by $2^k$ is a left inverse of the Jordan scaling $m \mapsto 2^k m$ produced by expansion by $2^k$, which is derived in [Even Expansion](#even-expansion-frobenius-squaring) below: $\lceil 2^k m / 2^k \rceil = m$. It is not a right inverse. A block of size 5 decimated by 2 has size at most 3, and expansion by 2 then gives at most 6, and in any case an even size, so not 5.

**Decimation by odd $d$, coprime to the period.** When $d$ is odd and coprime to the period, so that no roots collide, the Jordan size is preserved: $m' = m$. The [Binomial Decimation Identity](Binomial%20Decimation%20Identity.md) (Step 4 of its proof) gives integers $c_i$ with

$$\binom{dt}{r} = \sum_{i=0}^{r} c_i \binom{t}{i}$$

as polynomials in $t$. Both sides have degree $r$ in $t$. The leading coefficient on the left is $d^r/r!$ and on the right is $c_r/r!$, so $c_r = d^r$, which is odd when $d$ is odd. Modulo 2, then, $\binom{dt}{r}$ is $\binom{t}{r}$ plus terms $\binom{t}{i}$ with $i < r$. Applied to a block of size $m$, the top term $\binom{t}{m-1}$ survives with its coefficient and no term above it appears.

> **OPEN:** What is the Jordan size after decimation by an odd $d$ with $\gcd(d, q) > 1$, when blocks at several colliding roots merge? — needs author input.
> The derivation above covers only the case with no collisions. For colliding roots, the collision rule for coefficients in the Root Rule section applies term by term, but no formula for the resulting Jordan size is established.

## Expansion

**Expansion by $d$** inserts $d-1$ zeros between consecutive samples:

$$s''[t] = \begin{cases} s[t/d] & d \mid t \\ 0 & \text{otherwise} \end{cases}$$

In the D-transform $S(D) = \sum_t s[t] D^t$, expansion is the substitution $S''(D) = S(D^d)$. This is the operation that underlies each term in the [polyphase identity](Interleaving%20and%20Clock-Controlled%20Registers.md#the-polyphase-identity): interleaving $d$ sequences is a sum of $d$ expansions with phase offsets.

**LC multiplier.** If $s$ is purely periodic, expansion by any $d \geq 1$ multiplies LC by exactly $d$: $L(s'') = d \cdot L(s)$. The hypothesis is needed: $s = 1, 0, 0, \ldots$ has $S(D) = 1$ and LC 1, and its expansion is the same sequence.

*Proof.* The proof uses the fact that the LC of a purely periodic sequence equals the degree of the denominator of its transform in lowest terms. Write $S(D) = N(D)/M(D)$ in lowest terms, with $\deg N < \deg M = L(s)$. The expanded sequence has transform $N(D^d)/M(D^d)$. Because $N$ and $M$ are coprime, there are polynomials $a, b$ with $aN + bM = 1$; substituting $D^d$ for $D$ gives $a(D^d)N(D^d) + b(D^d)M(D^d) = 1$, so $N(D^d)/M(D^d)$ is still in lowest terms. Its denominator has degree $\deg M(D^d) = d \cdot L(s)$. $\square$

The rest of this section describes where the factor $d$ comes from at the level of roots. The [odd-denominator constraint](Roots%20and%20the%20Algebraic%20Closure.md) — every element of $\overline{\mathbb{F}_2}$ has odd multiplicative order — forces the even and odd factors of $d$ to work through different mechanisms. The next two subsections treat odd $d$ and $d = 2$; the third combines them.

### Odd Expansion: Root Fans

When $d$ is odd, $dq$ is odd, and each root fans into $d$ new roots with valid labels:

$$\boxed{\alpha^{p/q} \;\longmapsto\; \left\{\, \alpha^{(p + kq)/(dq)} \;:\; k = 0, \ldots, d-1 \,\right\}}$$

Each fan element $\beta_k = \alpha^{(p+kq)/(dq)}$ satisfies $\beta_k^d = \alpha^{(p+kq)/q} = \alpha^{p/q}$ (since integer shifts are trivial in $\mathbb{Q}/\mathbb{Z}$). So the fan consists of the $d$-th roots of the original root.

**Fans are disjoint.** Distinct roots produce disjoint fans. *Proof:* every element $\beta$ of the fan of $\alpha^{p/q}$ satisfies $\beta^d = \alpha^{p/q}$, so a fan element determines the root it came from. An element cannot lie in the fans of two distinct roots. $\square$

**Frobenius closure.** A single fan need not be a complete Frobenius orbit — squaring can carry fan elements across fans. But the union of fans over a complete original orbit is always a union of complete orbits at the expanded level.

*Proof.* Write $C_q(p) = \{p, 2p, 4p, \ldots\} \subseteq \mathbb{Z}/q\mathbb{Z}$ for the Frobenius orbit of $p$, as a set of exponent numerators modulo $q$. Let $U = \pi^{-1}(C_q(p))$ where $\pi: \mathbb{Z}/(dq)\mathbb{Z} \to \mathbb{Z}/q\mathbb{Z}$ is reduction mod $q$. Since $C_q(p)$ is closed under doubling and $\pi$ commutes with doubling, $U$ is closed under doubling. Since $dq$ is odd, doubling permutes $\mathbb{Z}/(dq)\mathbb{Z}$, so $U$ decomposes into complete orbits. $\square$

**Example.** Take a period-3 sequence with roots $\alpha^{1/3}, \alpha^{2/3}$ (LC 2) and expand it by 5. The fan roots have denominator 15; the table lists their exponent numerators and the Frobenius orbits modulo 15 that they fall in.

| Original | Fan roots (mod 15) | Orbits hit |
|---|---|---|
| $\alpha^{1/3}$ | $1, 4, 7, 10, 13$ | $\{1,2,4,8\}$, $\{7,14,13,11\}$, $\{5,10\}$ |
| $\alpha^{2/3}$ | $2, 5, 8, 11, 14$ | same three orbits |

The expanded sequence has LC $= 4 + 4 + 2 = 10 = 5 \times 2$. Fan elements can reduce to the original roots: here $\alpha^{10/15} = \alpha^{2/3}$ and $\alpha^{5/15} = \alpha^{1/3}$.

### Even Expansion: Frobenius Squaring

When $d = 2$, the characteristic-2 identity $(f(D))^2 = f(D^2)$ gives $S(D^2) = (S(D))^2$. If $M(D)$ is the minimal polynomial:

$$M(D^2) = (M(D))^2$$

The roots are unchanged — squaring is a bijection on $\overline{\mathbb{F}_2}$. Instead, each **Jordan block size doubles**. By induction, expansion by $2^a$ gives

$$M(D^{2^a}) = (M(D))^{2^a}, \qquad \text{LC} = 2^a \cdot L(s), \qquad \text{Jordan: } m \mapsto 2^a m$$

**Example.** For a period-3 sequence with $M(D) = 1+D+D^2$:

| $2^a$ | Minimal polynomial | Jordan | Period | LC |
|---|---|---|---|---|
| 1 | $1+D+D^2$ | 1 | 3 | 2 |
| 2 | $(1+D+D^2)^2$ | 2 | 6 | 4 |
| 4 | $(1+D+D^2)^4$ | 4 | 12 | 8 |
| 8 | $(1+D+D^2)^8$ | 8 | 24 | 16 |

### General Expansion: Odd–Even Decomposition

For $d = 2^a \cdot d'$ ($d'$ odd), the two mechanisms combine independently:

$$\boxed{L(s'') = d \cdot L(s) = \underbrace{d'}_{\text{fan}} \cdot \underbrace{2^a}_{\text{Jordan}} \cdot L(s)}$$

The odd part creates $d'$ new eigenvalues per root. The even part scales each Jordan block by $2^a$. If $s$ has least period $q$, the period of $s''$ is $d'q \cdot 2^a = dq$.

**Example.** Expand a period-3, LC-2 sequence by $6 = 2 \cdot 3$. The odd part fans 2 roots into 6 (in the group of order 9), then the even part doubles each Jordan block. The result has 6 eigenvalues with Jordan size 2, period $9 \cdot 2 = 18$, and LC $= 12$.

## The Duality

Decimation by $d$ is a left inverse of expansion by $d$: expanding and then decimating recovers the original sequence, but decimating and then expanding does not. The table sets the results of the two preceding sections side by side. In the Odd part row $q$ is a root's denominator; in the LC and Period rows it is the least period of $s$.

| | Decimation by $d$ | Expansion by $d = 2^a d'$ |
|---|---|---|
| **Operation** | $s'[t] = s[dt]$ | $s''[t] = s[t/d]$ if $d \mid t$, else 0 |
| **Odd part** | $\alpha^{p/q} \mapsto \alpha^{dp/q}$ (can collide) | $\alpha^{p/q} \mapsto d'$ fan roots (no overlap) |
| **Even part** | Jordan: $m \mapsto$ at most $\lceil m/2^a \rceil$ | Jordan: $m \mapsto 2^a m$ |
| **LC** | $\leq L(s)$, with equality if $\gcd(d,q)=1$ | $= d \cdot L(s)$ for purely periodic $s$ |
| **Period** | divides $q/\gcd(d,q)$ | exactly $dq$ |

**Asymmetry.** Expansion of a purely periodic sequence always multiplies LC by $d$, as proved in [Expansion](#expansion). The cases quoted in this document (odd decimation coprime to the period, the period-15 m-sequence decimated by 3, the two period-3 collisions, $t \bmod 2$ decimated by 2, and expansion by $d \le 5$) are asserted in `test/Tools/test_decimation_expansion.py`. Decimation preserves LC when $\gcd(d,q)=1$; when $\gcd(d,q) > 1$ the LC may drop.

**Round trip.** Expanding and then decimating recovers the original sequence. Decimating and then expanding does not: it creates a zero-stuffed version with higher LC.

## Primitive Polynomial Enumeration

The library's `generate_primitive_polynomials` (in `PyPR.Tools.AlgClosure.NumberTheory`) uses the coprime-decimation property. Decimating a maximal-length sequence by $d$ with $\gcd(d, 2^n - 1) = 1$ permutes roots, producing a different primitive polynomial of the same degree. Two decimations $d_1, d_2$ give the same polynomial when $d_1 \equiv 2^k d_2 \pmod{2^n-1}$ (Frobenius already generates the orbit). The function `fast_decimate` runs the register, decimates, and recovers the polynomial via Berlekamp-Massey.

## Connections

- **[Interleaving and Clock-Controlled Registers](Interleaving%20and%20Clock-Controlled%20Registers.md)** uses expansion (root fans and Jordan scaling) for interleaving LC analysis, and decimation for shrinking-generator phase roots.
- **[Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md)** defines root labels $\alpha^{p/q}$, Frobenius orbits, and the odd-denominator constraint.
- **[Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md)** develops Jordan blocks, multiplicity brackets, and the period staircase.
- **[Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md)** defines coset classes. Expansion fans into larger groups; decimation does not preserve coset weight.
- **[Coset Regions](../conventions/Coset%20Regions.md)** gives a grid picture that applies to expansion but not to decimation.
- **[Binomial Decimation Identity](Binomial%20Decimation%20Identity.md)** decomposes $\binom{dn}{k}$ in the binomial basis; its top coefficient $d^k$ gives the preservation of Jordan size under odd decimation derived above.
