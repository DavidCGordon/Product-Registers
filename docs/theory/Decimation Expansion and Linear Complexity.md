# Decimation, Expansion, and Linear Complexity

When you [interleave](Interleaving%20and%20Clock-Controlled%20Registers.md) sequences or [shrink](Interleaving%20and%20Clock-Controlled%20Registers.md#shrinking) one by another, the root-level mechanics are **expansion** (interleaving introduces new roots via fans and Jordan blocks) and **decimation** (extracting a phase maps roots by multiplication). This document develops those mechanics: what each operation does to [roots](Roots%20and%20the%20Algebraic%20Closure.md), [Frobenius orbits](Roots%20and%20the%20Algebraic%20Closure.md), and [Jordan structure](Root%20Multiplicities%20and%20Jordan%20Decomposition.md), and how the effects compose.

## Decimation

**Decimation by $d$** extracts every $d$-th sample from a sequence:

$$s'[t] = s[dt].$$

This is what happens when you extract a single phase from a $d$-phase interleaved stream, or when the shrinking generator's phase sequences sample the data at stride $q_A$.

### Root Rule

A root $\alpha^{p/q}$ in the original maps to $\alpha^{dp \bmod q\, /\, q}$ in the decimated sequence — the exponent is multiplied by $d$:

$$\boxed{\alpha^{p/q} \;\longmapsto\; \alpha^{(dp \bmod q)/q}}$$

Whether this map is injective determines whether LC is preserved:

- **$\gcd(d, q) = 1$:** The map $p \mapsto dp$ is a bijection on $\mathbb{Z}/q\mathbb{Z}$. Roots permute, nothing collides. **LC and period are preserved.** This is the case in the shrinking generator when $\gcd(q_A, q_B) = 1$: the data roots $\rho \mapsto \rho^{q_A}$ stay distinct across phases.

- **$\gcd(d, q) = g > 1$:** Roots can collide — distinct $\alpha^{p_1/q}$ and $\alpha^{p_2/q}$ map to the same value when $p_1 \equiv p_2 \pmod{q/g}$. Over $\mathbb{F}_2$, even-multiplicity roots cancel. **LC can drop**, sometimes to $q/g$. The worst case is $d = q$: every root maps to $\alpha^0 = 1$.

### Frobenius Orbits Under Decimation

Frobenius (squaring) commutes with decimation: the orbit $\{p, 2p, 4p, \ldots\}$ maps to $\{dp, 2dp, 4dp, \ldots\}$. When $\gcd(d,q)=1$, each orbit maps bijectively to another of the same size.

### Coset Weight

Decimation does not preserve coset weight. The map $p \mapsto dp \bmod q$ scrambles binary representations, so decimation has no clean description in terms of [coset classes](Root%20Expressions%20and%20LC%20Estimation.md#coset-classes).

### Jordan Blocks Under Decimation

**By $2^k$:** A Jordan term $\binom{t}{r}\alpha^t$ becomes $\binom{2^k t}{r}(\alpha^{2^k})^t$. By Lucas' theorem, $\binom{2^k t}{r} \equiv 0 \pmod{2}$ for $1 \leq r < 2^k$. Jordan size compresses:

$$m' = \lceil m / 2^k \rceil$$

This exactly inverts even expansion.

| Original $m$ | $\div 2$ | $\div 4$ | $\div 8$ |
|---|---|---|---|
| 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 |
| 4 | 2 | 1 | 1 |
| 5 | 3 | 2 | 1 |
| 8 | 4 | 2 | 1 |

**By odd $d > 1$:** No closed form. The interaction of $\binom{dt}{r} \bmod 2$ with $d$, $t$, $r$ via Lucas' theorem is complex. Best determined empirically.

## Expansion

**Expansion by $d$** inserts $d-1$ zeros between each sample:

$$s''[t] = \begin{cases} s[t/d] & d \mid t \\ 0 & \text{otherwise} \end{cases}$$

In the D-transform: $S''(D) = S(D^d)$. This is the operation that underlies each term in the [polyphase identity](Interleaving%20and%20Clock-Controlled%20Registers.md#the-polyphase-identity): interleaving $d$ sequences is a sum of $d$ expansions with phase offsets.

The [odd-denominator constraint](Roots%20and%20the%20Algebraic%20Closure.md) — every element of $\overline{\mathbb{F}_2}$ has odd multiplicative order — forces even and odd factors to work through different mechanisms.

### Odd Expansion: Root Fans

When $d$ is odd, $dq$ is odd, and each root fans into $d$ new roots with valid labels:

$$\boxed{\alpha^{p/q} \;\longmapsto\; \left\{\, \alpha^{(p + kq)/(dq)} \;:\; k = 0, \ldots, d-1 \,\right\}}$$

Each fan element $\beta_k = \alpha^{(p+kq)/(dq)}$ satisfies $\beta_k^d = \alpha^{(p+kq)/q} = \alpha^{p/q}$ (since integer shifts are trivial in $\mathbb{Q}/\mathbb{Z}$). So the fan consists of the $d$-th roots of the original root.

**Fans are disjoint.** Distinct roots produce disjoint fans. *Proof:* if $(p_1 + k_1 q) \equiv (p_2 + k_2 q) \pmod{dq}$ with $|k_i| < d$ and $|p_1 - p_2| < q$, then $k_1 = k_2$ and $p_1 = p_2$. $\square$

**Frobenius closure.** A single fan need not be a complete Frobenius orbit — squaring can carry fan elements across fans. But the union of fans over a complete original orbit is always a union of complete orbits at the expanded level.

*Proof.* Let $U = \pi^{-1}(C_q(p))$ where $\pi: \mathbb{Z}/(dq)\mathbb{Z} \to \mathbb{Z}/q\mathbb{Z}$ is reduction mod $q$. Since $C_q(p)$ is closed under doubling and $\pi$ commutes with doubling, $U$ is closed under doubling. Since $dq$ is odd, doubling permutes $\mathbb{Z}/(dq)\mathbb{Z}$, so $U$ decomposes into complete orbits. $\square$

**LC multiplier:** expansion by odd $d$ multiplies LC by exactly $d$.

**Example.** Period-3 sequence (roots $\alpha^{1/3}, \alpha^{2/3}$, LC = 2), expanded by 5:

| Original | Fan roots (mod 15) | Orbits hit |
|---|---|---|
| $\alpha^{1/3}$ | $1, 4, 7, 10, 13$ | $\{1,2,4,8\}$, $\{7,14,13,11\}$, $\{5,10\}$ |
| $\alpha^{2/3}$ | $2, 5, 8, 11, 14$ | same three orbits |

LC $= 4 + 4 + 2 = 10 = 5 \times 2$. Note $\alpha^{10/15} = \alpha^{2/3}$ and $\alpha^{5/15} = \alpha^{1/3}$ — fan elements can reduce to the original roots.

### Even Expansion: Frobenius Squaring

When $d = 2$, the characteristic-2 identity $(f(D))^2 = f(D^2)$ gives $S(D^2) = (S(D))^2$. If $m(D)$ is the minimal polynomial:

$$m(D^2) = (m(D))^2$$

The roots are unchanged — squaring is a bijection on $\overline{\mathbb{F}_2}$. Instead, each **Jordan block size doubles**. By induction, expansion by $2^a$:

$$m(D^{2^a}) = (m(D))^{2^a}, \qquad \text{LC} = 2^a \cdot L(s), \qquad \text{Jordan: } m \mapsto 2^a m$$

**Example.** Period-3, $m(D) = 1+D+D^2$:

| $2^a$ | Minimal polynomial | Jordan | Period | LC |
|---|---|---|---|---|
| 1 | $1+D+D^2$ | 1 | 3 | 2 |
| 2 | $(1+D+D^2)^2$ | 2 | 6 | 4 |
| 4 | $(1+D+D^2)^4$ | 4 | 12 | 8 |
| 8 | $(1+D+D^2)^8$ | 8 | 24 | 16 |

### General Expansion: Odd–Even Decomposition

For $d = 2^a \cdot d'$ ($d'$ odd), the two mechanisms combine independently:

$$\boxed{L(s'') = d \cdot L(s) = \underbrace{d'}_{\text{fan}} \cdot \underbrace{2^a}_{\text{Jordan}} \cdot L(s)}$$

The odd part creates $d'$ new eigenvalues per root. The even part scales each Jordan block by $2^a$. Period: $d'q \cdot 2^a = dq$.

**Example.** Period-3, LC-2, expanded by $6 = 2 \cdot 3$: the odd part fans 2 roots into 6 (in the group of order 9), then the even part doubles each Jordan block. Result: 6 eigenvalues with Jordan size 2, period $9 \cdot 2 = 18$, LC $= 12$.

## The Duality

Decimation and expansion invert each other at the root level:

| | Decimation by $d$ | Expansion by $d = 2^a d'$ |
|---|---|---|
| **Operation** | $s'[t] = s[dt]$ | $s''[t] = s[t/d]$ if $d \mid t$, else 0 |
| **Odd part** | $\alpha^{p/q} \mapsto \alpha^{dp/q}$ (can collide) | $\alpha^{p/q} \mapsto d'$ fan roots (no overlap) |
| **Even part** | Jordan: $m \mapsto \lceil m/2^a \rceil$ | Jordan: $m \mapsto 2^a m$ |
| **LC** | $\leq L(s)$, equality iff $\gcd(d,q)=1$ | $= d \cdot L(s)$ always |
| **Period** | divides $q/\gcd(d,q)$ | exactly $dq$ |

**Asymmetry:** expansion always multiplies LC by $d$ (disjoint fans, no cancellation). Decimation preserves LC only when coprime.

**Round-trip:** expand then decimate recovers the original. Decimate then expand does not — it creates a zero-stuffed version with higher LC.

## Primitive Polynomial Enumeration

The library's `generate_primitive_polynomials` uses the coprime-decimation property. Decimating a maximal-length sequence by $d$ with $\gcd(d, 2^n - 1) = 1$ permutes roots, producing a different primitive polynomial of the same degree. Two decimations $d_1, d_2$ give the same polynomial when $d_1 \equiv 2^k d_2 \pmod{2^n-1}$ (Frobenius already generates the orbit). The function `fast_decimate` runs the register, decimates, and recovers the polynomial via Berlekamp-Massey.

## Connections

- **[Interleaving and Clock-Controlled Registers](Interleaving%20and%20Clock-Controlled%20Registers.md):** Uses expansion (root fans + Jordan) for interleaving LC analysis, and decimation for shrinking-generator phase roots.
- **[Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md):** Root labels $\alpha^{p/q}$, Frobenius orbits, odd-denominator constraint.
- **[Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md):** Jordan blocks, multiplicity brackets, period staircase.
- **[Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md):** Coset classes — expansion fans into larger groups, decimation scrambles weights.
- **[Coset Regions](../conventions/Coset%20Regions.md):** Grid picture applies to expansion but not decimation.
