# Interleaving and Clock-Controlled Registers

Two operations generate every clock-controlled construction: **interleaving** (merging several sequences into one) and **shrinking** (keeping a subsequence selected by a control). This document develops both, then derives the LC formulas for the major circuit types by analyzing what these operations do to [roots](Roots%20and%20the%20Algebraic%20Closure.md).

Prerequisites: [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md) (root labels $\alpha^{p/q}$, Frobenius orbits), [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md) (Jordan blocks, multiplicity brackets), [Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md) (root fans, Jordan scaling — the root-level mechanics behind interleaving).

## Interleaving

### Definition

Given $d$ sequences $s_0, \ldots, s_{d-1}$, their **$d$-phase interleaving** writes them into alternating positions of a single output:

$$y[dt + r] = s_r[t], \qquad 0 \leq r < d.$$

The output cycles through the $d$ phases: $s_0[0], s_1[0], \ldots, s_{d-1}[0], s_0[1], s_1[1], \ldots$

Conversely, extracting phase $r$ recovers the source: $s_r[t] = y[dt + r]$.

### The Polyphase Identity

In the D-transform ($S(D) = \sum_t s[t] D^t$):

$$\boxed{Y(D) = \sum_{r=0}^{d-1} D^r\, S_r(D^d)}$$

Each $S_r(D^d)$ is the [expansion](Decimation%20Expansion%20and%20Linear%20Complexity.md#definition-1) of phase $r$ by $d$ — samples placed at every $d$-th position, zeros elsewhere. The shift $D^r$ slides it into the correct slot. Summing fills all positions.

This identity says: **interleaving = sum of expansions with phase offsets**. Every construction below reduces to specifying what the phase sequences $s_r$ are.

### What Interleaving Does to Roots

Since interleaving is a sum of expansions, the [root fan and Jordan scaling](Decimation%20Expansion%20and%20Linear%20Complexity.md#general-expansion-oddeven-decomposition) from expansion apply phase by phase. Write $d = 2^a d'$ ($d'$ odd):

- **Odd part** ($d'$): each phase root $\alpha^{p/q}$ fans into $d'$ new roots $\alpha^{(p+kq)/(d'q)}$.
- **Even part** ($2^a$): each Jordan block scales by $2^a$ (no new roots).

The candidate root set of the output is the union across all phases:

$$R_Y \subseteq \bigcup_{r=0}^{d-1} \bigcup_{\rho \in R_r} F_{d'}(\rho)$$

This is an upper bound — roots from different phases can cancel in the sum. The shift $D^r$ changes coefficients but does not introduce new roots.

**Frobenius orbits are preserved.** The fan of a complete Frobenius orbit is a union of complete orbits at the expanded level ([proof](Decimation%20Expansion%20and%20Linear%20Complexity.md#frobenius-closure-of-fan-unions)), so the LC is always a sum of orbit sizes.

### Phase Cancellation

When phases share roots, the fan union overcounts. Shared roots can cancel in $Y(D) = \sum_r D^r S_r(D^d)$.

**Identical phases.** If all phases carry the same $s$:

$$Y(D) = \left(\sum_{r=0}^{d-1} D^r\right) S(D^d) = P(D) \cdot S(D^d)$$

where $P(D) = 1 + D + \cdots + D^{d-1}$ is the phase polynomial. Its roots can cancel with roots of $S(D^d)$.

**Two-phase, $d = 2$.** $Y(D) = A(D^2) + D \cdot B(D^2)$. If $a = b$, this is $(1+D) A(D^2)$: the factor $(1+D)$ adds root $\alpha^0 = 1$, and $A(D^2) = (A(D))^2$ doubles all Jordan blocks. So $\text{LC} = 2L(a) + 1$.

### LC Bound Procedure

Given interleaving rate $d = 2^a d'$:

1. **Phase roots.** Find root sets $R_r$ and Jordan sizes for each phase $s_r$.
2. **Fan by odd part.** Each root fans into $d'$ roots: $R_r^{\text{fan}} = \bigcup_{\rho \in R_r} F_{d'}(\rho)$.
3. **Scale by even part.** Multiply all Jordan sizes by $2^a$.
4. **Union.** $R_Y = \bigcup_r R_r^{\text{fan}}$.
5. **Cancel (optional).** Remove roots that cancel pairwise across phases.
6. **Count.** Partition into Frobenius orbits: $\text{LC} \leq \sum_{\text{orbits}} |\text{orbit}| \cdot (\text{max Jordan size})$.
7. **Period.** $\text{lcm}$ over surviving roots of $q \cdot 2^{\lceil \log_2 m \rceil}$.

### Worked Example

Interleave two distinct period-3 sequences ($R_0 = R_1 = \{\alpha^{1/3}, \alpha^{2/3}\}$, Jordan size 1) at rate $d = 2 = 2^1 \cdot 1$:

1. Both phases: roots $\{\alpha^{1/3}, \alpha^{2/3}\}$, size 1.
2. Odd $d' = 1$: no fanning.
3. Even $2^1$: Jordan sizes double to 2.
4. $R_Y = \{\alpha^{1/3}, \alpha^{2/3}\}$, each with block size 2.
5. Phases share roots. The minimal polynomial divides $(1+D+D^2)^2 \cdot (1+D)$.
6. One orbit of size 2 with Jordan length 2, plus root $\alpha^0$ from $(1+D)$. LC $\leq 2 \cdot 2 + 1 = 5$.
7. Period: $\text{lcm}(3 \cdot 2, 2) = 6$.

## Shrinking

### Definition

Given a **control** sequence $a[t]$ and a **data** sequence $b[t]$, the **shrinking** of $b$ by $a$ keeps only the data samples where the control is 1:

$$y[j] = b[\tau(j)], \qquad \tau(j) = j\text{-th time } t \text{ with } a[t] = 1.$$

The output has variable rate — its length per control period equals $w$, the number of ones in one period of $a$.

### Reduction to Interleaving

Group the output by control periods. If $a$ has period $q_A$ with $w$ ones at positions $t_0 < \cdots < t_{w-1}$, then in the $k$-th control period:

$$y[kw + j] = b[t_j + kq_A], \qquad 0 \leq j < w.$$

Define the **phase-$j$ slow-time sequence** $s_j[k] = b[t_j + kq_A]$. The shrunk output is the $w$-phase interleaving of $s_0, \ldots, s_{w-1}$.

### Phase Roots

Each phase $s_j$ samples $b$ at a fixed offset $t_j$ within each control period, advancing by $q_A$ steps per period. If $b[t] = \sum_\rho c_\rho \rho^t$:

$$s_j[k] = \sum_{\rho \in R_B} \underbrace{c_\rho \rho^{t_j}}_{\text{coefficient}} \cdot \underbrace{(\rho^{q_A})^k}_{\text{slow-time root}}$$

The slow-time roots are $\{\rho^{q_A} : \rho \in R_B\}$ — the data roots [decimated](Decimation%20Expansion%20and%20Linear%20Complexity.md#how-decimation-transforms-roots) by $q_A$. When $\gcd(q_A, q_B) = 1$, this decimation is a bijection (distinct roots stay distinct), so each phase has $L_B$ roots.

The interleaving of these $w$ phases then fans/scales each root according to the odd–even decomposition of $w$, giving the LC formula.

## Circuit Types

### Clock-Controlled Registers (Periodic Schedule)

A register with $d$ transition maps $U_0, \ldots, U_{d-1}$ applied in a repeating cycle. Augment the state to include the phase counter: $\widetilde{V} = V \times \mathbb{Z}/d\mathbb{Z}$. After one full cycle, the state advances by the **monodromy** $M = U_{d-1} \circ \cdots \circ U_0$.

The output at fast time $t = dk + r$ is $y[dk+r] = h(U_{r-1} \circ \cdots \circ U_0(M^k(v_0)))$. This is exactly a $d$-phase interleaving of slow-time sequences governed by $M$.

**Uniform clock** ($U_r = U$ for all $r$): $M = U^d$, and reconstructing the fast-time output from the $d$ phases is equivalent to expansion — root fans (odd $d$) or Jordan scaling (even $d$).

**Non-uniform clock**: $M$ is a product of distinct operators. The eigenstructure of $M$ must be computed directly.

### Controlled Sampling vs. Controlled Clocking

- **Controlled sampling** (shrinking): the register advances uniformly, but the output is recorded only at selected times. A decimation/interleaving problem.
- **Controlled clocking**: the transition itself changes per step. Requires the monodromy analysis.

## LC Formulas for Classical Constructions

### Shrinking Generator (SG)

**Setup.** Control LFSR $A$ (period $q_A$, $w$ ones per period). Data LFSR $B$ (period $q_B$, LC $= L_B$). Independent.

**LC formula.** From the reduction above — shrinking produces a $w$-phase interleaving of decimated data — with $\gcd(q_A, q_B) = 1$:

$$\boxed{L(\text{SG}) = w \cdot L_B \qquad (\gcd(q_A, q_B) = 1)}$$

This is an **equality**, not a bound. The $w$ phases have coefficients $c_\rho \rho^{t_j}$ where the positions $t_j$ are all distinct, so no pair of phases contributes identical terms that would cancel in the XOR sum. Verified exhaustively for small cases.

**The coprimality condition is essential.** When $\gcd(q_A, q_B) \neq 1$, the decimation $\rho \mapsto \rho^{q_A}$ collapses distinct roots, and LC drops below $w \cdot L_B$.

With the odd–even decomposition $w = 2^a w'$ ($w'$ odd): $w'$ fans roots (new eigenvalues), $2^a$ scales Jordan blocks.

**Example.** m-sequence control of degree $n_A$: period $q_A = 2^{n_A}-1$, $w = 2^{n_A-1}$. With data of degree $n_B$ and $\gcd(n_A, n_B) = 1$ (ensures $\gcd(q_A, q_B) = 1$):

$$L(\text{SG}) = 2^{n_A-1} \cdot n_B$$

| $n_A$ | $n_B$ | $w$ | LC |
|---|---|---|---|
| 3 | 4 | 4 | 16 |
| 3 | 5 | 4 | 20 |
| 4 | 3 | 8 | 24 |
| 4 | 5 | 8 | 40 |

### Self-Shrinking Generator (SSG)

The SSG shrinks a sequence by itself: pair consecutive bits, use the first as control for the second.

$$c[t] = a[2t], \qquad x[t] = a[2t+1], \qquad y[j] = x[\tau(j)] \text{ where } c[\tau(j)] = 1.$$

Both $c$ and $x$ are [decimations](Decimation%20Expansion%20and%20Linear%20Complexity.md#how-decimation-transforms-roots) of $a$ by 2 with different phase offsets. Since $q = 2^n - 1$ is odd, $\gcd(2, q) = 1$, so both have period $q$ and LC $n$. The control has $w = 2^{n-1}$ ones (balance property).

The SG formula gives $L \leq 2^{n-1} \cdot n$, but this overestimates because the control and data are **correlated** — both derived from one LFSR. The $w$ phase sequences are $w$ correlated views of the same $n$-dimensional state, not $w$ independent observations.

The known bound (Meier–Staffelbach 1994):

$$L(\text{SSG}) \leq 2^{n-1}$$

a factor of $n$ tighter. The SSG output is a nonlinear filter of the LFSR state, and the tighter bound follows from the filter-function degree bound.

Empirically, the $2^{n-1}$ bound is valid but not tight.

### Alternating-Step Generator (ASG)

**Setup.** Control LFSR $A$ (period $q_A$, $w$ ones). Data LFSRs $B$ (LC $= L_B$) and $C$ (LC $= L_C$). When $A[t] = 1$, advance $B$; when $A[t] = 0$, advance $C$. Output:

$$y[t] = B[N_1(t)] \oplus C[N_0(t)]$$

where $N_1(t)$ counts ones and $N_0(t)$ counts zeros of $A$ up to time $t$.

**Why this is not two shrinking generators.** The SG keeps only the $w$ accepted samples; the ASG outputs at every step. At positions where $A = 0$, register $B$ is not absent — it **holds its last value**. This sample-and-hold effect means $B$'s contribution spans all $q_A$ output positions, not just the $w$ where $B$ advances. The held values make $b[t] = B[N_1(t)]$ a nonlinear function of the combined (control, data) state, not a simple interleaving.

**Bound** (Günther 1988):

$$\boxed{L(\text{ASG}) \leq q_A \cdot (L_B + L_C)}$$

For m-sequence control ($q_A = 2^{n_A} - 1$): $L(\text{ASG}) \leq (2^{n_A}-1)(L_B + L_C)$.

Verified empirically: $n_A = 3$, $n_B = 4$, $n_C = 5$ gives LC $= 63 = 7 \times 9$, matching the bound exactly.

## Connections

- **[Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md):** The root-level mechanics behind interleaving (expansion = root fans + Jordan scaling) and phase extraction (decimation = root mapping).
- **[Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md):** Root labels, Frobenius orbits, odd-denominator constraint.
- **[Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md):** Jordan blocks and multiplicity brackets for even interleaving rates.
- **[Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md):** The XOR combination rule used in phase cancellation and ASG decomposition.
- **[Resolvent Analysis](Resolvent%20Analysis.md):** For clock-controlled registers with chaining.
- **[Coset Regions](../conventions/Coset%20Regions.md):** Geometric vocabulary for cross-field phase contributions.
