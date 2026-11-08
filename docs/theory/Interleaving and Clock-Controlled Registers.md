# Interleaving and Clock-Controlled Registers

A clock-controlled construction is one in which a control sequence decides either which outputs of a register are recorded or which transition the register takes at each step. This document analyses such constructions through two operations on sequences: **interleaving**, which merges several sequences into one, and **shrinking**, which keeps the subsequence of one sequence that a control selects. It develops both, shows that shrinking reduces to interleaving, and then states the linear complexity formulas and bounds for three classical generators — shrinking, self-shrinking and alternating-step — in terms of what the two operations do to [roots](Roots%20and%20the%20Algebraic%20Closure.md).

Prerequisites: [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md) (root labels $\alpha^{p/q}$, Frobenius orbits), [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md) (Jordan blocks, multiplicity brackets), and [Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md) (root fans, Jordan scaling — the root-level mechanics behind interleaving). The definitions this document relies on are restated here so that it can be read on its own.

All sequences are binary. The **linear complexity** (LC) of a sequence $s$, written $L(s)$, is the length of the shortest LFSR that generates it. A sequence satisfying a linear recurrence over $\mathbb{F}_2$ is a linear combination of terms $\binom{t}{r}\rho^t$ with $\rho$ in the algebraic closure $\overline{\mathbb{F}_2}$; the elements $\rho$ that occur are the **roots** of the sequence, and $R_s$ denotes the set of them. A root is written $\alpha^{p/q}$, where, for $p/q$ in lowest terms, the odd integer $q$ is its multiplicative order. The **Jordan size** of a root is the least $m$ such that the terms $\binom{t}{r}\rho^t$ occurring have $0 \leq r < m$; a root of Jordan size $1$ contributes only $\rho^t$. The **Frobenius orbit** of $\rho$ is $\{\rho, \rho^2, \rho^4, \ldots\}$, and in a sequence over $\mathbb{F}_2$ the roots of an orbit occur together.

## Interleaving

### Definition

Given $d$ sequences $s_0, \ldots, s_{d-1}$, their **$d$-phase interleaving** writes them into alternating positions of a single output:

$$y[dt + r] = s_r[t], \qquad 0 \leq r < d.$$

The output cycles through the $d$ phases: $s_0[0], s_1[0], \ldots, s_{d-1}[0], s_0[1], s_1[1], \ldots$ Conversely, extracting phase $r$ recovers the source: $s_r[t] = y[dt + r]$.

### The Polyphase Identity

The **D-transform** of a sequence $s$ is the power series $S(D) = \sum_t s[t] D^t$. In terms of it, interleaving reads

$$\boxed{Y(D) = \sum_{r=0}^{d-1} D^r\, S_r(D^d)}$$

Each $S_r(D^d)$ is the [expansion](Decimation%20Expansion%20and%20Linear%20Complexity.md#expansion) of phase $r$ by $d$: the samples of $s_r$ placed at every $d$-th position, with zeros elsewhere. The factor $D^r$ shifts that expansion into slot $r$, and the sum over $r$ fills every position. The identity therefore says that **interleaving is a sum of expansions with phase offsets**, and every construction below is analysed by specifying what its phase sequences $s_r$ are.

### What Interleaving Does to Roots

Since interleaving is a sum of expansions, the [root fan and Jordan scaling](Decimation%20Expansion%20and%20Linear%20Complexity.md#general-expansion-oddeven-decomposition) rules for expansion apply phase by phase. Write $d = 2^a d'$ with $d'$ odd. The two factors act differently.

The odd part $d'$ fans each phase root $\rho = \alpha^{p/q}$ into the $d'$ roots

$$F_{d'}(\rho) = \left\{\, \alpha^{(p+kq)/(d'q)} \;:\; k = 0, \ldots, d'-1 \,\right\},$$

the **fan** of $\rho$, which consists of the $d'$-th roots of $\rho$. The even part $2^a$ introduces no new roots; it multiplies the Jordan size of every root by $2^a$, before any cancellation between phases.

The candidate root set of the output is the union of the fans across all phases:

$$R_Y \subseteq \bigcup_{r=0}^{d-1} \bigcup_{\rho \in R_r} F_{d'}(\rho)$$

where $R_r$ is the root set of phase $s_r$. The inclusion can be strict, because roots contributed by different phases can cancel in the sum. The shift $D^r$ changes the coefficients with which the roots appear but does not introduce new roots.

The roots of the output come in whole Frobenius orbits. The output is a sequence over $\mathbb{F}_2$, so its minimal polynomial has coefficients in $\mathbb{F}_2$, and the roots of such a polynomial are closed under squaring. The candidate set is consistent with this: the fan of a complete Frobenius orbit is a union of complete orbits at the expanded level ([proof](Decimation%20Expansion%20and%20Linear%20Complexity.md#odd-expansion-root-fans)). The LC is therefore a sum over orbits, $\sum_{O} |O| \cdot m_O$, where $m_O$ is the Jordan size on the orbit $O$; when every Jordan size is $1$ it is a sum of orbit sizes.

### Phase Cancellation

When phases share roots, the union of fans overcounts: shared roots can cancel in $Y(D) = \sum_r D^r S_r(D^d)$. Two cases show the effect.

**Identical phases.** If all phases carry the same sequence $s$, the polyphase identity factors:

$$Y(D) = \left(\sum_{r=0}^{d-1} D^r\right) S(D^d) = P(D) \cdot S(D^d)$$

where $P(D) = 1 + D + \cdots + D^{d-1}$ is the **phase polynomial**. The roots of a sequence are the poles of its D-transform, and $P(D)$ multiplies the numerator of $Y(D)$. Its zeros can therefore cancel roots of $S(D^d)$, and they never add any. For a purely periodic $s$, write $S(D) = N(D)/C(D)$ in lowest terms, with $\deg N < \deg C = L(s)$. Then

$$Y(D) = \frac{P(D)\, N(D^d)}{C(D^d)}, \qquad \text{LC} = d \cdot L(s) - \deg \gcd\big(P(D),\, C(D^d)\big).$$

Only $P$ can share a factor with the denominator: $\gcd(N, C) = 1$ gives polynomials $u, v$ with $uN + vC = 1$, and substituting $D^d$ for $D$ in that identity shows $\gcd(N(D^d), C(D^d)) = 1$.

**Two phases.** For $d = 2$ with phase sequences $a$ and $b$, the identity is $Y(D) = A(D^2) + D \cdot B(D^2)$. If $a = b$, this is $(1+D) A(D^2)$, and $A(D^2) = (A(D))^2$. Let $a$ be purely periodic and write $A(D) = N(D)/C(D)$ in lowest terms, with $\deg N < \deg C = L(a)$, so that

$$Y(D) = \frac{(1+D)\, N(D)^2}{C(D)^2}.$$

Since $\gcd(N^2, C^2) = 1$, the common factor of numerator and denominator is $\gcd(1+D,\, C^2)$, which is $1+D$ if $(1+D) \mid C$ and $1$ otherwise. There are two cases.

- If $(1+D) \nmid C$, the fraction is already reduced. Its denominator has degree $2L(a)$, so $\text{LC} = 2L(a)$.
- If $C = (1+D)^e C'$ with $e \geq 1$ and $(1+D) \nmid C'$, one factor cancels: $Y = N^2 / \big((1+D)^{2e-1} C'^2\big)$, which is reduced because $(1+D) \nmid N$. So $\text{LC} = 2L(a) - 1$.

In words: the LC is $2L(a) - 1$ if $1$ is a root of $a$, and $2L(a)$ otherwise. Every root other than $1$ doubles its Jordan size, and the root $1$, of Jordan size $e$, goes to Jordan size $2e - 1$. The factor $(1+D)$ never raises the LC above $2L(a)$. Small cases: the all-ones sequence gives LC $1$; the sequence $1, 0$ repeating gives $3$; $1, 1, 0$ repeating gives $4$; an m-sequence of degree 3 gives $6$. Measured by exhaustive computation over every sequence of period at most 8 (502 cases), the LC was $2L$ in 96 of them and $2L - 1$ in 406. The rule is asserted in `test/Tools/test_clock_controlled.py`; the full sweep is not part of the tests.

### LC Bound Procedure

The rules above combine into a procedure for bounding the LC and period of an interleaving at rate $d = 2^a d'$.

1. **Phase roots.** Find the root set $R_r$ and the Jordan sizes of each phase $s_r$.
2. **Fan by the odd part.** Replace each root by its fan of $d'$ roots: $R_r^{\text{fan}} = \bigcup_{\rho \in R_r} F_{d'}(\rho)$.
3. **Scale by the even part.** Multiply all Jordan sizes by $2^a$.
4. **Union.** Form $R_Y = \bigcup_r R_r^{\text{fan}}$.
5. **Cancel (optional).** Write $Y = N/C$ over the common denominator $C$, with combined numerator $N = \sum_r D^r N_r(D^d)$, and reduce by $\gcd(N, C)$. Each cancelled factor lowers the Jordan size of its root by one, and a root disappears only when its size reaches zero. Cancellation is decided by the combined numerator, not by pairs of phases. Skipping this step leaves an upper bound.
6. **Count.** Partition $R_Y$ into Frobenius orbits; then $\text{LC} \leq \sum_{\text{orbits}} |\text{orbit}| \cdot (\text{max Jordan size})$.
7. **Period.** Take the $\text{lcm}$, over the surviving roots, of $q \cdot 2^{\lceil \log_2 m \rceil}$, where $q$ is the order of the root and $m$ its Jordan size.

### Worked Example

Interleave two distinct period-3 sequences $a$ and $b$ at rate $d = 2 = 2^1 \cdot 1$. The three candidates are the sequences repeating $1,1,0$, $0,1,1$ and $1,0,1$, with D-transforms $P(D)/C(D)$ for $C = 1 + D + D^2$ and $P = 1$, $D$, $1+D$ respectively. All have root set $\{\alpha^{1/3}, \alpha^{2/3}\}$ with Jordan size 1. The seven steps give:

1. Both phases have roots $\{\alpha^{1/3}, \alpha^{2/3}\}$, of Jordan size 1.
2. The odd part is $d' = 1$, so there is no fanning.
3. The even part is $2^1$, so the Jordan sizes double to 2.
4. $R_Y = \{\alpha^{1/3}, \alpha^{2/3}\}$, each with Jordan size 2.
5. The phases share roots, so cancellation is possible. With $A = P_a/C$ and $B = P_b/C$, the polyphase identity gives $Y(D) = N(D)/C(D)^2$ with $N(D) = P_a(D^2) + D \cdot P_b(D^2)$ of degree at most 3. The minimal polynomial divides $(1+D+D^2)^2$.
6. There is one orbit of size 2 with Jordan size at most 2, so LC $\leq 2 \cdot 2 = 4$. The root $\alpha^0$ does not occur.
7. The period divides $3 \cdot 2 = 6$.

Whether the bound is attained depends on the pair. For three of the six ordered pairs, $N$ is divisible by $C$ and one factor cancels: $(P_a, P_b) = (1, D)$ gives $N = 1 + D^3$, $(D, 1+D)$ gives $N = D(1+D+D^2)$, and $(1+D, 1)$ gives $N = 1+D+D^2$. These outputs have LC 2 and period 3; for instance $1,1,0$ interleaved with $0,1,1$ is $1,0,1,1,0,1,\ldots$, the sequence repeating $1,0,1$. The other three ordered pairs give LC 4 and period 6.

## Shrinking

### Definition

Given a **control** sequence $a[t]$ and a **data** sequence $b[t]$, the **shrinking** of $b$ by $a$ keeps only the data samples at which the control is 1:

$$y[j] = b[\tau(j)], \qquad \tau(j) = j\text{-th time } t \text{ with } a[t] = 1.$$

The output has variable rate: each period of the control contributes $w$ output samples, where $w$ is the number of ones in one period of $a$.

### Reduction to Interleaving

Shrinking becomes an interleaving once the output is grouped by control periods. Suppose $a$ has period $q_A$, with its $w$ ones at positions $t_0 < \cdots < t_{w-1}$ within a period. Then in the $k$-th control period

$$y[kw + j] = b[t_j + kq_A], \qquad 0 \leq j < w.$$

Define the **phase-$j$ slow-time sequence** $s_j[k] = b[t_j + kq_A]$. Comparing with the definition of interleaving, the shrunk output is the $w$-phase interleaving of $s_0, \ldots, s_{w-1}$.

### Phase Roots

Each phase $s_j$ samples $b$ at the fixed offset $t_j$ within a control period and advances by $q_A$ steps from one period to the next. Suppose every root of $b$ has Jordan size 1, so that $b[t] = \sum_\rho c_\rho \rho^t$ over its root set $R_B$. Then

$$s_j[k] = \sum_{\rho \in R_B} \underbrace{c_\rho \rho^{t_j}}_{\text{coefficient}} \cdot \underbrace{(\rho^{q_A})^k}_{\text{slow-time root}}$$

The slow-time roots are $\{\rho^{q_A} : \rho \in R_B\}$, the data roots [decimated](Decimation%20Expansion%20and%20Linear%20Complexity.md#root-rule) by $q_A$. Write $q_B$ for the period of $b$ and $L_B$ for its LC. When $\gcd(q_A, q_B) = 1$ this decimation is a bijection — distinct roots stay distinct — so each phase has $L_B$ roots.

Interleaving these $w$ phases then fans and scales each root according to the odd–even decomposition of $w$, as in the procedure above. This is the source of the shrinking-generator bound below.

## Circuit Types

Clock control enters a circuit in two ways: the transition applied at each step can vary, or the transition can stay fixed while only some outputs are recorded. This section treats the first and then contrasts the two.

### Clock-Controlled Registers (Periodic Schedule)

Consider a register with $d$ transition maps $U_0, \ldots, U_{d-1}$ on a state space $V$, applied in a repeating cycle. Augment the state with the phase counter, $\widetilde{V} = V \times \mathbb{Z}/d\mathbb{Z}$. After one full cycle the state has advanced by the **monodromy** $M = U_{d-1} \circ \cdots \circ U_0$.

With $h$ the output map and $v_0$ the initial state, the output at fast time $t = dk + r$ is $y[dk+r] = h(U_{r-1} \circ \cdots \circ U_0(M^k(v_0)))$. For each fixed $r$ this is a slow-time sequence in $k$ governed by $M$, and $y$ is the $d$-phase interleaving of those $d$ sequences.

The schedule is a **uniform clock** when $U_r = U$ for all $r$. Then $M = U^d$, and reconstructing the fast-time output from the $d$ phases is equivalent to expansion: root fans for odd $d$, Jordan scaling for even $d$. Under a **non-uniform clock**, $M$ is a product of distinct operators, and its eigenstructure must be computed directly.

### Controlled Sampling vs. Controlled Clocking

In **controlled sampling**, of which shrinking is the example, the register advances uniformly and the output is recorded only at selected times. Analysing it is a decimation and interleaving problem. In **controlled clocking** the transition itself changes from step to step, and the analysis requires the monodromy.

## LC Formulas for Classical Constructions

The three generators below are implemented, for arbitrary input sequences, in `Tools/ClockControlled.py`, and the statements checked computationally are in `test/Tools/test_clock_controlled.py`.

### Shrinking Generator (SG)

The shrinking generator shrinks the output of a data LFSR $B$ (period $q_B$, LC $L_B$) by the output of an independent control LFSR $A$ (period $q_A$, with $w$ ones per period).

By the reduction above, the output is a $w$-phase interleaving of decimated copies of the data, and its LC is bounded by

$$\boxed{L(\text{SG}) \leq w \cdot L_B}$$

The bound follows from §Phase Roots and the polyphase identity, under the hypotheses of that section: the data sequence is purely periodic and every root of it has Jordan size 1. Every phase $s_j$ has its roots in $\{\rho^{q_A} : \rho \in R_B\}$, a set of at most $L_B$ elements, so each D-transform can be written $S_j(D) = N_j(D)/C'(D)$ over one common denominator $C'$ with $\deg N_j < \deg C' \leq L_B$. Then

$$Y(D) = \sum_{j=0}^{w-1} \frac{D^j\, N_j(D^w)}{C'(D^w)}$$

has a denominator of degree $w \cdot \deg C' \leq w \cdot L_B$, and each numerator term has degree $\deg D^j N_j(D^w) \leq (w - 1) + w(\deg C' - 1) < w \cdot \deg C'$, below that of the denominator. This gives the bound. Coprimality of $q_A$ and $q_B$ is not needed for it.

The bound is not an equality in general, even when the periods are coprime. Take the control $1,1,1,0,0,0,0,0,0,0$ repeating, with $q_A = 10$ and $w = 3$, and for the data an m-sequence of degree 3, with $q_B = 7$ and $L_B = 3$. The selected times are $10k + j$ for $j = 0, 1, 2$, and $10k + j \equiv 3k + j \pmod 7$, so $y[3k + j] = b[10k + j] = b[3k + j]$: the output is the data sequence itself, of LC 3 against $w \cdot L_B = 9$. This case is asserted in `test/Tools/test_clock_controlled.py`.

When $\gcd(q_A, q_B) \neq 1$, the decimation $\rho \mapsto \rho^{q_A}$ can collapse distinct roots, and the LC can then fall further below $w \cdot L_B$. It need not collapse them: with $q_A = 3$ and a data m-sequence of degree 4, the roots $\rho, \rho^2, \rho^4, \rho^8$ cube to four distinct elements of order 5.

In the odd–even decomposition $w = 2^a w'$ with $w'$ odd, the factor $w'$ fans roots, producing new eigenvalues, and the factor $2^a$ scales Jordan sizes.

**Example.** Let the control be an m-sequence of degree $n_A$, so that $q_A = 2^{n_A}-1$ and $w = 2^{n_A-1}$, and let the data be an m-sequence of degree $n_B$ with $\gcd(n_A, n_B) = 1$, which ensures $\gcd(q_A, q_B) = 1$. Then

$$L(\text{SG}) \leq 2^{n_A-1} \cdot n_B.$$

The bound is attained in the four cases below, which `test_clock_controlled.py` asserts:

| $n_A$ | $n_B$ | $w$ | LC |
|---|---|---|---|
| 3 | 4 | 4 | 16 |
| 3 | 5 | 4 | 20 |
| 4 | 3 | 8 | 24 |
| 4 | 5 | 8 | 40 |

It is not always attained, even for m-sequences, and whether it is attained depends on the pair of primitive polynomials and not on the phases of the two sequences. Shifting control and data together by one step shifts the output by $a[0] \in \{0, 1\}$ positions, and a shift of a purely periodic sequence has the same LC. When $\gcd(q_A, q_B) = 1$ the joint shift $(1, 1)$ generates $\mathbb{Z}/q_A \times \mathbb{Z}/q_B$, so every pair of phases is reached by joint shifts, and every phase pair gives a shift of one and the same output sequence.

Measured by exhaustive computation over every pair of primitive polynomials and every relative phase, the bound was attained in all cases for $(n_A, n_B) = (2,3)$, $(3,2)$, $(3,4)$, $(4,3)$, $(2,5)$, $(3,5)$ and $(5,2)$. For the other three parameter pairs measured, some polynomial pairs fall short:

| $n_A$ | $n_B$ | bound | polynomial pairs below the bound | their LC |
|---|---|---|---|---|
| 5 | 3 | 48 | 2 of 12 | 45 |
| 4 | 5 | 40 | 2 of 12 | 35 |
| 5 | 4 | 64 | 6 of 12 | 60 |

The row $(4, 5)$ of the first table is therefore a statement about the polynomial pair the test uses, not about every such pair. The shortfall at $(5, 3)$ is asserted in `test/Tools/test_clock_controlled.py`; the full sweep is not part of the tests.

### Self-Shrinking Generator (SSG)

The self-shrinking generator shrinks a sequence $a$ by itself: consecutive bits are paired, and the first of each pair serves as the control for the second.

$$c[t] = a[2t], \qquad x[t] = a[2t+1], \qquad y[j] = x[\tau(j)] \text{ where } c[\tau(j)] = 1.$$

Take $a$ to be an m-sequence of degree $n$, with period $q = 2^n - 1$. Both $c$ and $x$ are [decimations](Decimation%20Expansion%20and%20Linear%20Complexity.md#root-rule) of $a$ by 2, with different phase offsets. Since $q$ is odd, $\gcd(2, q) = 1$, so both have period $q$ and LC $n$. The control has $w = 2^{n-1}$ ones per period, by the balance property of m-sequences.

The shrinking-generator bound gives $L \leq 2^{n-1} \cdot n$, but this overestimates. Control and data here have the same period, $q_A = q_B = q$, so $\gcd(q_A, q_B) = q$: every data root $\rho$ satisfies $\rho^{q} = 1$, and the decimation of §Phase Roots sends every root to $1$.

The known bound (Meier–Staffelbach 1994, `MeierStaffelbach1994`) is

$$L(\text{SSG}) \leq 2^{n-1},$$

a factor of $n$ tighter, and it follows from the reduction of §Reduction to Interleaving. The phase sequences are $s_j[k] = x[t_j + kq] = x[t_j]$, because $x$ has period $q$: every phase is constant. So $y[kw + j] = x[t_j]$ does not depend on $k$, and $y$ has period dividing $w = 2^{n-1}$. In terms of the shift, $(1 + D^{2^{n-1}})\,y = 0$, and in characteristic 2 the Frobenius identity gives $1 + D^{2^{n-1}} = (1+D)^{2^{n-1}}$. The minimal polynomial of $y$ divides a polynomial of degree $2^{n-1}$, hence $L \leq 2^{n-1}$.

The bound holds but was not attained in any case measured. Measured by exhaustive computation over every m-sequence of degree $n$ and every phase, the LC was:

| $n$ | measured LC | bound $2^{n-1}$ |
|---|---|---|
| 3 | 2 to 3 | 4 |
| 4 | 5 | 8 |
| 5 | 10 to 13 | 16 |
| 6 | 25 to 28 | 32 |
| 7 | 54 to 59 | 64 |

That the LC lies below the bound is asserted in `test/Tools/test_clock_controlled.py` for $n = 3, 4, 5$; the full sweep is not part of the tests.

### Alternating-Step Generator (ASG)

The alternating-step generator uses a control LFSR $A$ (period $q_A$, $w$ ones per period) to clock two data LFSRs, $B$ with LC $L_B$ and $C$ with LC $L_C$. When $A[t] = 1$ the register $B$ advances; when $A[t] = 0$ the register $C$ advances. The output is

$$y[t] = B[N_1(t)] \oplus C[N_0(t)]$$

where $N_1(t)$ counts the ones and $N_0(t)$ the zeros of $A$ before time $t$.

This is not a pair of shrinking generators. The SG keeps only the $w$ accepted samples, whereas the ASG produces an output at every step. At positions where $A = 0$, register $B$ is not absent — it **holds its last value**. Because of this sample-and-hold effect, the contribution of $B$ spans all $q_A$ output positions of a control period, not just the $w$ at which $B$ advances. Each held sequence is nonetheless an interleaving, with $q_A$ phases rather than $w$, and that gives the bound

$$\boxed{L(\text{ASG}) \leq q_A \cdot (L_B + L_C)}$$

The derivation assumes $B$ and $C$ purely periodic with every root of Jordan size 1. Write $u[t] = B[N_1(t)]$ for the held contribution of $B$. Over one control period $N_1$ increases by $w$, so $N_1(q_A k + r) = wk + N_1(r)$ for $0 \leq r < q_A$, and with $B[t] = \sum_\rho c_\rho \rho^t$,

$$u[q_A k + r] = B[wk + N_1(r)] = \sum_{\rho \in R_B} c_\rho\, \rho^{N_1(r)}\, (\rho^{w})^k.$$

So $u$ is the $q_A$-phase interleaving of decimations of $B$ by $w$, taken at the offsets $N_1(r)$. Every phase has its roots in $\{\rho^w : \rho \in R_B\}$, at most $L_B$ of them, and the polyphase identity, used as in the shrinking-generator bound, gives $L(u) \leq q_A L_B$. The same argument applied to $C[N_0(t)]$, which advances $q_A - w$ times per control period, gives a bound of $q_A L_C$. The LC of an XOR of two sequences is at most the sum of their LCs, so $L(y) \leq q_A (L_B + L_C)$.

For m-sequence control, $q_A = 2^{n_A} - 1$ and the bound reads $L(\text{ASG}) \leq (2^{n_A}-1)(L_B + L_C)$. Günther ("Alternating Step Generators Controlled by De Bruijn Sequences", EUROCRYPT '87, LNCS 304, 1988) treats the case of a de Bruijn control sequence of period $2^k$ and gives the upper bound $L \leq (m_1 + m_2)\,2^k$ for data registers of lengths $m_1$ and $m_2$; the derivation above is the corresponding upper bound for a general control period $q_A$.

At $n_A = 3$, $n_B = 4$, $n_C = 5$ the LC, measured by exhaustive computation over all 24 triples of primitive polynomials, is $63 = 7 \times 9$, which meets the bound exactly. The value 63 at these parameters is asserted in `test/Tools/test_clock_controlled.py`; the full sweep is not part of the tests.

## Connections

- **[Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md):** The root-level mechanics behind interleaving (expansion = root fans + Jordan scaling) and phase extraction (decimation = root mapping).
- **[Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md):** Root labels, Frobenius orbits, odd-denominator constraint.
- **[Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md):** Jordan blocks and multiplicity brackets for even interleaving rates.
- **[Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md):** The XOR combination rule used in phase cancellation and ASG decomposition.

## Source

| | |
|---|---|
| `Tools/ClockControlled.py` | `shrinking_generator`, `self_shrinking_generator`, `alternating_step_generator` |
| `test/Tools/test_clock_controlled.py` | LC checks for the three generators |
