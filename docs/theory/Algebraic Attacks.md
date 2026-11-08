# Algebraic Attacks — Mathematical Foundations

This document covers the shared mathematical framework underlying the three algebraic attacks implemented in PyPR: the Naive Algebraic Attack (NAA), the Reduced Algebraic Attack (RAA), and the Fast Algebraic Attack (FAA). It also covers the cube attack, which derives its equations differently but solves them with the same stores and solvers.

The document assumes the [Algebraic Normal Form](Algebraic%20Normal%20Form.md) of a Boolean function. For implementation details (offline/online phases, output dict keys, code paths), see [Components Architecture](../architecture/Components_Architecture.md) and [Attack Compatibilities](../architecture/Attack_Compatibilities.md).

## Core Idea

All three algebraic attacks recover the initial state $s_0$ of an $n$-bit feedback register by generating a system of polynomial equations over $\mathbb{F}_2$ that relate the unknown state variables to the observed keystream, then solving that system.

At clock cycle $t$, the register is in state $s_t$, which is a deterministic polynomial function of $s_0$ (determined entirely by the feedback function). Writing $f$ for the output function of the register, the output at time $t$ is:

$$k_t = f(s_t)$$

Expanding $f(s_t)$ as a polynomial in $s_0$, each clock cycle produces one equation in the initial-state variables.

## Attack Hierarchy

### NAA — Direct Approach

The Naive Algebraic Attack uses $f$ directly. How many monomials (unknowns) the system has depends on the state update.

When the state update is linear, as in an LFSR or a filter generator built on a single MPR, each coordinate of $s_t$ is a linear function of $s_0$. Then $f(s_t)$ has degree at most $\deg(f)$ in $s_0$, and the number of monomials is at most $\sum_{k \leq \deg(f)} \binom{n}{k}$, the number of monomials of degree at most $\deg(f)$ in $n$ variables.

For a nonlinear update this bound fails. Take $f = x_0$ and an update with $(s_{t+1})_0 = (s_t)_1 (s_t)_2$: then $\deg(f) = 1$, but $k_1 = (s_0)_1 (s_0)_2$ has degree 2. For CMPRs the monomial set is bounded instead by the monomial profile ([Monomial Profile Theory](Monomial%20Profile%20Theory.md)).

### RAA — Degree Reduction via Annihilators

The Reduced Algebraic Attack works with a **low-degree pair** $(g, h)$ satisfying $f \cdot g = h$, where both $g$ and $h$ have lower degree than $f$. Evaluating $f \cdot g = h$ at $s_t$ and substituting $f(s_t) = k_t$ gives the identity $h(s_t) = k_t \cdot g(s_t)$, which allows combining the equations for $g$ and $h$ at each clock cycle using the observed keystream.

When $h = 0$ ($g$ is a true annihilator of $f$), the identity reads $0 = k_t \cdot g(s_t)$. A clock cycle with $k_t = 1$ gives the equation $g(s_t) = 0$; a clock cycle with $k_t = 0$ gives $0 = 0$ and no equation. The attack still runs, but equations are produced only when $k_t = 1$.

See [Monomial Profile Theory](Monomial%20Profile%20Theory.md) and [Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md) for how the degree reduction affects the monomial space and linear complexity bounds.

### FAA — Linear Recurrence Exploitation

The Fast Algebraic Attack extends RAA by exploiting the linear recurrence of $h = f \cdot g$. Because $h$ satisfies a recurrence of length $L = \text{LC}(h)$, the linear complexity of $h$, the attacker can form combined equations from linear combinations of adjacent keystream-weighted $g$ equations. A recurrence of length $L$ relates $L + 1$ consecutive terms, so each combined equation uses a window of $L + 1$ consecutive keystream bits. This typically reduces the required keystream length compared to RAA.

> **OPEN:** Under what conditions does FAA require less keystream than RAA, and by how much? — needs author input.
> The claim is not sourced in this documentation, and each combined equation consumes $L + 1$ consecutive keystream bits. Courtois's CRYPTO 2003 paper on fast algebraic attacks is the place to look for a citation; it has not been checked against this statement.

The FAA equation matrix is subject to a rank loss: the set of "equal-or-better" low-degree pairs forms a subspace $W$ that removes exactly $\dim(W)$ columns from the equation matrix, regardless of keystream. See [Trajectory and Column LC](Trajectory%20and%20Column%20LC.md) for how this interacts with CMPR structure.

> **OPEN:** The rank-loss statement is not defined or proved anywhere in this documentation — needs author input. Four questions are open:
> 1. "Equal-or-better" than what, and in which order on pairs $(g, h)$?
> 2. $W$ is a subspace of which vector space?
> 3. Which matrix is "the equation matrix", and does "removes columns" mean zero columns, linearly dependent columns, or a rank deficit of $\dim(W)$?
> 4. Where is "exactly $\dim(W)$" proved?

### Cube Attack — Superpolys over Tweakable Bits

The cube attack applies when some bits of $s_0$ are **tweakable**: set freely by the attacker, as an initialization vector (IV) is. Split the state bits into the known bits $K$ (every tweakable bit, at its public value) and the unknown bits $U$, and write $p_t$ for the output at keystream position $t$ as a polynomial in $s_0$.

For a **cube** $I$ of tweakable bits, sum $p_t$ over the $2^{|I|}$ assignments of the bits in $I$, holding every other bit fixed. A monomial $x^m$ of $p_t$ contributes $x^{m \setminus I} \cdot \bigoplus_{v} \prod_{i \in m \cap I} v_i$, and that inner sum counts the assignments with $v_i = 1$ on $m \cap I$: there are $2^{|I| - |m \cap I|}$, odd exactly when $m \supseteq I$. So the sum is

$$S_{I,t} = \bigoplus_{m \supseteq I} x^{m \setminus I} = P_{I,t}(x),$$

the **superpoly** of $I$ in the [cube identity](Cube%20Equation%20Generation.md#background-the-cube-attack-identity) $p_t = T_I \cdot P_{I,t} + Q_t$, evaluated at the non-cube bits. With the known bits at their values, $P_{I,t}$ is a polynomial in $x_U$ alone. The offline phase computes it on a simulated register; the online phase measures $S_{I,t}$ on the target; each $(I, t)$ whose $P_{I,t}$ is not constant gives the equation $P_{I,t}(x_U) = S_{I,t}$. The known bits enter the system as the equations $x_i = v_i$, so the solver recovers the whole state.

**Degree bound from the monomial profile.** A profile term with per-block counts $c$ covers monomials drawing at most $c_b$ variables from block $b$ ([Monomial Profile Theory](Monomial%20Profile%20Theory.md#what-is-a-monomial-profile)). Let a cube draw $k_b$ variables from block $b$. A monomial of the term can contain $T_I$ only if $c_b \geq k_b$ in every block, and then $x^{m \setminus I}$ draws at most $c_b - k_b$ variables from block $b$. So the term contributes superpoly monomials of degree at most its **excess** $e = \sum_b (c_b - k_b)$, over the blocks where $c_b > k_b$, and the superpoly's degree is at most the largest excess among the terms containing the cube. `MonomialProfile.get_cube_candidates` forms candidates by removing one variable from one block of a profile term and returns this bound, and those blocks, with each. The comparison ranges over every block either side mentions: a containing term with variables in a block the candidate does not use contributes those variables too. A comparison restricted to the candidate's own blocks would drop $c_b$ for those other blocks and so understate the degree bound. The recovery described next, run over the points $|s| \leq d$ with too small a $d$, would then omit nonzero coefficients and return a wrong superpoly.

**Exact recovery up to the bound.** Evaluating $P$ at the indicator $e_s$ of a set $s$ of unknown bits (all other unknown bits $0$) sums the coefficients of the monomials inside $s$: $P(e_s) = \bigoplus_{u \subseteq s} \mathrm{coef}(u)$. [Möbius inversion](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md) inverts this exactly for every $s$: $\mathrm{coef}(s) = \bigoplus_{u \subseteq s} P(e_u)$, each $P(e_u)$ being one cube sum. The degree bound $d$ is what makes the points $|s| \leq d$ over the target blocks' unknown bits sufficient: every coefficient outside them is zero, so the recovered coefficients are all of $P$.

**One-to-next chaining.** In one-to-next chaining, each block's chaining function takes its inputs from the single block immediately upstream of it; this is the arrangement `arman_template` produces, except where the upstream block is a single bit, a case the template handles separately. On CMPRs chained by `arman_template` (block sizes 7/5/3, 7/5/3/2 and 13/7/5/3, four seeds each) no candidate's cube lies outside every block its superpoly involves, at any degree. A cube taken over blocks of tweakable (IV) bits in $K$ therefore yields no equation on a separate block of unknown (key) bits in $U$. The mechanism is the trading of degree into the next block described in [Mesh Optimization](../architecture/Mesh%20Optimization.md). Linear candidates do occur, but their cube shares a block with its own superpoly. The library's author states that the same property holds for one-to-next chaining in general, as a consequence of that trading of degree; that claim is recorded here without proof.

> **OPEN:** The derivation of this property for one-to-next chaining in general, from the degree-trading mechanism, has not been written out.
> What is established here is the computation on the `arman_template` configurations listed above.

## Equation Generation

Three generators produce equations from register clocking:

- **`CubeEqGenerator`** computes coefficients by a Numba-JIT cube-sum algorithm. It requires the monomial profiles in advance. On CMPRs at steady state it was measured against symbolic composition at about 20× the speed in one run and about 6× in another for a 15-bit register (the two runs disagree; see the open question at the first table there), and at about 150× for a 17-bit register; the measurements are in [Cube Equation Generation](Cube%20Equation%20Generation.md#practical-comparison).
- **`SubstitutionEqGenerator`** composes the current per-bit equations with the feedback functions, keeping only the bits the output function uses, and discovers monomials dynamically.
- **`SymbolicEqGenerator`** composes in the opposite direction, the feedback functions with the current equations, and supports partial initialization: symbolic starting expressions, for instance when part of the key is known.

See [Components Architecture](../architecture/Components_Architecture.md) §4 for the generator taxonomy and store compatibility.

## Equation Solving

Once enough independent equations are collected, the system is solved via:

1. **Reduction** — LU back-substitution, Gaussian elimination (RREF), or Groebner basis computation.
2. **Guess-and-prune** — Exhaustive search over unpivoted variables, pruned by checking partial solutions against additional keystream.

See [Components Architecture](../architecture/Components_Architecture.md) §5 for the solver taxonomy and [Attack Compatibilities](../architecture/Attack_Compatibilities.md) for per-attack solver constraints.
