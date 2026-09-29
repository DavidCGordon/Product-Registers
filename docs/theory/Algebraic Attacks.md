# Algebraic Attacks — Mathematical Foundations

This document covers the shared mathematical framework underlying the three algebraic attacks implemented in PyPR: the Naive Algebraic Attack (NAA), the Reduced Algebraic Attack (RAA), and the Fast Algebraic Attack (FAA). It also covers the cube attack, which derives its equations differently but solves them with the same stores and solvers.

For implementation details (offline/online phases, output dict keys, code paths), see [architecture/Components Architecture](../architecture/Components_Architecture.md) and [architecture/Attack Compatibilities](../architecture/Attack_Compatibilities.md).

## Core Idea

All three attacks recover the initial state $s_0$ of a feedback register by generating a system of polynomial equations over GF(2) that relate the unknown state variables to the observed keystream, then solving that system.

At clock cycle $t$, the register is in state $s_t$, which is a deterministic polynomial function of $s_0$ (determined entirely by the feedback function). The output at time $t$ is:

$$k_t = f(s_t)$$

Expanding $f(s_t)$ as a polynomial in $s_0$, each clock cycle produces one equation in the initial-state variables.

## Attack Hierarchy

### NAA — Direct Approach

Uses $f$ directly. The system has degree $\deg(f)$, so the number of monomials (unknowns) is at most $\binom{n}{\leq \deg(f)}$.

### RAA — Degree Reduction via Annihilators

Works with a **low-degree pair** $(g, h)$ satisfying $f \cdot g = h$, where both $g$ and $h$ have lower degree than $f$. The identity $h(s_t) = k_t \cdot g(s_t)$ allows combining the equations for $g$ and $h$ at each clock cycle using the observed keystream.

When $h = 0$ ($g$ is a true annihilator of $f$), the attack degrades gracefully: equations are produced only when $k_t = 1$.

See [Monomial Profile Theory](Monomial%20Profile%20Theory.md) and [Root Expressions and LC Estimation](Root%20Expressions%20and%20LC%20Estimation.md) for how the degree reduction affects the monomial space and linear complexity bounds.

### FAA — Linear Recurrence Exploitation

Extends RAA by exploiting the linear recurrence of $h = f \cdot g$. Because $h$ satisfies a recurrence of length $L = \text{LC}(h)$, the attacker can form combined equations from $L$-term linear combinations of adjacent keystream-weighted $g$ equations. This typically reduces the required keystream length compared to RAA.

The FAA has a subtle rank-loss phenomenon: the set of "equal-or-better" low-degree pairs forms a subspace $W$ that removes exactly $\dim(W)$ columns from the equation matrix, regardless of keystream. See [Trajectory and Column LC](Trajectory%20and%20Column%20LC.md) for how this interacts with CMPR structure.

### Cube Attack — Superpolys over Tweakable Bits

The cube attack applies when some bits of $s_0$ are **tweakable**: set freely by the attacker, like an IV. Split the state bits into the known bits $K$ (every tweakable bit, at its public value) and the unknown bits $U$, and write $p_t$ for the output at keystream position $t$ as a polynomial in $s_0$.

For a **cube** $I$ of tweakable bits, sum $p_t$ over the $2^{|I|}$ assignments of the bits in $I$, holding every other bit fixed. A monomial $x^m$ of $p_t$ contributes $x^{m \setminus I} \cdot \bigoplus_{v} \prod_{i \in m \cap I} v_i$, and that inner sum counts the assignments with $v_i = 1$ on $m \cap I$: there are $2^{|I| - |m \cap I|}$, odd exactly when $m \supseteq I$. So the sum is

$$S_{I,t} = \bigoplus_{m \supseteq I} x^{m \setminus I} = P_{I,t}(x),$$

the **superpoly** of $I$ in the [cube identity](Cube%20Equation%20Generation.md#background-the-cube-attack-identity) $p_t = T_I \cdot P_{I,t} + Q_t$, evaluated at the non-cube bits. With the known bits at their values, $P_{I,t}$ is a polynomial in $x_U$ alone. The offline phase computes it on a simulated register; the online phase measures $S_{I,t}$ on the target; each $(I, t)$ whose $P_{I,t}$ is not constant gives the equation $P_{I,t}(x_U) = S_{I,t}$. The known bits enter the system as the equations $x_i = v_i$, so the solver recovers the whole state.

**Degree bound from the monomial profile.** A profile term with per-block counts $c$ covers monomials drawing at most $c_b$ variables from block $b$ ([Monomial Profile Theory](Monomial%20Profile%20Theory.md#what-is-a-monomial-profile)). Let a cube draw $k_b$ variables from block $b$. A monomial of the term can contain $T_I$ only if $c_b \geq k_b$ in every block, and then $x^{m \setminus I}$ draws at most $c_b - k_b$ variables from block $b$. So the term contributes superpoly monomials of degree at most its **excess** $e = \sum_b (c_b - k_b)$, over the blocks where $c_b > k_b$, and the superpoly's degree is at most the largest excess among the terms containing the cube. `MonomialProfile.get_cube_candidates` forms candidates by removing one variable from one block of a profile term and returns this bound, and those blocks, with each. The comparison ranges over every block either side mentions: a containing term with variables in a block the candidate does not use contributes those variables too. Comparing only the candidate's own blocks classified such cubes as linear and produced wrong equations.

**Exact recovery up to the bound.** Evaluating $P$ at the indicator $e_s$ of a set $s$ of unknown bits (all other unknown bits $0$) sums the coefficients of the monomials inside $s$: $P(e_s) = \bigoplus_{u \subseteq s} \mathrm{coef}(u)$. [Möbius inversion](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md) inverts this exactly for every $s$: $\mathrm{coef}(s) = \bigoplus_{u \subseteq s} P(e_u)$, each $P(e_u)$ being one cube sum. The degree bound $d$ is what makes the points $|s| \leq d$ over the target blocks' unknown bits sufficient: every coefficient outside them is zero, so the recovered coefficients are all of $P$.

**One-to-next chaining.** On CMPRs chained by `arman_template` (block sizes 7/5/3, 7/5/3/2 and 13/7/5/3, four seeds each) no candidate's cube lies outside every block its superpoly involves, at any degree, so a cube over IV blocks yields no equation on a separate key block. The mechanism is the trading of degree into the next block described in [Mesh Optimization](../architecture/Mesh%20Optimization.md). Linear candidates do occur, but their cube shares a block with its own superpoly. The author confirms this holds for one-to-next chaining generally, as a consequence of that trading logic. **OPEN:** the derivation has not been written out here.

## Equation Generation

Three generators produce equations from register clocking:

- **CubeEqGenerator** — Numba-JIT cube-sum algorithm; 2–3 orders of magnitude faster for CMPRs when monomial profiles are known. See [Cube Equation Generation](Cube%20Equation%20Generation.md).
- **SubstitutionEqGenerator** — Composes per-bit equations; discovers monomials dynamically.
- **SymbolicEqGenerator** — Composes in the opposite direction; supports partial initialization.

See [architecture/Components Architecture](../architecture/Components_Architecture.md) §4 for the generator taxonomy and store compatibility.

## Equation Solving

Once enough independent equations are collected, the system is solved via:

1. **Reduction** — LU back-substitution, Gaussian elimination (RREF), or Groebner basis computation.
2. **Guess-and-prune** — Exhaustive search over unpivoted variables, pruned by checking partial solutions against additional keystream.

See [architecture/Components Architecture](../architecture/Components_Architecture.md) §5 for the solver taxonomy and [architecture/Attack Compatibilities](../architecture/Attack_Compatibilities.md) for per-attack solver constraints.
