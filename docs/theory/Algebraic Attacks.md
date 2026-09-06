# Algebraic Attacks — Mathematical Foundations

This document covers the shared mathematical framework underlying the three algebraic attacks implemented in PyPR: the Naive Algebraic Attack (NAA), the Reduced Algebraic Attack (RAA), and the Fast Algebraic Attack (FAA).

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
