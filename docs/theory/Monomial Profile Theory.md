# Monomial Profile Theory

## What is a Monomial Profile?

The ANF of a signal coming out of a CMPR can contain an enormous number of monomials. Rather than tracking the exact set, we track an **upper bound** on which monomials can appear. This is analogous to the root expression approach for linear complexity: instead of computing which roots are actually present, we bound the set of possible roots.

A **monomial profile** is a set of pairs $\langle e \cdot w \rangle$, where each pair means: "this signal can contain monomials that draw 1 to $w$ variables from the $e$-bit component register."

A full monomial profile for a signal is a collection of such products (one product per "term" in the expression), where each product combines contributions from different registers.

**Example:** Let $C$ be a 12-bit CMPR with a 7-bit MPR (bits 5--11) and a 5-bit MPR (bits 0--4), using the [high-to-low block ordering](../conventions/Notation%20and%20Terminology.md) standard in the library. The profile $\langle 7 \cdot 4 \rangle \langle 5 \cdot 2 \rangle$ covers all monomials with 1–4 variables from the 7-bit register and 1–2 variables from the 5-bit register — for example, $c_0 c_5$, $c_0 c_3 c_8 c_9 c_{11}$, $c_0 c_5 c_6 c_7 c_8$.

## The Five Core Propositions

These give the rules for counting and combining monomial profiles. They parallel the rules for [root expressions](Root%20Expressions%20and%20LC%20Estimation.md) almost exactly.

---

**Proposition 1** (Counting): $|\langle e \cdot w \rangle| = \sum_{i=1}^{w} \binom{e}{i}$

*Proof:* For each degree $i$ from 1 to $w$, there are $\binom{e}{i}$ ways to choose $i$ variables from $e$ bits. Summing over all possible degrees gives the total count.

---

**Proposition 2** (Same-register product): If all pairs share the same register size $e$, then

$$\prod_{i=1}^{k} \langle e \cdot w_i \rangle = \left\langle e \cdot \sum_{i=1}^{k} w_i \right\rangle$$

*Proof:* All variables come from the same register, so the product of monomials of lengths $\leq w_1, \ldots, w_k$ is just a monomial of length at most $\sum w_i$ (with equality when all chosen variables are distinct).

---

**Proposition 3** (Cross-register product): If all register sizes $e_i$ are distinct, then

$$\left| \prod_{i=1}^{k} \langle e_i \cdot w_i \rangle \right| = \prod_{i=1}^{k} |\langle e_i \cdot w_i \rangle|$$

*Proof:* Variables from different registers are always distinct from each other, so there are no shared variables between terms. The counts multiply independently.

---

**Proposition 4** (Sum / inclusion-exclusion): For any two monomial expressions $E_1$ and $E_2$,

$$|E_1 + E_2| = |E_1| + |E_2| - |E_1 \cap E_2|$$

*Proof:* Directly from the principle of inclusion-exclusion on the sets of monomials.

---

**Proposition 5** (Intersection of single products): Let $E_1$ and $E_2$ each be a single product of pairs.

- If they involve different sets of registers, $E_1 \cap E_2 = \emptyset$ — any monomial in $E_1$ must contain at least one variable from a register not present in $E_2$, so no overlap is possible.
- If they involve the same set of registers $\{e_1, \ldots, e_k\}$, then:

$$E_1 \cap E_2 = \prod_{i=1}^{k} \langle e_i \cdot \min(w_{1,i}, w_{2,i}) \rangle$$

*Proof:* Shared monomials can include at most $\min(w_{1,i}, w_{2,i})$ variables from register $i$, since they must satisfy both bounds simultaneously.

---

## Similarity to Root Expressions

The algebra of monomial profiles and [root expressions](Root%20Expressions%20and%20LC%20Estimation.md) is parallel at the counting level. Both have five combination rules and represent upper bounds on a set (monomials vs. roots), but root expressions additionally carry full-coset embedding and Jordan-length data. The [root multiplicity note](Root%20Multiplicities%20and%20Jordan%20Decomposition.md) explains why the parallel is a qualified analogy rather than an identity.

There are two concrete differences, which is why the two have separate implementations (`RootExpression` and `MonomialProfile` classes):

### 1. Full-Coset Degeneracy

For a root expression with register size $e$, the highest-weight pair $\langle e \cdot e \rangle$ includes the exponent $2^e - 1$. In the field $\mathbb{F}_{2^e}$, $\alpha^{2^e - 1} = 1$ — this is the identity, not a useful root. This creates random cancellations and also means there are only $2^e - 2$ usable roots outside the base field.

For monomial profiles, the full-weight monomial (all $e$ variables) is a perfectly valid term with no degeneracy. There are genuinely $2^e - 1$ distinct nonzero monomials in such a set. Although $2^e - 2$ vs $2^e - 1$ seems like a small difference, it compounds significantly when taking products across many registers.

### 2. Multiple Registers of the Same Size

When a CMPR contains multiple registers with the same size (and therefore the same minimal polynomial), root expressions require a Jordan-length mechanism to handle the algebraic interaction between the identical fields. This is implemented in the library (necessary for T-Functions) but is complicated, fragile, and excluded from the paper because it has more degenerate edge cases. See [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md) for the decomposition behind it.

For monomial profiles, the proofs only require that variables come from *different registers*, not differently-sized registers. In code, each pair tracks the index of the MPR it came from (not just its size). No multiplicity mechanism is needed.

In the simple case where all component sizes are distinct, neither multiplicities nor block indices matter, and the two systems are completely equivalent.

**Overall:** monomial profile upper bounds are higher than (but usually close to) root expression bounds, and the system is simpler to reason about and implement.

## Statistical Model for Expected Behavior

The upper bound is not always tight — some monomials in the profile may cancel out of the actual ANF. The question is: how much do they degenerate?

For **root expressions**, cancellation is unlikely because roots appear in Frobenius cosets. For a coset of size $n$, all $n$ roots must cancel simultaneously for the bound to be loose. Assuming statistical independence, this happens with probability only $1/2^n$. This makes the root expression bound very tight almost always.

For **monomial profiles**, monomials don't come in cosets. There is no structure that makes simultaneous cancellation rare. The natural model is that each monomial independently cancels with probability $1/2$.

Under this model, for a profile predicting $M$ possible monomials, the expected number of surviving monomials is $M/2$, with standard deviation $\sqrt{M}/2$.

*Example:* For the C17 CMPR used in the paper, the profile predicts $M = 12515$ monomials. The expected surviving count is $12515/2 = 6257.5$ with standard deviation $\sqrt{12515}/2 \approx 55.93$. Empirical ANF iteration confirmed this model matches the observed behavior well.

This statistical model also appears in the analysis of [algebraic attacks](Algebraic%20Attacks.md) and [cube-based equation generation](Cube%20Equation%20Generation.md), where the expected number of surviving monomials determines the difficulty of equation generation.

## API

- `C.monomial_profiles()` — returns a list of `MonomialProfile` objects, one per bit.
- `mp.upper()` — the upper bound on the number of monomials.
- `mp.lower()` — a lower estimate (holds with high probability).
- `BooleanFunction.eval_ANF(profiles)` — evaluate a function symbolically with monomial profiles as inputs, giving the output's profile. Works for any function using only AND, XOR, and CONST(1). The profiles also drive [cube-based equation generation cost estimates](Cube%20Equation%20Generation.md), where each profile's degree vector determines the per-monomial cube-sum cost.
