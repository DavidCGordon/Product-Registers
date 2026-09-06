# Polynomial Types

Three types of polynomial arise naturally when working with field elements and feedback registers: minimal, characteristic, and primitive. They are related but distinct, and each plays a specific role in register analysis. This document defines them and traces their relationships.

For the algebraic closure and field structure these polynomials live over, see [Roots and the Algebraic Closure](Roots%20and%20the%20Algebraic%20Closure.md). For how polynomials relate to sequences (primal vs dual), see [Polynomial Conventions](../conventions/Polynomial%20Conventions.md). For notation conventions, see [Notation and Terminology](../conventions/Notation%20and%20Terminology.md).

## Minimal Polynomial

The **minimal polynomial** $m_\alpha(x)$ of an element $\alpha \in \overline{\mathbb{F}_2}$ is the unique monic polynomial of smallest degree over $\mathbb{F}_2$ that has $\alpha$ as a root. Its key properties:

- It is irreducible over $\mathbb{F}_2$.
- Its degree equals $\text{Ord}(q)$, the size of the Frobenius orbit of $\alpha$ (see [Roots — Frobenius](Roots%20and%20the%20Algebraic%20Closure.md#the-frobenius-automorphism-and-conjugates)).
- Its roots are exactly the conjugates of $\alpha$: $\{\alpha, \alpha^2, \alpha^{2^2}, \ldots, \alpha^{2^{d-1}}\}$.
- It divides every polynomial over $\mathbb{F}_2$ that has $\alpha$ as a root.

For the generator $\alpha$ of the representation $\mathbb{F}_2[x]/\langle P \rangle$, the minimal polynomial is $P$ itself — this is the tautological case where representation polynomial and minimal polynomial coincide.

## Characteristic Polynomial

The **characteristic polynomial** of an element $\alpha$, viewed as the linear map $M_\alpha$ on $\mathbb{F}_2^n$, is $\det(xI - M_\alpha)$. Its degree is $n$ (the ambient vector space dimension), regardless of the algebraic degree of $\alpha$.

When $\alpha$ generates $\mathbb{F}_{2^n}$ (i.e., its minimal polynomial has degree $n$), the characteristic and minimal polynomials coincide. When $\alpha$ lives in a proper subfield $\mathbb{F}_{2^d} \subset \mathbb{F}_{2^n}$ with $d < n$, the characteristic polynomial is $m_\alpha(x)^{n/d}$ — the minimal polynomial raised to the power $n/d$, accounting for the repeated eigenvalue structure across the copies of $\mathbb{F}_{2^d}$ inside $\mathbb{F}_{2^n}$.

In register terms, the characteristic polynomial of the update matrix governs the recurrence that the full register state satisfies. For a simple LFSR with a primitive polynomial, the characteristic and minimal polynomials are the same. For a CMPR with multiple blocks, the characteristic polynomial of the full system is the product of the blocks' individual polynomials — possibly with shared factors when blocks have the same size.

## Primitive Polynomial

A **primitive polynomial** is an irreducible polynomial $P$ of degree $n$ whose roots are primitive elements of $\mathbb{F}_{2^n}$ — generators of the multiplicative group $\mathbb{F}_{2^n}^\times$. Equivalently, $P$ is primitive iff the element $x$ in $\mathbb{F}_2[x]/\langle P \rangle$ has multiplicative order exactly $2^n - 1$.

Key facts about primitive polynomials:

- Every finite field $\mathbb{F}_{2^n}$ has primitive polynomials of degree $n$; there are $\phi(2^n - 1)/n$ of them, where $\phi$ is Euler's totient.
- The set of primitive polynomials is closed under reversal (reciprocal): if $P(x)$ is primitive, so is $x^n P(1/x)$ (see [Polynomial Conventions — Symmetry Relations](../conventions/Polynomial%20Conventions.md#symmetry-relations)).
- LFSRs with primitive polynomials achieve the maximum period $2^n - 1$.
- MPR blocks in a CMPR are always constructed from primitive polynomials, ensuring maximal period for each component.

The library's `NumberTheory` module provides `num_primitive_polynomials(n)` to count them and `generate_primitive_polynomials(base)` to enumerate them by decimation from a known primitive polynomial (see [Roots — Arithmetic](Roots%20and%20the%20Algebraic%20Closure.md#arithmetic-in-exponential-form) for how decimation works).

In the [Notation and Terminology](../conventions/Notation%20and%20Terminology.md) conventions: "primitive polynomial" is the default noun for the polynomial associated with a register. Use "primal" / "dual" only when disambiguating which exact relationship that association entails. For many properties there is no difference, but it can be important when trying to replicate exact circuits or hardware.
