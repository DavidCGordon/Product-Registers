# Notation and Terminology

Quick-reference for library-wide conventions. For matrix and vector orientation, see [Matrix Indexing](Matrix%20Indexing.md). For the mathematical framework behind polynomial conventions (primal/dual, generating functions, hardware), see [Polynomial Conventions](Polynomial%20Conventions.md). For prose standards in theory docs and docstrings, see [Writing Standards](Writing%20Standards.md).

## Coefficient Lists

Polynomial coefficient lists are always `[a_0, a_1, ..., a_n]` representing $a_0 + a_1 x + \cdots + a_n x^n$. There is only one storage convention in the codebase. The "reversal" issue discussed in [Polynomial Conventions](Polynomial%20Conventions.md) concerns the primal/dual *interpretation*, not a different storage format.

## CMPR Block Ordering and T-Functions

CMPR chaining propagates **from upstream blocks into downstream ones** (high index to low index). Block 0 is the "source" with no chaining inputs; successive blocks layer chaining on top. `CMPR.blocks` returns bit indices ordered high-to-low.

PyPR's `TFunction` uses bit $n-1$ as the LSB, opposite the usual little-endian counter convention. This aligns the T-function's triangular dependency order with CMPR's high-to-low chaining flow. This highlights the connection between T-functions and CMPRs, and unifies different analyses


## Vocabulary

- **Primitive polynomial** is the default noun for the polynomial associated with a register. Use "primal" / "dual" only when disambiguating which relationship a polynomial has to a generated sequence.
- **Shift direction**: use "shift up" / "shift down" in user-facing descriptions; use "values progress toward bits with lower (higher) indices" in technical sections.
