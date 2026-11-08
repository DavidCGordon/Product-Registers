# Notation and Terminology

This document is a quick reference for the conventions used throughout the library and its documentation: how polynomial coefficients are stored, how CMPR blocks and T-function bits are ordered, and which terms are reserved for which objects. For matrix and vector orientation, see [Matrix Indexing](Matrix%20Indexing.md). For the mathematical framework behind polynomial conventions (primal/dual, generating functions, hardware), see [Polynomial Conventions](Polynomial%20Conventions.md). For prose standards in theory docs and docstrings, see [Writing Standards](Writing%20Standards.md).

## Coefficient Lists

Polynomial coefficient lists are always `[a_0, a_1, ..., a_n]` representing $a_0 + a_1 x + \cdots + a_n x^n$. There is only one storage convention in the codebase. The "reversal" issue discussed in [Polynomial Conventions](Polynomial%20Conventions.md) concerns the primal/dual *interpretation*, not a different storage format.

## CMPR Block Ordering and T-Functions

CMPR chaining propagates **from upstream blocks into downstream ones** (high index to low index). Block 0 is the "source" with no chaining inputs; successive blocks layer chaining on top. `CMPR.blocks` returns bit indices ordered high-to-low.

PyPR's `TFunction` uses bit $n-1$ as the LSB, opposite the usual little-endian counter convention. This aligns the T-function's triangular dependency order with CMPR's high-to-low chaining flow. `TFunction` is implemented as a subclass of `CMPR` with one component per bit, so with this alignment the two share their analysis code.

## Vocabulary

- **Primitive polynomial** is the default noun for the polynomial associated with a register. "Primal" and "dual" are used only when disambiguating which relationship a polynomial has to a generated sequence.
- **Shift direction** is described as "shift up" / "shift down" in user-facing descriptions, and as "values progress toward bits with lower (higher) indices" in technical sections.
- **Frame** is the length-$m$ object: the register state, equivalently the sliding substring $\text{seq}[t..t+m-1]$ under the identification of states with windows of the sequence. **Window** is reserved for the length-$(m-1)$ object that appears once the discarded bit is moved to the other side of the generation equation — the input to $B$ in $f(s) = s_0 \oplus B$.

  The two differ by exactly the bit the shift discards, which is the bit every invertibility argument turns on, so the distinction is worth a word each. "Window" is never used for the state, even though "sliding window" is the natural phrase for it and the state genuinely is one.

  Where a substring of some other stated length is meant, it is qualified — "an $L$-window", "length-3 windows" — rather than left as a bare "window".

- **2-adic complexity** is the number of value cells $d$ in the smallest FCSR generating a sequence — the `diadic_complexity` parameter of `FCSR.__init__`. `FCSR_size` computes such a count from a fraction, though the count it returns is sufficient rather than least. It is the FCSR counterpart of an LFSR's linear complexity.
- **Size of the fraction** is $\Phi(p, q) = \max(|p|, |q|)$ — `phi` in `fcsrSynthesis.py`, and the quantity `BM_FCSR` compares at each discrepancy. It is not the 2-adic complexity: the two differ by roughly a logarithm, $d$ being about $\log_2 \Phi$.

  The FCSR literature often attaches "2-adic complexity", and "2-adic span", to a quantity on the $\log_2 \Phi$ scale rather than to a cell count. In this library the name always means the cell count $d$; a figure quoted across that boundary should state which of the two scales it is on. See [2-adic Integers and Rational Sequences](../theory/2-adic%20Integers%20and%20Rational%20Sequences.md) §Two measures of size for the full treatment.
