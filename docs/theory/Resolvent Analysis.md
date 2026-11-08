# Resolvent Analysis for Chaining Propagation

This document derives the resolvent matrix of a CMPR block: the matrix of formal power series that expresses a chained block's state sequence in terms of its chaining input and its initial state. It gives the derivation in both the D-transform and the Z-transform conventions, relates the two resulting matrices, and defines the propagation matrix obtained from the resolvent. It assumes linear algebra over $\mathbb{F}_2$ and generating functions of sequences; the conventions for the latter are in [Polynomial Conventions](../conventions/Polynomial%20Conventions.md).

## Chaining Between Blocks

A CMPR has multiple MPR blocks. Block 0, the source, receives no chaining input and runs independently as a plain MPR. Every other block receives **chaining input** from the blocks upstream of it: extra nonlinear terms added to its feedback functions. Chaining runs from blocks with higher bit indices into blocks with lower bit indices (see [Notation and Terminology](../conventions/Notation%20and%20Terminology.md)). The question this document answers is how a disturbance or signal in the upstream blocks propagates into the state sequence of a block that receives chaining input.

## Setup

Let $B[t]$ be the state vector of one chained block at time $t$, and let $U$ be its linear update matrix (the underlying MPR's companion matrix). Let $\mathcal{C}(A[t])$ be the chaining input, a function of the upstream state $A[t]$. The block's recurrence is:

$$B[t+1] = UB[t] \oplus \mathcal{C}(A[t])$$

This is a linear system driven by an external input. We want to solve for $B$ as a function of the chaining input and the initial state $B[0]$. In the derivations below, $C(z)$ and $C(D)$ denote the transform of the chaining-input sequence $t \mapsto \mathcal{C}(A[t])$ in the respective convention.

## Derivation: Z-Transform Convention

The choice between D-transform and Z-transform conventions corresponds to the primal/dual polynomial distinction of [Polynomial Conventions](../conventions/Polynomial%20Conventions.md): under the D-transform the dual polynomial of a sequence convolves it to zero, and under the Z-transform the primal polynomial does. In the z-transform convention, a one-step delay is $z^{-1}$. If $B(z) = \sum_t B[t] z^{-t}$, then a one-step advance is multiplication by $z$.

Taking the z-transform of the recurrence $B[t+1] = UB[t] \oplus \mathcal{C}(A[t])$:

$$
\begin{align*}
z(B(z) \oplus B[0]) &= UB(z) \oplus C(z) \\
zB(z) \oplus UB(z) &= C(z) \oplus zB[0] \\
(zI \oplus U)B(z) &= C(z) \oplus zB[0] \\
B(z) &= (zI \oplus U)^{-1}(C(z) \oplus zB[0])
\end{align*}
$$

## Derivation: D-Transform Convention

In the D-transform convention (common in coding theory / algebraic approaches), a one-step delay is $D$. The generating function is $B(D) = \sum_t B[t] D^t$, so advancing by one step means dividing by $D$ (i.e., multiplying by $D^{-1}$).

Applying $D^{-1}$ to both sides of $B[t+1] = UB[t] \oplus \mathcal{C}(A[t])$:

$$
\begin{align*}
D^{-1}(B(D) \oplus B[0]) &= UB(D) \oplus C(D) \\
B(D) \oplus B[0] &= UDB(D) \oplus DC(D) \\
(I \oplus UD)B(D) &= DC(D) \oplus B[0] \\
B(D) &= (I \oplus UD)^{-1}(DC(D) \oplus B[0])
\end{align*}
$$

Alternatively, without multiplying through by $D$ first:

$$
\begin{align*}
D^{-1}B(D) \oplus UB(D) &= C(D) \oplus D^{-1}B[0] \\
(D^{-1}I \oplus U)B(D) &= C(D) \oplus D^{-1}B[0] \\
B(D) &= (D^{-1}I \oplus U)^{-1}(C(D) \oplus D^{-1}B[0])
\end{align*}
$$

## The Two Resolvent Conventions

The two derivations give two related resolvent matrices:

- **D-convention:** $(I \oplus UD)^{-1}$. This is the matrix the library computes (`CMPR.resolvent_matrices`).
- **Z-convention:** $(D^{-1}I \oplus U)^{-1} = (zI \oplus U)^{-1}$. The library produces this form only in `generate_resolvent_example(..., use_z_convention=True)`, described under [API](#api).

The two differ by a factor of $D$. Since

$$D^{-1}I \oplus U = D^{-1}(I \oplus UD),$$

inverting both sides gives

$$(D^{-1}I \oplus U)^{-1} = D\,(I \oplus UD)^{-1}:$$

the Z-convention resolvent is $D$ (equivalently $z^{-1}$) times the D-convention resolvent. The vectors they multiply differ by the inverse factor, $C(D) \oplus D^{-1}B[0] = D^{-1}\bigl(DC(D) \oplus B[0]\bigr)$, so both products are the same $B(D)$. The choice of convention is cosmetic; it matters only for matching the notation of a particular proof when generating an example for it.

The entries of these matrices are elements of $\text{GF}(2)[[D]]$ — formal power series over GF(2). In code, the entries are `BooleanGF` objects (`PyPR.BooleanLogic.BooleanGF`), which represent them as quotients of polynomials in $D$. Their eigenstructure feeds the root/multiplicity calculation described in [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md).

## The Propagation Matrix

The **propagation matrix** is a binary mask derived from the resolvent: entry $(i, j)$ is 1 if the resolvent entry $(i, j)$ is nonzero, and 0 otherwise. A 1 in position $(i, j)$ means that bit $j$'s chaining input eventually affects bit $i$ within the same block. The resolvent's eigenstructure feeds into the [root expression](Root%20Expressions%20and%20LC%20Estimation.md) computation, which bounds the linear complexity of the propagated signal, and the propagation pattern determines which [monomial profiles](Monomial%20Profile%20Theory.md) grow nontrivially as disturbances propagate through the chaining.

The resolvent $(I \oplus UD)^{-1}$ depends only on the block's update matrix $U$, so the propagation matrix does not depend on the chaining functions.

**Proposition.** For an MPR block with the default update, in which $U$ is multiplication by $x$ on $\mathbb{F}_2[x]/P$, the propagation matrix is all ones: every bit of the block affects every other bit.

*Proof.* Expanding the inverse as a power series,

$$(I \oplus UD)^{-1} = \sum_{t \geq 0} U^t D^t,$$

so the resolvent entry $(i, j)$ is nonzero if and only if $(U^t)_{i,j} = 1$ for some $t \geq 0$. Identify the block state with $\mathbb{F}_2[x]/P$, so that the coordinate vector $e_j$ is the element $x^j$. Then $U^t e_j = x^{t+j}$. The element $x$ is invertible because $P(0) = 1$, which holds for the primitive polynomials that MPR blocks use, and it has finite multiplicative order $N$. For any $i$ and $j$, choose $t \geq 0$ with $t \equiv i - j \pmod{N}$; then $U^t e_j = x^{i} = e_i$, so $(U^t)_{i,j} = 1$. $\square$

**Other update polynomials.** An MPR built with a non-default `update_poly` has $U$ equal to multiplication by some other element $u$, and the propagation matrix need not be all ones. The first step of the proof still applies: entry $(i, j)$ is 1 exactly when $x^i$ appears in $u^t x^j$ for some $t \geq 0$. With $u = 1$, $U$ is the identity and the mask is the identity matrix. With $P = x^4 + x + 1$ and $u = x + x^3$ the mask has 13 ones out of 16. Both cases are asserted in `test/Tools/linear_complexity_test.py`.

`CMPR.propagation_matrices` computes the masks from the resolvent matrices, whose entries are `BooleanGF` fractions; an entry counts as nonzero when its numerator is nonzero.

## API

In the entries below, `C` is a `CMPR` instance.

- `C.update_matrices` — list of update matrices $U$, one per block.
- `C.resolvent_matrices` — list of resolvent matrices $(I \oplus UD)^{-1}$, one per block. It is a property and takes no arguments.
- `generate_resolvent_example(cmpr, use_z_convention=False)` in `PyPR.Tools.ResolventSolving` — `cmpr` is a `FeedbackRegister` whose feedback function is a CMPR. The function returns a dictionary of the transforms that appear in the derivations above. With `use_z_convention=True` it multiplies each resolvent matrix by $D$ and the combined input vector by $z = D^{-1}$, which is the Z-convention form.
- `C.propagation_matrices` — the support masks of the resolvent matrices, one per block, as `uint8` arrays.
