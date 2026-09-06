# Resolvent Analysis for Chaining Propagation

## What problem does this solve?

A CMPR has multiple MPR blocks. The "lower" blocks (higher bit indices) run independently — they are plain MPRs. The "upper" blocks receive **chaining input** from lower blocks: extra nonlinear terms added to their feedback functions. We want to understand how a disturbance or signal in the lower blocks propagates through to the output of an upper block.

## Setup

Let $B[t]$ be the state vector of one upper block at time $t$, and let $U$ be its linear update matrix (the underlying MPR's companion matrix). Let $\mathcal{C}(A[t])$ be the chaining input, a function of the lower block state $A[t]$. The block's recurrence is:

$$B[t+1] = UB[t] \oplus \mathcal{C}(A[t])$$

This is a linear system driven by an external input. We want to solve for $B$ as a function of $C$ and the initial state $B[0]$.

## Derivation: Z-Transform Convention

The choice between D-transform and Z-transform conventions mirrors the primal/dual polynomial distinction discussed in [Polynomial Conventions](../conventions/Polynomial%20Conventions.md). In the z-transform convention, a one-step delay is $z^{-1}$. If $B(z) = \sum_t B[t] z^{-t}$, then a one-step advance is multiplication by $z$.

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

- **D-convention:** $(I \oplus UD)^{-1}$ — this is the default in the library.
- **Z-convention:** $(D^{-1}I \oplus U)^{-1} = (zI \oplus U)^{-1}$ — this is the alternate flag.

The two matrices are the same up to a factor of $D$ (equivalently $z^{-1}$): the D-convention resolvent is $D$ times the Z-convention resolvent. The choice of convention is cosmetic and only matters for matching the notation of a particular proof you are trying to generate an example for.

The entries of these matrices are elements of $\text{GF}(2)[[D]]$ — formal power series over GF(2). In code, these entries are `SequenceTransform` objects. Their eigenstructure feeds the root/multiplicity calculation described in [Root Multiplicities and Jordan Decomposition](Root%20Multiplicities%20and%20Jordan%20Decomposition.md).

## The Propagation Matrix

The **propagation matrix** is a binary mask derived from the resolvent: entry $(i, j)$ is 1 if the resolvent entry $(i, j)$ is nonzero, and 0 otherwise. A 1 in position $(i, j)$ means that bit $j$'s chaining input eventually affects bit $i$ within the same block. The resolvent's eigenstructure feeds into the [root expression](Root%20Expressions%20and%20LC%20Estimation.md) computation, which bounds the linear complexity of the propagated signal, and the propagation pattern determines which [monomial profiles](Monomial%20Profile%20Theory.md) grow nontrivially as disturbances propagate through the chaining.

In practice, the propagation matrix is observed to be all-1s for any reasonable chaining — i.e., every bit affects every other bit. There is no formal proof that this is always the case, but it holds in all computed examples.

**Known limitation:** The propagation matrix currently cannot be computed in code because `SequenceTransform.__eq__` raises a `ValueError` when compared to an integer (the `!= 0` test used to build the mask fails). This is a known bug.

## API

- `C.update_matrices` — list of update matrices $U$, one per block.
- `C.resolvent_matrices` — list of resolvent matrices $(I \oplus UD)^{-1}$, one per block. Pass `alternate=True` for the Z-convention.
- `C.propagation_matrices` — intended binary masks, but currently broken (see above).
