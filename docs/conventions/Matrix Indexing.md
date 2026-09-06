# Matrix Indexing and Vector Orientation

This document fixes the matrix and vector convention used throughout PyPR. It is deliberately concrete: displayed matrices use the same coordinates as NumPy arrays.

## Matrix Coordinates

For an $r \times c$ matrix $M$:

- `M[i, j]` means row `i`, column `j`, with both indices zero-based.
- Rows are displayed from top to bottom: row 0 is the top row.
- Columns are displayed from left to right: column 0 is the leftmost column.
- The first coordinate is the vertical coordinate; the second is the horizontal coordinate.

Thus the displayed entry in row `i` and column `j` is exactly the value returned by `M[i, j]` in NumPy. Do not silently reinterpret a displayed matrix as if its rows were written bottom-to-top.

## Column Vectors

Vectors are columns, and their coordinates use the same top-to-bottom order:

$$
\mathbf{x} =
\begin{pmatrix}
 x_0\\x_1\\\vdots\\x_{n-1}
\end{pmatrix}.
$$

A linear update is written on the left as

$$
\mathbf{y} = M\mathbf{x}, \qquad y_i = \bigoplus_j M_{i,j}x_j.
$$

Consequently, row `i` describes the computation of output coordinate `i`, while column `j` records where input coordinate `j` contributes. In particular, a register update matrix has the contract

$$
\texttt{next\_state} = M\,\texttt{current\_state},
$$

with `M[i, j] = 1` exactly when the update for state coordinate `i` reads coordinate `j`. This is also the contract of `FeedbackFunction.update_matrix`.

For a one-bit vector, for example,

$$
\begin{pmatrix}1&0\\1&1\end{pmatrix}
\begin{pmatrix}a_0\\a_1\end{pmatrix}
=
\begin{pmatrix}a_0\\a_0\mathbin\oplus a_1\end{pmatrix}.
$$

The lower-triangular matrix therefore performs the downward zeta sum on a column vector indexed by subsets in increasing integer order. Its transpose performs the upward zeta sum:

$$
\begin{pmatrix}1&1\\0&1\end{pmatrix}
\begin{pmatrix}v_0\\v_1\end{pmatrix}
=
\begin{pmatrix}v_0\mathbin\oplus v_1\\v_1\end{pmatrix}.
$$

## Binary and Bit Indices

When an integer index represents a bit-vector, PyPR uses

$$
 n = \sum_{b=0}^{m-1} n_b 2^b.
$$

Thus bit 0 is the least-significant binary bit, but coordinate 0 is still the top entry of a column vector. These are different ideas: numerical significance does not reverse NumPy's row order. The natural vector ordering is therefore

$$
0,1,2,\ldots,2^m-1,
$$

and entry `i` corresponds to the bit-set `bits(i)` wherever the surrounding construction uses integer-indexed subsets.

For a register state, coordinate `i` means state bit `i`; it is not implicitly the bottom row merely because bit 0 is the least-significant bit. In particular, a shift toward lower bit indices moves values upward in the displayed column.

## Reversal and Complementation

The exchange matrix

$$
J = \begin{pmatrix}0&1\\1&0\end{pmatrix}
$$
acts on a column vector by swapping coordinates. At scale $m$, $J^{\otimes m}$ maps integer index $i$ to

$$
\bar{i} = (2^m-1)-i,
$$

which is bitwise complement under the natural integer ordering. This statement depends on the row and column convention above; transposing or reversing the displayed vector changes which operation a written matrix appears to perform.

## Related Conventions

- [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md) uses rows for evaluation indices and columns for ANF coefficients, so its transform satisfies $\mathbf{f}=T\mathbf{a}$ under this convention.
- [Polynomial Conventions](Polynomial%20Conventions.md) uses the same column-vector orientation for Fibonacci, Galois, and MPR update matrices. Its “shift up” and “shift down” descriptions refer to bit indices, not to an unstated reversal of displayed rows.
- [Notation and Terminology](Notation%20and%20Terminology.md) records the library's coefficient-list and bit-order vocabulary.
