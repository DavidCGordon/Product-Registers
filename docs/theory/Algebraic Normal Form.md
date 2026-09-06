# Algebraic Normal Form

The Algebraic Normal Form (ANF) is the canonical way to represent a Boolean function as a polynomial over $\mathbb{F}_2$. This document develops the definition from first principles, connects it to binomial coefficients via Lucas' theorem, and shows that the truth-table-to-ANF conversion is a matrix operation — specifically, multiplication by Pascal's triangle mod 2, which decomposes as a Kronecker power of a single $2 \times 2$ matrix. Matrix entries and column-vector orientation follow [Matrix Indexing](../conventions/Matrix%20Indexing.md).

## Boolean Functions

A **Boolean function** on $m$ variables is a function $f: \mathbb{F}_2^m \to \mathbb{F}_2$ — it takes $m$ binary inputs and returns a binary output. We freely identify an input vector $(x_0, x_1, \ldots, x_{m-1}) \in \{0,1\}^m$ with the integer $n = \sum_{b=0}^{m-1} x_b \, 2^b$, so $f$ becomes a function from $\{0, 1, \ldots, 2^m - 1\}$ to $\{0, 1\}$.

The **truth table** of $f$ is the vector of all $2^m$ output values, indexed by $n$:

$$\mathbf{f} = \bigl(f(0),\; f(1),\; f(2),\; \ldots,\; f(2^m - 1)\bigr)$$

There are $2^{2^m}$ Boolean functions on $m$ variables — one for each possible truth table.

## Multilinear Polynomials over $\mathbb{F}_2$

Over any field, a function from a finite set to the field can be represented by a polynomial. Over $\mathbb{F}_2$, the constraint $x^2 = x$ (every element of $\mathbb{F}_2$ is idempotent) forces a strong simplification: any polynomial in $x_0, \ldots, x_{m-1}$ can be reduced so that each variable appears with exponent at most 1. The result is a **multilinear** polynomial — a sum of terms, each a product of a distinct subset of the variables.

Concretely, for each subset $S \subseteq \{0, \ldots, m-1\}$, define the **monomial**

$$\mu_S(x_0, \ldots, x_{m-1}) = \prod_{b \in S} x_b$$

(When $S = \emptyset$, the empty product is 1 — the constant monomial.) There are $2^m$ such monomials, one per subset.

**Definition.** The **Algebraic Normal Form** of a Boolean function $f$ is the unique representation

$$f(x_0, \ldots, x_{m-1}) = \bigoplus_{S \subseteq \{0,\ldots,m-1\}} a_S \prod_{b \in S} x_b$$

with coefficients $a_S \in \{0,1\}$.

**Uniqueness.** There are $2^m$ monomials and $2^{2^m}$ possible coefficient vectors $(a_S)$, giving $2^{2^m}$ candidate polynomials — the same count as the number of Boolean functions. So uniqueness follows if the monomials are linearly independent (over $\mathbb{F}_2$) as functions. We prove this below via the evaluation matrix.

The **algebraic degree** of $f$ is the size of the largest $S$ with $a_S = 1$ — the highest-order monomial that appears.

## The Binomial Basis

### Encoding subsets as integers

Each subset $S \subseteq \{0, \ldots, m-1\}$ corresponds to a unique integer $i = \sum_{b \in S} 2^b$ in $\{0, \ldots, 2^m - 1\}$. Under this encoding, we can index monomials by integers instead of sets: the monomial $\mu_i$ is the product of the variables whose bit positions are set in $i$.

### Monomials are binomial coefficients mod 2

What does the monomial $\mu_i$ evaluate to on input $n$? It computes $\prod_{b \in \text{bits}(i)} n_b$, which is 1 if and only if every bit set in $i$ is also set in $n$ — the bit-subset condition $i \preccurlyeq n$.

Lucas' theorem provides a second characterization of the same condition:

$$\binom{n}{i} \equiv \prod_{b \geq 0} \binom{n_b}{i_b} \pmod{2}$$

Each factor $\binom{n_b}{i_b}$ is 0 only when $i_b = 1$ and $n_b = 0$. So $\binom{n}{i} \equiv 1 \pmod{2}$ exactly when $i \preccurlyeq n$, giving the identification:

$$\mu_i(n) \;=\; \prod_{b \in \text{bits}(i)} n_b \;=\; \binom{n}{i} \pmod{2}$$

The monomial indexed by $i$ and the binomial coefficient $\binom{n}{i} \pmod{2}$ compute the same function of $n$.

### The ANF in binomial form

Substituting this into the ANF definition, the representation becomes:

$$\boxed{f(n) \;=\; \bigoplus_{i=0}^{2^m - 1}\; a_i\,\binom{n}{i} \pmod{2}}$$

The binomial coefficients $\binom{n}{0}, \binom{n}{1}, \ldots, \binom{n}{2^m - 1}$, read mod 2, are a basis for the space of Boolean functions on $m$ bits. The ANF coefficients $a_i$ are the coordinates in this basis.

## The Transform Matrix

### Evaluation as matrix multiplication

Evaluating the ANF at every input $n$ is a matrix-vector product. Define the $2^m \times 2^m$ matrix $T$ by

$$T_{n,i} \;=\; \binom{n}{i} \pmod{2}$$

Then the truth table $\mathbf{f}$ and the coefficient vector $\mathbf{a}$ are related by

$$\mathbf{f} \;=\; T\,\mathbf{a}$$

(Each row $n$ of $T$ dots with $\mathbf{a}$ to give $f(n) = \bigoplus_i a_i \binom{n}{i}$.)

For $m = 2$, the matrix is:

$$T = \begin{pmatrix}
\binom{0}{0} & \binom{0}{1} & \binom{0}{2} & \binom{0}{3} \\[4pt]
\binom{1}{0} & \binom{1}{1} & \binom{1}{2} & \binom{1}{3} \\[4pt]
\binom{2}{0} & \binom{2}{1} & \binom{2}{2} & \binom{2}{3} \\[4pt]
\binom{3}{0} & \binom{3}{1} & \binom{3}{2} & \binom{3}{3}
\end{pmatrix} \equiv \begin{pmatrix}
1 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 0 & 1 & 0 \\
1 & 1 & 1 & 1
\end{pmatrix} \pmod{2}$$

This is lower-triangular with 1s on the diagonal (since $\binom{n}{i} = 0$ for $i > n$ and $\binom{n}{n} = 1$), confirming the monomials are linearly independent and the ANF is unique.

### The one-bit matrix

For a single variable ($m = 1$), the matrix is

$$T_1 = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$$

This encodes the two evaluations: $f(0) = a_0$ (the constant term) and $f(1) = a_0 + a_1$ (constant plus linear term).

### The Kronecker product

The full $m$-bit matrix is the $m$-fold **Kronecker product** (tensor product) of the one-bit matrix:

$$T \;=\; T_1^{\otimes m} \;=\; \underbrace{T_1 \otimes T_1 \otimes \cdots \otimes T_1}_{m \text{ copies}}$$

This factorization follows from the bit-by-bit nature of the Lucas characterization. The entry $T_{n,i} = \binom{n}{i} \pmod{2}$ factors as:

$$\binom{n}{i} \;\equiv\; \prod_{b=0}^{m-1} \binom{n_b}{i_b} \;=\; \prod_{b=0}^{m-1} (T_1)_{n_b,\, i_b} \pmod{2}$$

This per-coordinate factorization is precisely the defining property of the Kronecker product: $(A \otimes B)_{(r_1 r_2),\,(c_1 c_2)} = A_{r_1,c_1} \cdot B_{r_2,c_2}$, extended to $m$ factors.

## Pascal's Triangle and the Sierpiński Pattern

The matrix $T$ whose entries are $\binom{n}{i} \pmod{2}$ is Pascal's triangle read mod 2. The pattern of 1s — the positions where $\binom{n}{i}$ is odd — forms the **Sierpiński triangle**, the well-known fractal.

The Kronecker product explains the self-similarity. The recursion $T^{\otimes (m+1)} = T_1 \otimes T^{\otimes m}$ expands as:

$$T^{\otimes(m+1)} = \begin{pmatrix} 1 \cdot T^{\otimes m} & 0 \cdot T^{\otimes m} \\[4pt] 1 \cdot T^{\otimes m} & 1 \cdot T^{\otimes m} \end{pmatrix} = \begin{pmatrix} T^{\otimes m} & 0 \\[4pt] T^{\otimes m} & T^{\otimes m} \end{pmatrix}$$

Each doubling of the number of bits produces a matrix with three copies of the previous level arranged in a triangle — the defining construction of the Sierpiński fractal. The top-right block is always zero because $\binom{n}{i} = 0$ when $i > n$, and the three nonzero blocks arise from the three nonzero entries of $T_1$.

## The Transform is Self-Inverse

The one-bit matrix satisfies $T_1^2 = I$ over $\mathbb{F}_2$:

$$\begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}^2 = \begin{pmatrix} 1 & 0 \\ 1+1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} \pmod{2}$$

Since the Kronecker product respects powers — $(A^{\otimes m})^2 = (A^2)^{\otimes m}$ — the full transform is also self-inverse:

$$(T^{\otimes m})^2 = (T_1^2)^{\otimes m} = I^{\otimes m} = I$$

This means the **same matrix** converts in both directions:

$$\mathbf{f} = T\,\mathbf{a} \qquad\text{and}\qquad \mathbf{a} = T\,\mathbf{f}$$

Written out entry-by-entry, the inversion direction is:

$$a_i = \bigoplus_{j \,\preccurlyeq\, i}\; f(j)$$

which is the [Möbius inversion formula](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md#the-möbius-transform). The self-inverse property explains it from the matrix side: the evaluation matrix (ANF $\to$ truth table) and the Möbius matrix (truth table $\to$ ANF) are the same object over $\mathbb{F}_2$, because $+1 = -1$ makes the zeta and Möbius functions identical. Over $\mathbb{Z}$, the Möbius inverse of Pascal's triangle has alternating signs; over $\mathbb{F}_2$, those signs collapse and the matrix is its own inverse.

## Efficient Computation: the Butterfly Decomposition

The Kronecker factorization $T = T_1 \otimes T_1 \otimes \cdots \otimes T_1$ means $T$ can be applied as $m$ independent passes, one per bit position, without ever forming the full $2^m \times 2^m$ matrix.

**One pass for bit $b$:** for every integer $i$ with bit $b$ set, XOR the value at $i$ with the value at $i$ with bit $b$ cleared:

$$a[i] \;\leftarrow\; a[i] \oplus a[i \oplus 2^b]$$

**The full transform** applies this pass for each $b = 0, 1, \ldots, m-1$. Each pass visits $2^{m-1}$ pairs, so the total work is $m \cdot 2^{m-1} = O(m \cdot 2^m)$ — compared to $O(4^m)$ for naive matrix multiplication.

This is the Boolean analog of the FFT butterfly: each Kronecker factor contributes one round of pairwise operations, and the rounds compose to apply the full matrix. The structure is identical to what competitive programmers call the "sum over subsets" (SOS) transform, and what hardware designers recognize as the building block of Reed-Muller encoders.

**Example** ($m = 2$). Start with truth table $\mathbf{f} = (0, 1, 0, 0)$.

Pass for bit 0 ($b = 0$): XOR each odd-indexed entry with the entry below it.
- $a[1] \leftarrow a[1] \oplus a[0] = 1 \oplus 0 = 1$
- $a[3] \leftarrow a[3] \oplus a[2] = 0 \oplus 0 = 0$
- State: $(0, 1, 0, 0)$

Pass for bit 1 ($b = 1$): XOR each entry with bit 1 set with the entry with bit 1 cleared.
- $a[2] \leftarrow a[2] \oplus a[0] = 0 \oplus 0 = 0$
- $a[3] \leftarrow a[3] \oplus a[1] = 0 \oplus 1 = 1$
- State: $(0, 1, 0, 1)$

Result: $\mathbf{a} = (0, 1, 0, 1)$, so $f(n) = \binom{n}{1} + \binom{n}{3} \pmod{2}$ — matching the [example](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md#example) in the Möbius inversion doc.

## Connections

- **[Möbius Inversion on the Bit-Subset Lattice](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md):** The lattice-theoretic perspective on the same transform — proves the inversion formula via interval counting rather than matrix algebra.
- **[Binomial Decimation Identity](Binomial%20Decimation%20Identity.md):** Uses the ANF decomposition to express $\binom{dn}{k} \pmod{2}$ in the binomial basis.
- **[Algebraic Attacks](Algebraic%20Attacks.md):** ANF is the representation underlying algebraic cryptanalysis — annihilator computation, equation generation, and low-degree relation search all operate in this ring.
