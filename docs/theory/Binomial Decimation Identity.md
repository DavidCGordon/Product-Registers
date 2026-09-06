# Binomial Decimation Identity

For positive integers $d, k$ and any non-negative integer $n$, the binomial coefficient $\binom{dn}{k}$ can be decomposed modulo 2 into a linear combination (over $\mathbb{F}_2$) of binomial coefficients in $n$:

$$\binom{dn}{k} \;\equiv\; \sum_{i=0}^{k}\; \left(\bigoplus_{j \,\preccurlyeq\, i}\; \binom{dj}{k}\right) \binom{n}{i} \pmod{2}$$

where $\preccurlyeq$ denotes the **bit-subset partial order**: $j \preccurlyeq i$ when every bit that is set in the binary representation of $j$ is also set in $i$ (equivalently, $j \mathop{\&} i = j$). By Lucas' theorem, this is the same condition as $\binom{i}{j} \equiv 1 \pmod{2}$. See [Möbius Inversion on the Bit-Subset Lattice](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md) for a full treatment of $\preccurlyeq$ and the transform used in the proof.

## Proof

### Step 1: Define the function to decompose

Fix positive integers $d$ and $k$, and define

$$f(n) \;=\; \binom{dn}{k} \pmod{2}$$

This is a function from non-negative integers to $\{0,1\}$: it returns the parity of the binomial coefficient $\binom{dn}{k}$. By Lucas' theorem, $\binom{dn}{k} \equiv 1 \pmod{2}$ if and only if $k \preccurlyeq dn$ — every bit set in $k$ must also be set in the product $dn$. So $f(n)$ asks a question about the binary digits of $dn$: do they cover all the bits of $k$?

### Step 2: $f$ is a Boolean function of finitely many bits of $n$

Although $n$ can be arbitrarily large, $f(n)$ depends on only finitely many of its bits. The reason is that the predicate $k \preccurlyeq dn$ inspects only bits $0$ through $\lfloor \log_2 k \rfloor$ of $dn$ (since $k$ has no bits above that position), and integer multiplication carries only upward: bit $p$ of $dn$ is determined by bits $0, 1, \ldots, p$ of $n$. So $f(n)$ depends on at most bits $0$ through $\lfloor \log_2 k \rfloor$ of $n$.

This makes $f$ a Boolean function of finitely many input bits, which means it has a unique [Algebraic Normal Form](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md#algebraic-normal-form). The ANF writes $f$ as a sum of binomial-coefficient monomials: there exist unique coefficients $a_i \in \{0,1\}$ such that

$$f(n) \;=\; \bigoplus_{i}\; a_i\,\binom{n}{i} \pmod{2}$$

holds for every $n \geq 0$. Each $\binom{n}{i} \pmod{2}$ is a monomial in the bits of $n$ — it equals 1 when $n$ has all the bits of $i$ set, and 0 otherwise — and the ANF is the unique XOR-combination of these monomials that reproduces $f$.

### Step 3: Möbius inversion gives the coefficients

The ANF representation says that $f(n) = \bigoplus_{i \preccurlyeq n} a_i$: evaluating $f$ at $n$ amounts to XORing together all coefficients $a_i$ for which $i \preccurlyeq n$. The [Möbius transform](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md#the-möbius-transform) inverts this relationship — to recover each $a_i$, XOR together all values of $f$ at bit-subsets of $i$:

$$a_i \;=\; \bigoplus_{j \,\preccurlyeq\, i}\; f(j) \;=\; \bigoplus_{j \,\preccurlyeq\, i}\; \binom{dj}{k} \pmod{2}$$

The correctness of this inversion rests on a cancellation argument: substituting $f(j) = \bigoplus_{l \preccurlyeq j} a_l$ into the right-hand side, each $a_l$ with $l \neq i$ appears an even number of times (once for each element of the interval $\{j : l \preccurlyeq j \preccurlyeq i\}$, which has $2^{\text{popcount}(i) - \text{popcount}(l)} \geq 2$ elements) and cancels over $\mathbb{F}_2$. Only $a_i$ appears exactly once and survives. See the [proof of correctness](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md#proof-of-correctness) for the full argument.

### Step 4: The sum truncates at $k$

Over $\mathbb{Q}$, $\binom{dn}{k}$ is a polynomial of degree $k$ in $n$, and the binomial coefficients $\binom{n}{0}, \binom{n}{1}, \ldots$ form a basis for all polynomials (with $\binom{n}{i}$ contributing degree $i$). So there exist unique integers $c_i$ with

$$\binom{dn}{k} \;=\; \sum_{i=0}^{k}\; c_i\,\binom{n}{i}$$

as an identity in $\mathbb{Z}[n]$, and $c_i = 0$ for $i > k$ because $\binom{dn}{k}$ has degree exactly $k$. Reducing mod 2, the coefficients $c_i \bmod 2$ satisfy the same mod-2 identity as the ANF coefficients $a_i$. Since the ANF representation is unique, $a_i = c_i \bmod 2$, and in particular $a_i = 0$ for all $i > k$. $\square$

## Connections

- **[Möbius Inversion on the Bit-Subset Lattice](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md):** The transform machinery underlying Step 3.
- **[Decimation, Expansion, and Linear Complexity](Decimation%20Expansion%20and%20Linear%20Complexity.md):** Decimation by $d$ maps roots via $\alpha^{p/q} \mapsto \alpha^{dp/q}$; the Jordan block behavior under odd decimation involves $\binom{dt}{r} \pmod{2}$, which is exactly the predicate this identity decomposes.
