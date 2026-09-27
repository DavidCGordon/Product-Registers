# 2-adic Integers and Rational Sequences

The arithmetic an FCSR realizes. An LFSR computes in $\mathbb{F}_2[[x]]$, where
carries are discarded; an FCSR computes in $\mathbb{Z}_2$, where they are kept
and propagated. This document covers the 2-adic side: which sequences are
rational, how clocking acts on the fraction, and which fractions give sequences
with no transient.

For how this is built in hardware terms — the interleaved carry cells, the tap
derivation, the encoding of a fraction as a bit pattern — see
[FCSR Implementation](../architecture/FCSR%20Implementation.md).

## The 2-adic integers

A **2-adic integer** is a formal sum $\sum_{i \geq 0} s_i 2^i$ with
$s_i \in \{0,1\}$, added and multiplied with carries propagating upward. Unlike
$\mathbb{Z}$, these sums need not terminate, and $\mathbb{Z}_2$ is exactly the
set of all of them.

Every binary sequence is therefore a 2-adic integer, read low bit first, and the
correspondence is a bijection. The one fact that makes this useful is that
negative integers have infinite expansions:

$$-1 = 1 + 2 + 4 + 8 + \cdots = \ldots1111_2,$$

since adding $1$ carries forever and leaves $0$. More generally $-n$ is the
two's complement pattern extended to infinity.

Odd integers are invertible, and the inverse is forced one bit at a time.
Suppose $q$ is odd and bits $u_0, \ldots, u_{n-1}$ have been fixed so that
$q \sum_{i < n} u_i 2^i \equiv 1 \pmod{2^n}$; write
$q \sum_{i < n} u_i 2^i = 1 + 2^n k$. Appending $u_n$ adds $q u_n 2^n$ to the
product, so the congruence extends to $2^{n+1}$ exactly when $k + q u_n$ is
even, and since $q$ is odd that is $u_n \equiv k \pmod 2$ — one admissible
choice at every step. The base case is the same computation with the empty sum:
$q u_0 \equiv 1 \pmod 2$ forces $u_0 = 1$. The bits are therefore determined,
and $u = \sum_i u_i 2^i$ is the unique 2-adic integer with $q u = 1$. Even
integers have no inverse, since every product with an even number has low bit
$0$.

So a fraction $p/q$ in lowest terms names an element of $\mathbb{Z}_2$ — an $x$
with $qx = p$ — exactly when $q$ is odd. If $q$ is odd, $x = p q^{-1}$ is such
an element, and it is the only one. If $q$ is even, any such $x$ would give
$p = qx = 2 \cdot (q/2)x$, a product with an even number, so $p$ would have low
bit $0$; but $p$ and $q$ both even contradicts lowest terms. These fractions are
what FCSRs generate.

## Reading off the expansion

Given $p/q$ with $q$ odd, the bits come out one at a time. The low bit is forced:
$s_0 \equiv p q^{-1} \equiv p \pmod 2$, because $q$ is odd. Subtracting it and
dividing by two leaves another fraction with the same denominator:

$$\frac{p}{q} = s_0 + 2\cdot\frac{p'}{q}, \qquad p' = \frac{p - s_0 q}{2}.$$

Write $\sigma$ for that map on numerators:

$$\boxed{\;\sigma(p) = \frac{p - s_0 q}{2}, \qquad s_0 = p \bmod 2.\;}$$

Then the emitted sequence is $s_0, s_1, s_2, \ldots$ where $s_t$ is the low bit
of $\sigma^t(p)$. **This is what an FCSR does**: one clock is one application of
$\sigma$, and the bit it emits is the parity of the current numerator.

$\sigma$ is well defined — $p - s_0 q$ is even by the choice of $s_0$ — and it is
the only step the rest of this document needs.

## Which sequences are rational

> **Proposition.** A binary sequence is eventually periodic if and only if it is
> the expansion of some $p/q$ with $q$ odd.

Take $q > 0$ for the rest of this document. Negating $p$ and $q$ together leaves
the fraction unchanged, so every fraction with odd denominator has such a
representative, and positivity is needed below: the bound
$|\sigma(p)| \leq (|p| + q)/2$ comes from
$|p - s_0 q| \leq |p| + s_0 q \leq |p| + q$, whose last step assumes $q \geq 0$.

The forward direction is what *synthesis* rests on — `BM_FCSR` can only return a
fraction for an observed sequence because one exists — and it takes three lines.
Let the sequence have preperiod $k$ and period $T \geq 1$; write
$A = \sum_{i<k} s_i 2^i$ for the transient, $c = \sum_{i<T} s_{k+i} 2^i$ for one
period, and $B \in \mathbb{Z}_2$ for the value of the tail
$s_k, s_{k+1}, \ldots$, so that $S = A + 2^k B$. Periodicity of the tail says
that deleting its low $T$ bits returns it unchanged, $(B - c)/2^T = B$, i.e.
$B(2^T - 1) = -c$. Substituting,

$$S = A - \frac{2^k c}{2^T - 1} = \frac{A(2^T - 1) - 2^k c}{2^T - 1},$$

whose denominator $2^T - 1$ is odd.

The reverse direction is what *generation* rests on — a register state gives a
numerator, and what it emits is that fraction's expansion. It follows from
$\sigma$ in two steps, using $|\sigma(p)| \leq (|p| + q)/2$ throughout.

*The numerators enter $[-q, q]$ and stay there.* If $|p| > q$ then
$(|p| + q)/2 < |p|$, so the magnitude strictly decreases; being a non-negative
integer it cannot decrease forever, so some iterate has $|p| \leq q$. From then
on $(|p| + q)/2 \leq q$, so the bound is preserved.

*Therefore the sequence cycles.* Only finitely many integers lie in $[-q, q]$,
so two iterates coincide; $\sigma$ is a function, so the numerators repeat from
there, and with them the bits. $\square$

## Two measures of size: $d$ and $\Phi$

The **2-adic complexity** of a sequence is the number of value cells in the
smallest FCSR generating it — the $d$ of
[FCSR Implementation](../architecture/FCSR%20Implementation.md). This is the
part linear complexity plays for LFSRs: the size of the smallest generator in
the relevant class, an LFSR length there and a cell count here.

A note on units when comparing against the literature. Papers on FCSRs use
"2-adic complexity" and "2-adic span" for a quantity on the $\log_2 \Phi$ scale
rather than for a cell count, and they do not always say which they mean. The
two differ by a small additive constant — a register realizing $p/q$ carries
about $\log_2 \Phi$ value cells — so they agree on which sequences are hard and
disagree only on the number attached. Within this repository the name always
means the cell count $d$. When quoting a figure across that boundary, say which
one it is: `fcsrSynthesis.py` benchmarks against `HuSha2019`, and a comparison
is only meaningful with both sides measured the same way.

The quantity synthesis actually tracks is not $d$ but the size of the fraction,

$$\Phi(p, q) = \max(|p|, |q|),$$

which is `phi` in `fcsrSynthesis.py` and what `BM_FCSR` compares — at each
discrepancy it tests
`phi(num_curr, den_curr) < scale * phi(num_prev, den_prev)`.

`FCSR_size` is what the implementation computes when it needs a cell count for
a given fraction, and what it returns is an upper bound rather than the
minimum: a cell count *sufficient* to realize $p/q$ under `state_from_frac`'s
encoding, which need not be the least $d$ over all FCSRs generating the
sequence. `FCSR_size(-1, 1)` returns $2$,
while `FCSR(1, 1)` is a well-formed register with a single value cell whose
state `[1]` emits $1, 1, 1, \ldots$, the expansion of $-1/1$.

Minimizing $\Phi$ and minimizing `FCSR_size` select the same fraction, but only
*within one sequence*, and that is all the reduce-or-not convention needs.

1. Fix a sequence and let $p_0/q_0$ generate it, with $q_0 > 0$ odd and
   $\gcd(p_0, q_0) = 1$ — reduce any generating fraction to get one. Any other
   $p/q$ with $q > 0$ odd generating the same
   sequence names the same $x \in \mathbb{Z}_2$, so $p = qx$ and $p_0 = q_0x$
   give $pq_0 = qxq_0 = p_0q$ as integers; with $\gcd(p_0,q_0) = 1$ this forces
   $(p, q) = (mp_0, mq_0)$ for a positive integer $m$, odd because $q$ is.
2. $\Phi(mp_0, mq_0) = m\,\Phi_0$ with $\Phi_0 = \Phi(p_0,q_0) \geq q_0 \geq 1$,
   so $\Phi$ is strictly increasing in $m$.
3. `FCSR_size`$(mp_0, mq_0)$ is minimized at $m = 1$. Scaling by $m \geq 1$ does
   not change the sign of the numerator, so the same branch of `FCSR_size`
   applies at every $m \geq 1$ except where an explicitly special-cased fraction
   intervenes — and inside a branch every quantity built from $|mp_0|$ and
   $mq_0$ is non-decreasing in $m$: $\max(|mp_0|, mq_0)$ in the positive branch,
   $1 + \lceil\log_2(|mp_0|/3 + 1)\rceil$ and $\lceil \log_2 mq_0\rceil$ in the
   other, as are the maxima and logarithms applied to them.

   The special cases are $p/q \in \{0/1,\ 1/1\}$, where `FCSR_size` returns the
   cell count the register actually needs rather than the one the logarithm
   gives. Neither disturbs the conclusion. $0/1$ scales to $0/m$, which is not
   in lowest terms for $m > 1$, so $m = 1$ is the only case. For $1/1$ the
   special case returns $2$ at $m = 1$ while the positive branch governs
   $m \geq 3$, giving $3, 4, 4, 5$ at $m = 3, 5, 7, 9$ — so the minimum is still
   at $m = 1$, now because $2$ is below the branch's values rather than because
   one branch runs throughout.

4. Both are therefore minimized at $m = 1$, the reduced fraction — the same
   argmin.

This does **not** say that $\Phi$ and `FCSR_size` order *different* sequences
the same way. `FCSR_size` is not a function of $\Phi$ at all: its non-positive
branch depends on $|p|$ and $q$ separately rather than through $\max(|p|,q)$.
At $-70/1$, $\Phi = 70$ and `FCSR_size` $= 6$; at $-1/65$, $\Phi = 65$ and
`FCSR_size` $= 7$ — the larger $\Phi$ with the smaller cell count.

### How many bits pin the minimal fraction

Minimality is not a property one can check against a prefix without a bound on
how long the prefix must be — two different fractions can agree for a while. The
bound proved below is $N > 2\log_2 \Phi + 1$, which under the correspondence
$L \leftrightarrow \log_2 \Phi$ between an LFSR length and a fraction size is
Berlekamp-Massey's $2L$ requirement with one bit added. It comes from one
observation.

> **Lemma.** If $p_1/q_1$ and $p_2/q_2$ have odd denominators and their
> expansions agree on the first $N$ bits, then $2^N \mid p_1q_2 - p_2q_1$.
>
> *Proof.* Agreeing on the first $N$ bits says the two elements of
> $\mathbb{Z}_2$ are congruent mod $2^N$, i.e.
> $p_1q_1^{-1} \equiv p_2q_2^{-1} \pmod{2^N}$. Both $q_i$ are odd, hence units,
> so multiplying through by $q_1q_2$ preserves the congruence and gives
> $p_1q_2 \equiv p_2q_1 \pmod{2^N}$ in $\mathbb{Z}_2$. The difference
> $m = p_1q_2 - p_2q_1$ is an integer, and for an integer the congruence
> descends to $\mathbb{Z}$: write $m = 2^jw$ with $w$ odd (the case
> $m = 0$ being trivial); if $j < N$ then $m = 2^Nu$ with $u \in \mathbb{Z}_2$
> would give $w = 2^{N-j}u$ with $N - j \geq 1$, a product with an even number
> and so of low bit $0$, contradicting $w$ odd. Hence $2^N \mid m$ in
> $\mathbb{Z}$. $\square$

> **Proposition.** If in addition $2\,\Phi_1\Phi_2 < 2^N$, where
> $\Phi_i = \Phi(p_i, q_i)$, then $p_1q_2 = p_2q_1$ — the two are equal as
> rationals, hence name the same element of $\mathbb{Z}_2$ and the same
> sequence.
>
> *Proof.* $|p_1q_2 - p_2q_1| \leq |p_1||q_2| + |p_2||q_1| \leq 2\,\Phi_1\Phi_2
> < 2^N$. The lemma says $2^N$ divides that quantity, and the only multiple of
> $2^N$ smaller than $2^N$ in absolute value is $0$. So $p_1q_2 = p_2q_1$. $\square$

Equality as rationals is weaker than equality of the pairs $(p_i, q_i)$, and the
gap is not idle: $-1/3$ and $-3/9$ both expand to $1,0,1,0,\ldots$, satisfy the
hypothesis at $N = 8$ ($2\Phi_1\Phi_2 = 54 < 256$), and have
$p_1q_2 = p_2q_1 = -9$, yet $(-1,3) \neq (-3,9)$ — and a fraction and its
reduction name *different registers*, as
[FCSR Implementation](../architecture/FCSR%20Implementation.md) §Encoding a
fraction shows. Pinning the pair takes the extra hypothesis of the corollary.

> **Corollary.** Let $p/q$ generate the observed sequence, with $q > 0$ odd and
> $\Phi = \Phi(p,q)$ minimal among all such fractions generating it, and
> suppose $N > 1 + 2\log_2 \Phi$, so that $2\Phi^2 < 2^N$. Then $(p,q)$ is the
> **only** pair with $q' > 0$ odd and $\Phi' \leq \Phi$ whose expansion agrees
> with the sequence on $N$ bits.
>
> *Proof.* First, $\Phi$-minimality forces $\gcd(p,q) = 1$.
>
> 1. $\Phi(ma, mb) = m\,\Phi(a,b)$ for any integer $m \geq 1$, since
>    $\max(|ma|,|mb|) = m\max(|a|,|b|)$; so with $(a,b)$ fixed and
>    $\Phi(a,b) \geq 1$, $\Phi(ma,mb)$ is strictly increasing in $m$.
> 2. Let $p_0/q_0$ be the reduced form of $p/q$ with $q_0 > 0$, so
>    $(p,q) = (mp_0, mq_0)$ with $m = \gcd(p,q) \geq 1$, and $\Phi = m\Phi_0$
>    for $\Phi_0 = \Phi(p_0,q_0)$ by (1). The pair $(p_0,q_0)$ names the same
>    element of $\mathbb{Z}_2$ and so generates the same sequence, so
>    minimality of $\Phi$ gives $\Phi \leq \Phi_0$, i.e. $m\Phi_0 \leq \Phi_0$;
>    with $\Phi_0 \geq q_0 \geq 1$ this forces $m = 1$, so $p/q$ is in lowest
>    terms.
>
> Now take a competitor $p'/q'$ with $q' > 0$ odd, $\Phi' \leq \Phi$, agreeing
> on $N$ bits.
>
> 3. $2\Phi'\Phi \leq 2\Phi^2 < 2^N$, so the proposition gives $p'q = pq'$,
>    i.e. $p'/q' = p/q$ as rationals. With $\gcd(p,q) = 1$ from (2), $q \mid pq'$
>    forces $q \mid q'$; writing $q' = m'q$ and substituting gives $p' = m'p$,
>    with $m' \geq 1$ because $q, q' > 0$ and $m'$ odd because $q'$ is.
> 4. $\Phi' = \Phi(m'p, m'q) = m'\Phi \leq \Phi$ with $\Phi \geq q \geq 1$
>    forces $m' = 1$, so $(p', q') = (p, q)$ as pairs. $\square$

Step 2 is where the $q > 0$ convention is spent: without it $(p,q)$ and
$(-p,-q)$ would be distinct pairs of equal $\Phi$ generating the same sequence,
and no hypothesis would separate them.

So past $N > 2\log_2\Phi + 1$ observed bits the minimal fraction is pinned: a
search that returns a consistent fraction no larger than the true minimum has
necessarily returned that minimum, because there is nothing else for it to
return. That is the $2L$ leading term of Berlekamp-Massey's requirement with one
extra bit, and the extra bit is not slack in this argument: the proposition needs
$2\Phi_1\Phi_2 < 2^N$, and the factor $2$ there comes from
$|p_1q_2 - p_2q_1| \leq |p_1||q_2| + |p_2||q_1|$, two terms each of size up to
$\Phi_1\Phi_2$.

What the corollary does not supply is the remaining step — that `BM_FCSR`'s
approximation loop does return a fraction with $\Phi' \leq \Phi$ rather than
merely a consistent one. That is an invariant of the algorithm rather than of
the arithmetic.

For that invariant, `KlapperGoresky1993` is the reference: `bibliography.bib`
records it as the primary source for the FCSR synthesis algorithm and notes that
PyPR's 2-adic Berlekamp-Massey variant follows it, with `KlapperGoresky1997`
supplying later modifications. One caveat when reading the guarantee against
this implementation: `fcsrSynthesis.py`'s `D()` helper departs from the
published version and follows `ArnaultBergerMinier2005` instead.

## Purely periodic sequences

A sequence is **purely periodic** when it repeats from the very first bit, with
no transient. Among the rational sequences these are exactly the ones whose
numerator sits in a bounded window.

> **Proposition.** With $q > 0$ odd, the expansion of $p/q$ is purely periodic if
> and only if $-q \leq p \leq 0$.
>
> *Proof.* ($\Rightarrow$) Suppose the expansion repeats with period $T$, so
> shifting $T$ times returns the same 2-adic integer. Writing
> $c = \sum_{i<T} s_i 2^i \in [0, 2^T - 1]$ for the first period, that says
> $S = (S - c)/2^T$ with $S = p/q$. Rearranging, $S(2^T - 1) = -c$, so
> $S = -c/(2^T-1) \in [-1, 0]$, and multiplying by $q > 0$ gives
> $-q \leq p \leq 0$.
>
> ($\Leftarrow$) Suppose $-q \leq p \leq 0$. Since $q$ is odd, $2$ is a unit
> modulo $q$; let $T = \operatorname{ord}_q(2)$, so $q \mid 2^T - 1$, say
> $2^T - 1 = mq$ with $m \geq 1$. Put $c = -pm$. From $0 \leq -p \leq q$ we get
> $0 \leq c \leq qm = 2^T - 1$, so $c$ is a legitimate $T$-bit integer, and
>
> $$S = \frac{p}{q} = \frac{pm}{mq} = \frac{-c}{2^T - 1},$$
>
> so $S(2^T-1) = -c$, i.e. $S - c = 2^T S$. That equation says two things.
> First, $S \equiv c \pmod{2^T}$, and with $0 \leq c \leq 2^T - 1$ this makes
> $c$ exactly the integer formed by the low $T$ bits of $S$ — so deleting those
> bits is the operation $(S-c)/2^T$. Second, that operation returns $S$. The
> expansion therefore repeats from its first bit, with period dividing $T$.
> $\square$

Call $[-q, 0]$ the **periodic core**. It holds $q+1$ numerators, both endpoints
included: the $q$ integers $-q \leq p \leq -1$, plus $0$. Both endpoints are
fixed points of $\sigma$, since $\sigma(0) = 0$ and $\sigma(-q) = (-q - q)/2 =
-q$ using that $-q$ is odd.

> **Proposition.** $\sigma$ restricts to a bijection of the core onto itself.
>
> *Proof.* Closure: for $-q \leq p \leq 0$, $\sigma(p)$ is $p/2$ or $(p-q)/2$,
> both in $[-q, 0]$ given the range of $p$. Injectivity: if
> $\sigma(p_1) = \sigma(p_2)$ then $p_1 - s_1 q = p_2 - s_2 q$, so
> $p_1 - p_2 = (s_1 - s_2)q \in \{0, \pm q\}$. Take the three cases on the
> parities.
>
> - $s_1 = s_2$: then $p_1 - p_2 = 0$, which is the conclusion.
> - $s_1 = 1$, $s_2 = 0$: then $p_1 - p_2 = q > 0$, and both lie in $[-q, 0]$,
>   an interval of width exactly $q$, so the only possibility is $p_1 = 0$ and
>   $p_2 = -q$. But $p_1 = 0$ is even, so $s_1 = 0$ — contradicting $s_1 = 1$.
> - $s_1 = 0$, $s_2 = 1$: then $p_1 - p_2 = -q$, forcing $p_1 = -q$ and
>   $p_2 = 0$ by the same width argument. But $-q$ is odd, so $s_1 = 1$ —
>   again a contradiction.
>
> So only the first case survives. Equivalently and more briefly: the only pair
> of core numerators differing by $\pm q$ is $\{0, -q\}$, and $\sigma$ does not
> identify those two, since $\sigma(0) = 0 \neq -q = \sigma(-q)$ by the two
> fixed-point computations above. $\square$

This is why pure periodicity and reversibility of the shift travel together: on
the core — a finite set — $\sigma$ is a bijection, so it decomposes the core
into cycles and every numerator lies on one, rather than on a path running into
one.

## Cycle structure on the core

The core splits into $\{0\}$ and the block $[-q, -1]$, each carried onto itself
by $\sigma$. Only $0$ maps to $0$: $\sigma(p) = 0$ means $p - s_0 q = 0$, so
$p \in \{0, q\}$, and of those only $0$ is in the core. Since $\sigma$ is a
bijection of the core, the rest of the core must go to the rest of the core.

On that block, reduction mod $q$ conjugates $\sigma$ to multiplication by
$2^{-1}$. The reduction map $\pi : [-q, -1] \to \mathbb{Z}/q$ is a bijection
because $[-q, -1]$ is a run of $q$ consecutive integers, hence a complete
residue system — which is also the cleanest reading of the count $q+1$: the $q$
residues, plus the extra point $0$. And $2\sigma(p) = p - s_0 q \equiv p
\pmod q$, so $\pi(\sigma(p)) = 2^{-1}\pi(p)$, with $2$ invertible mod $q$
because $q$ is odd.

The cycle through $p$ therefore has the length of the cycle through $\pi(p)$
under multiplication by $2^{-1}$: the least $k \geq 1$ with
$2^{-k} p \equiv p \pmod q$, equivalently with $q \mid p(2^k - 1)$. Write
$g = \gcd(p, q)$, $p = g p'$, $q = g q'$. Then $q \mid p(2^k-1)$ iff
$q' \mid p'(2^k - 1)$ iff $q' \mid 2^k - 1$, the last step because
$\gcd(p', q') = 1$. The exponents $k$ with $2^k \equiv 1 \pmod{q'}$ are closed
under subtraction, so they are exactly the multiples of the least one, and the
cycle length is

$$\operatorname{ord}_{q/\gcd(p,q)}(2).$$

At $p = -q$ this reads $\operatorname{ord}_1(2) = 1$, recovering that endpoint
as a fixed point. Every cycle length divides $\operatorname{ord}_q(2)$: from
$q' \mid q \mid 2^{\operatorname{ord}_q(2)} - 1$ and the same least-element
argument, $\operatorname{ord}_{q'}(2) \mid \operatorname{ord}_q(2)$. The bound
is attained, at $p = -1$: it lies in the core for every $q \geq 1$ and has
$\gcd(-1, q) = 1$, so its cycle has length $\operatorname{ord}_q(2)$.

When $q$ is prime with $2$ of maximal order, the core is two fixed points plus a
single $(q-1)$-cycle. When $q$ is composite the divisors of $q$ split the core
further: for $q = 9$, the six units mod $9$ form one $6$-cycle
($\operatorname{ord}_9(2) = 6$) while $\{3, 6\}$ form a $2$-cycle
($\operatorname{ord}_3(2) = 2$).

## What this does not say

The core is invariant and $\sigma$ is injective on it, but the core is **not**
saturated: numerators outside it still map into it.

Solving $\sigma(x) = r$ over all of $\mathbb{Z}$ gives $x = 2r + sq$ for
$s \in \{0,1\}$, so the preimages are $2r$ and $2r + q$ — the first even, the
second odd since $q$ is, which is consistent with $s$ being read off as the
parity of $x$. Both are genuine preimages, and they are distinct because
$q \neq 0$. For a core numerator, **exactly one** of the two lies in the core:
at least one because $\sigma$ maps the core onto the core, at most one because
$\sigma$ is injective there. So exactly one preimage lies outside — at $q = 5$,
the core numerator $-3$ has preimages $-6$ and $-1$, of which only $-1$ is in
$[-5, 0]$. That outside preimage is how a sequence with a transient reaches its
periodic part.

The transient folds inward rather than outward. For $p < -q$, both branches of
$\sigma$ are at least $(p-q)/2$, and

$$\frac{p-q}{2} > p \iff p - q > 2p \iff p < -q,$$

so numerators below the core strictly increase toward it and cannot escape
downward.

Above the core the fold runs the other way, and it does not overshoot. For any
$p > 0$, both branches satisfy $\sigma(p) \leq p/2 < p$, so positive numerators
strictly decrease, and being integers they cannot do so forever. Both branches
are also at least $(p - q)/2$: when $p > q$ that is positive, so the numerator
is still positive after the step; when $0 < p \leq q$ it is at least
$(1-q)/2 > -q$. The step that leaves the positive side is therefore taken from
$0 < p \leq q$ with $p$ odd, and it lands at $(p-q)/2 \in [-q/2, 0]$, inside the
core.

## Source

| | |
|---|---|
| `Tools/RegisterSynthesis/fcsrSynthesis.py` | `BM_FCSR`, `FCSR_size` |
| `FeedbackFunctions/FCSR.py` | `state_from_frac` |

## References

`KlapperGoresky1993` — "2-Adic Shift Registers"; recorded in
`bibliography.bib` as the primary reference for the synthesis algorithm
`BM_FCSR` follows.

`KlapperGoresky1997` — extended treatment of the 2-adic theory of FCSRs,
including the equivalence of eventual periodicity and rationality proved above.

`ArnaultBergerMinier2005` — the corrected `D()` helper that `fcsrSynthesis.py`
uses in place of the published version.

`HuSha2019` — bounds on the 2-adic complexity of T-function output sequences;
recorded in `bibliography.bib` as the benchmark `fcsrSynthesis.py` is compared
against. §Two measures of size above is where the units of that comparison are
discussed.
