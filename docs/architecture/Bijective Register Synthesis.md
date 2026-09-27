# Bijective Register Synthesis

What the `bijective=True` option means in `Tools/RegisterSynthesis`, why the
condition it imposes is the right one, and how the same condition is applied on
the linear and nonlinear paths.

For Berlekamp-Massey itself — which polynomial it returns, and in which
coefficient convention — see [Polynomial Conventions](../conventions/Polynomial%20Conventions.md).

## Two maps, kept apart

A Fibonacci-configured register of length $m$ carries two functions, and the
index arithmetic below is unreadable unless they are held apart:

- the **feedback function** $f : \mathbb{F}_2^m \to \mathbb{F}_2$, which
  produces one bit;
- the **update map** $F : \mathbb{F}_2^m \to \mathbb{F}_2^m$, which produces the
  next state, $F(s) = (s_1, \ldots, s_{m-1}, f(s))$.

Synthesis searches for $f$. Invertibility is a property of $F$. The whole point
of this option is that the second can be imposed as a condition on the first.

Seed the register with the first $m$ bits of a sequence and the state at time
$t$ is the length-$m$ **frame** $s_t = (\text{seq}[t], \ldots, \text{seq}[t+m-1])$.
Under that identification $F(s_t) = s_{t+1}$ holds automatically in the shifted
coordinates — both sides carry $\text{seq}[t+1..t+m-1]$ — and only the last
coordinate is a real constraint:

$$\text{the register generates seq} \iff f(\text{seq}[t], \ldots, \text{seq}[t+m-1]) = \text{seq}[t+m] \quad\text{for every } t.$$

## When is the update map a bijection?

$F$ discards $s_0$ and shifts the rest, so $s_1, \ldots, s_{m-1}$ survive into
$F(s)$ unchanged. Inverting $F$ is therefore free in every coordinate but one:
the question is whether $s_0$ can be recovered from $F(s)$.

Since $x^2 = x$ over $\mathbb{F}_2$, every Boolean function has a unique
multilinear representative — its ANF (see
[Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md)) — so $f$ splits
uniquely on whether a monomial contains $s_0$:

$$f(s) = s_0 \cdot A(s_1, \ldots, s_{m-1}) \;\oplus\; B(s_1, \ldots, s_{m-1}).$$

> **Proposition.** $F$ is a bijection if and only if $A \equiv 1$.
>
> *Proof.* Fix $u = (s_1, \ldots, s_{m-1})$. The two states $(0, u)$ and $(1, u)$
> have the same image under the shift, so they collide under $F$ exactly when
> $f(0,u) = f(1,u)$, i.e. when $B(u) = A(u) \oplus B(u)$, i.e. when $A(u) = 0$.
>
> If $A \equiv 1$ no such collision exists, and distinct states differing
> anywhere else already differ in the shifted coordinates, so $F$ is injective —
> hence bijective, the domain being finite. If $A(u) = 0$ for some $u$, then
> $(0,u)$ and $(1,u)$ collide and $F$ is not injective. $\square$

So bijectivity is exactly the requirement that the discarded bit enter the
feedback **linearly and alone**:

$$f(s) = s_0 \oplus B(s_1, \ldots, s_{m-1}).$$

A linear register satisfies this whenever bit $0$ is a tap, since all its
monomials are singletons and $A$ is then the constant $1$.

## The reduction: search over the internal window

Substituting the constrained shape into the generation condition gives

$$\text{seq}[t+m] = \text{seq}[t] \;\oplus\; B(\text{seq}[t+1], \ldots, \text{seq}[t+m-1]).$$

Both occurrences of $\text{seq}$ outside $B$ are observed values, so move them
together:

$$B(\underbrace{\text{seq}[t+1], \ldots, \text{seq}[t+m-1]}_{\textstyle w_t,\ m-1 \text{ bits}}) \;=\; \underbrace{\text{seq}[t+m] \oplus \text{seq}[t]}_{\textstyle d_t}.$$

Two things changed, and one did not.

- The argument is the **internal window** $w_t$, of $m-1$ bits: the frame with
  the discarded bit removed, that bit having moved to the other side.
- The target $d_t$ is a lag-$m$ difference rather than a single observed bit;
  the pair $(w_t, d_t)$ spans $\text{seq}[t..t+m]$, one position wider than a
  state.
- $B$ is **unconstrained**. The bijectivity requirement is now carried entirely
  by the *shape* of the equation, not by any restriction on $B$.

That is the whole change. The search no longer ranges over feedback functions on
the full frame subject to an invertibility side-condition; it ranges over $B$ on
the internal window, with no side-condition at all. Every $B$ yields a bijective
register, and every bijective register arises this way.

## Imposing the existing condition on that set

What remains is a consistency question, and it is the same question each path
already answered — asked about $(w_t, d_t)$ instead of about frames and next
bits. A fit exists at length $m$ exactly when the pairs are consistent; the two
paths differ only in what $B$ is allowed to be.

**Nonlinear** (`BM_NL`). $B$ may be any Boolean function, so the pairs are
consistent iff no two agree on the window and disagree on the target:

$$w_{t_1} = w_{t_2} \implies d_{t_1} = d_{t_2}.$$

This is a lookup table. When it is consistent, $B$ is read off it directly as
the XOR of the minterms selecting the windows where $d_t = 1$, and $f$ is
assembled as $s_0 \oplus B$. Note this is the same obstruction the unconstrained
`BM_NL` grows its frame to avoid — a repeated argument demanding two different
outputs — applied to the internal window and the lag-$m$ difference.

**Linear** (`berlekamp_massey`). $B$ must itself be linear, so
$B(w_t) = d_t$ is not a table but a system of $\mathbb{F}_2$-linear equations in
the $m-1$ unknown tap coefficients. Consistency is solvability, which
elimination decides; a free variable is a genuine choice among fits.

In both cases the search takes the least $m$ for which the condition holds, so
the result is the shortest bijective register generating the sequence. Because
the bijective feedbacks are a subset of those the unconstrained search ranges
over, that length is at least the unconstrained one, and typically equal or one
greater.

## The precondition: the sequence must be purely periodic

A bijection on a finite set puts every element on a cycle, so a bijective
register emits a **purely periodic** sequence from any seed — one with no
transient. This is a condition on the input, not merely a property of the
output, and it is the one way to use the flag wrongly.

If the sequence has a genuine transient, no bijective register reproduces it,
and the search does not fail cleanly — it degenerates. Write $\rho$ for the least
period and $\tau$ for the minimal preperiod. Positions $t_1 = \tau - 1$ and
$t_2 = \tau - 1 + \rho$ have equal windows, since every index in them is at least
$\tau$ and so lies in the periodic part; and their targets differ, because
$d_{t_1} \oplus d_{t_2} = \text{seq}[\tau-1] \oplus \text{seq}[\tau-1+\rho] = 1$
by minimality of $\tau$. Both positions are in range whenever
$m \leq N - \tau - \rho$, so every such length is inconsistent and the returned
length is at least $N - \tau - \rho + 1$ — growing with the length of the
observation rather than reporting a property of the source.

The symptom is therefore a returned length that tracks $N$. Fitting two prefixes
of different lengths distinguishes the cases: a genuine fit is stable, a
degenerate one grows with slope $1$.

## API

`bijective` is accepted uniformly across the batch fits, their iterators, and
both register constructors. It defaults to `False` everywhere, and `False` is
the original behaviour in every case.

```python
berlekamp_massey(seq, bijective=False)                          # -> (L, polynomial)
BM_NL(seq, bijective=False)                                     # -> (m, f)

berlekamp_massey_iterator(seq, yield_rate, bijective=False)     # -> yields (L, polynomial)
BM_NL_iterator(seq, yield_rate, yield_corrected, bijective=False)  # -> yields (m, f)

Fibonacci.fromSeq(seq, nonlinear=False, bijective=False)        # -> (state, register)
Galois.fromSeq(seq, bijective=False)                            # -> (state, register)
```

Each iterator agrees with its batch counterpart on the full sequence. With
`bijective=True` the fit is a rescan over $m$ rather than an incremental update,
so the iterators have no streaming core to reuse there: each yield reports the
batch fit for the prefix consumed so far, at the batch cost.

`Galois.fromSeq` delegates to the same constrained search. Bijectivity is a
property of the polynomial, not of the realization — the Galois and Fibonacci
registers built from one polynomial have update matrices that are transposes of
each other, so one is invertible exactly when the other is.

A register returned under `bijective=True` reproduces the sequence, passes
`invert()`, and has preperiod $0$ from every seed. When no bijective register
shorter than the sequence generates it, the search raises `ValueError` rather
than returning a register that does not.

## Source

| | |
|---|---|
| `Tools/RegisterSynthesis/lfsrSynthesis.py` | `berlekamp_massey`, `berlekamp_massey_iterator` |
| `Tools/RegisterSynthesis/nlfsrSynthesis.py` | `BM_NL`, `BM_NL_iterator`, `KMP_table` |
| `FeedbackFunctions/Fibonacci.py` | `fromSeq`, `invert`, `_inverse_feedback` |
| `FeedbackFunctions/Galois.py` | `fromSeq`, `invert` |
