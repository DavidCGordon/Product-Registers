# FCSR Implementation

How `FeedbackFunctions/FCSR.py` realizes a 2-adic rational as a shift register:
the cell layout, why the carries are stored explicitly, how the connection
integer becomes tap positions, and which bit pattern corresponds to which
fraction.

For the arithmetic this implements — the shift on numerators, which fractions
give purely periodic sequences, and the cycle structure — see
[2-adic Integers and Rational Sequences](../theory/2-adic%20Integers%20and%20Rational%20Sequences.md).
That document is assumed here; this one describes only the realization.

## Why the carries are cells

An LFSR's feedback is a XOR of taps, which is addition in $\mathbb{F}_2$ and
discards carries. An FCSR adds *with* carry, so a carry produced at one clock
must survive to the next. It is state, and it has to live somewhere.

There are two places to put it, mirroring the LFSR situation.

- **Fibonacci form** sums all tapped bits at once: at each clock it forms
  $\Sigma = m + \sum_{\text{taps}} v_i$, emits $\Sigma \bmod 2$, and retains
  $m' = \lfloor \Sigma/2 \rfloor$ as the memory for the next clock. With $k$
  taps, $m \leq k-1$ is invariant — $m' \leq \lfloor (k+m)/2 \rfloor \leq
  \lfloor (2k-1)/2 \rfloor = k-1$ — so the memory ranges over $0, \ldots, k-1$
  and holding it takes a carry *register* of $\lceil \log_2 k \rceil$ bits,
  feeding back into itself. One multi-bit adder, whose size grows with the
  number of taps.
- **Galois form** distributes the addition. Each tapped position gets its own
  one-bit full adder and its own one-bit carry cell, with the output bit
  broadcast to all of them. No adder is wider than one bit, and each carry path
  is local to its position.

This implementation is the Galois form. The interleaved layout is what that
buys: a carry is stored next to the position that produced it.

## Layout

A register of 2-adic complexity $d$ occupies $2d - 1$ cells, interleaved:

```
index:   0       1       2       3       4      ...   2d-2
        val_0   car_0   val_1   car_1   val_2         val_{d-1}
```

**Even indices hold value bits, odd indices hold carry bits.** There are $d$
values and $d-1$ carries — the top value cell has no carry, because nothing is
added past it. `size` is $2d-1$, and the `values` and `carries` properties
expose the paired positions.

Bit $0$ is the output, as elsewhere in the library, and the register shifts
toward index $0$.

**2-adic complexity** here means the value-cell count $d$ — the
`diadic_complexity` parameter of `FCSR.__init__`. The size of the fraction,
$\Phi(p,q) = \max(|p|,|q|)$, is a different quantity and is written $\Phi$
throughout. `FCSR_size` computes a $d$ from a fraction, but it returns a cell
count *sufficient* for `state_from_frac`'s encoding rather than the least one
possible. At the single input $p = q = 1$ the logarithmic sizing would give
$1$ where $2$ is required, and both `FCSR_size` and `state_from_frac`
special-case it to $2$ — §The two special cases below works through why. See
[2-adic Integers and Rational Sequences](../theory/2-adic%20Integers%20and%20Rational%20Sequences.md)
§Two measures of size, which shows that reducing a fraction minimizes $\Phi$
and `FCSR_size` together, while across *different* sequences the two can order
fractions oppositely.

## Which positions are tapped

The taps are the binary digits of $(q+1)/2$, low bit first: position $i$ is
tapped when bit $i$ of $(q+1)/2$ is set. That $q$ is odd is what makes
$(q+1)/2$ an integer.

This is reading tap positions off a coefficient list, with the list stored as
the digits of an integer rather than as coefficients. Since $q$ is odd, $q+1$ is even and its binary expansion has no
constant term, so it can be written $q + 1 = \sum_{i \geq 1} q_i 2^i$ — which
both defines a coefficient list $(q_1, q_2, \ldots)$ and presents the connection
integer as $q = -1 + \sum_{i \geq 1} q_i 2^i$. Bit $i$ of $(q+1)/2$ is bit $i+1$
of $q+1$, so cell $i$ is tapped exactly when $q_{i+1} = 1$. That is the same
index offset a Galois LFSR uses: with constructor polynomial $P$, "bit $i$
receives the departing bit's value iff $P_{i+1} \neq 0$"
([Polynomial Conventions](../conventions/Polynomial%20Conventions.md), §Galois
LFSR). The correspondence is between the coefficient list $(P_i)$ of the LFSR's
polynomial and the base-2 digits $(q_i)$ of $q+1$, with the polynomial's
constant term corresponding to the $-1$ in $q = -1 + \sum_{i \geq 1} q_i 2^i$;
halving $q+1$ is how the implementation strips that offset and recovers the
list.

The top position, $i = d-1$, is special: there is no cell above it to add, so if
it is tapped its value cell simply receives the output bit. That is the feedback
closing the loop. Whether it is tapped at all is the subject of
§Closed and open loop below.

## What each cell computes

Position $i$ takes the value above it and its own carry. Untapped, that is a
**half adder**; tapped, the output bit joins as a third input and it becomes a
**full adder**. Writing $\text{out}$ for the output bit:

| | value cell $2i$ | carry cell $2i+1$ |
|---|---|---|
| untapped | $v_{i+1} \oplus c_i$ | $v_{i+1} \wedge c_i$ |
| tapped | $v_{i+1} \oplus c_i \oplus \text{out}$ | $\mathrm{maj}(v_{i+1},\, c_i,\, \text{out})$ |

The carry cell is exactly the carry-out of the adder whose sum is the value cell
beside it — `AND` for the half adder, majority for the full adder. That is the
whole mechanism.

`FCSR(3, 7)` has $(q+1)/2 = 4 = 100_2$, so only the top position is tapped:

```
[0] value: XOR(VAR(2), VAR(1))          [1] carry: AND(VAR(2), VAR(1))
[2] value: XOR(VAR(4), VAR(3))          [3] carry: AND(VAR(4), VAR(3))
[4] value: VAR(0)
```

`FCSR(3, 5)` has $(q+1)/2 = 3 = 11_2$, tapping positions 0 and 1, so both are
full adders:

```
[0] value: XOR(VAR(2), VAR(1), VAR(0))
[1] carry: XOR(AND(VAR(2),VAR(1)), AND(VAR(2),VAR(0)), AND(VAR(1),VAR(0)))
[2] value: XOR(VAR(4), VAR(3), VAR(0))
[3] carry: XOR(AND(VAR(4),VAR(3)), AND(VAR(4),VAR(0)), AND(VAR(3),VAR(0)))
```

## Closed and open loop

Whether the feedback reaches the top cell is decided by how $q$ sits against
$d$, and the first thing that decides is whether the register exists at all. A
tap at index $i$ needs a value cell at index $i$, and there are only $d$ of
them, indices $0$ through $d-1$. So every set bit of $(q+1)/2$ must sit below
index $d$, i.e. $(q+1)/2 < 2^d$, i.e. $q + 1 < 2^{d+1}$; since $q+1$ is even
this is

$$q \leq 2^{d+1} - 3.$$

Past that bound `FCSR(d, q)` raises `IndexError` — the constructor writes
`fn_list[2*i]` for every tap index $i$, and the list has only $2d-1$ entries.
The admissible denominators are therefore the odd $q$ with
$1 \leq q \leq 2^{d+1} - 3$, and this bound is what the rest of the document
assumes.

Inside that range no bit of $(q+1)/2$ sits above index $d-1$, so the inequality
$(q+1)/2 \geq 2^{d-1}$ can only be carried by bit $d-1$ itself. The top position
is therefore tapped if and only if $(q+1)/2 \geq 2^{d-1}$, i.e. if and only if

$$q \geq 2^d - 1.$$

The upper bound is what makes this an equivalence rather than a one-way
implication: with a set bit above index $d-1$ the inequality would hold while
bit $d-1$ stayed clear.

- **Closed loop, $q \geq 2^d - 1$.** The top value cell receives the output bit,
  and the register is the shift-with-carry described above.
- **Open loop, $q < 2^d - 1$.** The tap indices stop short of the top, and that
  cell holds its own value — a self-loop. `FCSR(3, 5)` above is such a case: its
  top cell reads `VAR(4)`, itself.

These are two regimes, not a correct case and a broken one. They read their
states differently, which is the subject of the next section.

## The state as a numerator

Each (value, carry) pair contributes $\text{value} + 2\cdot\text{carry}$: by
§What each cell computes the carry cell holds the carry-out of the adder whose
sum bit is the value cell beside it, and a carry-out is worth twice its sum bit.
The integer identity in step 1 below is that statement written arithmetically,
and it is what settles the weighting; this sentence only says where the factor
$2$ comes from. Write $v_i$ for the value bits and $c_i$ for the carry bits, and
collect them into the two words

$$a = \sum_{i=0}^{d-1} v_i 2^i, \qquad c = \sum_{i=0}^{d-2} c_i 2^i$$

— the ranges differ because the top value cell has no carry beside it, and a
subscripted $c_i$ is a single bit while a bare $c$ is the whole word. A state is
then read as the integer

$$r = a + 2c,$$

and the register emits the 2-adic expansion of $p/q$ with

$$\boxed{\;p = -r.\;}$$

The numerator is the **negation** of the reading — not a rescaling of it, and
not the reading itself. The formula $r = a + 2c$ is the closed-loop reading; the
open loop flips one weight inside it, and the derivation below is what forces
the flip rather than a separate convention imposed on top.

### One clock as an integer identity

Both claims — that the numerator is $-r$, and that the open loop weights its top
value cell negatively — come out of a single computation: track what one clock
does to $r$ as an integer, and compare it with the numerator map
$\sigma(p) = (p - s_0 q)/2$ of
[2-adic Integers and Rational Sequences](../theory/2-adic%20Integers%20and%20Rational%20Sequences.md).

Write $T = (q+1)/2 = \sum_{i=0}^{d-1} t_i 2^i$ for the tap word (no set bit sits
above index $d-1$, by the bound $q \leq 2^{d+1}-3$), $s = \text{out} = v_0$ for
the output bit, and $x_i = t_i s$ for what the feedback injects at position $i$.
Primes mark the next state.

1. For $i \leq d-2$ the cell pair is a full adder on $v_{i+1}$, $c_i$ and
   $x_i$, so **as integers**
   $$v_{i+1} + c_i + x_i = v_i' + 2c_i'.$$
   This is the $\oplus$ / $\mathrm{maj}$ table of §What each cell computes read
   arithmetically: the tapped row is the full adder on the three inputs, and the
   untapped row is the same identity at $x_i = 0$, where XOR is the sum bit and
   AND is the carry.

2. Weight each identity by $2^i$ and sum over $i = 0, \ldots, d-2$:
   $$\sum_{i=0}^{d-2} 2^i (v_i' + 2c_i') = \frac{a - v_0}{2} + c + s\,(T - t_{d-1}2^{d-1}),$$
   using $\sum_{i=0}^{d-2} 2^i v_{i+1} = \sum_{j=1}^{d-1} 2^{j-1} v_j = (a - v_0)/2$.

3. **Closed loop** ($t_{d-1} = 1$, so $v_{d-1}' = s$). The new reading is
   $r' = \sum_{i=0}^{d-2} 2^i(v_i' + 2c_i') + v_{d-1}' 2^{d-1}$, so with
   $r = a + 2c$ and $(a - v_0)/2 + c = (r - s)/2$,
   $$r' = \frac{a - v_0}{2} + c + s(T - 2^{d-1}) + s\,2^{d-1}
        = \frac{r - s}{2} + s\,\frac{q+1}{2} = \frac{r + sq}{2}.$$

4. **Open loop** ($t_{d-1} = 0$, so $v_{d-1}' = v_{d-1}$), reading the state
   with the top weight flipped from $+2^{d-1}$ to $-2^{d-1}$, i.e.
   $r = a + 2c - 2^{d} v_{d-1}$. Now
   $(a - v_0)/2 + c = (r - s)/2 + 2^{d-1} v_{d-1}$, and
   $$r' = \left[\frac{a - v_0}{2} + c + sT + v_{d-1}2^{d-1}\right] - 2^{d} v_{d-1}
        = \frac{r - s}{2} + sT = \frac{r + sq}{2},$$
   the three $v_{d-1}$ terms cancelling as $2^{d-1} + 2^{d-1} - 2^d = 0$.

5. Both regimes satisfy $r' = (r + sq)/2$, and in both $s = v_0 = r \bmod 2$,
   since $2c$ and $2^d v_{d-1}$ are even. Setting $p = -r$ turns the first into
   $p' = (p - sq)/2$, which is $\sigma$, and the second into $s = p \bmod 2$
   (negation preserves parity), which is $\sigma$'s bit rule.

6. So one clock carries $-r$ exactly as $\sigma$ carries a numerator, and the
   bit the register puts out is the bit $\sigma$ emits. By induction the output
   sequence is the expansion of $p/q$ with $p$ the initial $-r$. A 2-adic
   integer is determined by its expansion and $q$ is fixed, so the numerator is
   determined too: $p = -r$. $\square$

Step 4 is where the negative weight is forced rather than chosen. The unsigned
reading does not survive the open loop: carrying $\hat{r} = a + 2c$ through
step 2 and the same assembly gives
$\hat{r}' = (\hat{r} + sq)/2 + 2^{d-1} v_{d-1}$, which is
$\sigma$'s recurrence only when $v_{d-1} = 0$. And the correction is unique
among weight changes: reading $r = a + 2c + \lambda v_{d-1}$ for a constant
$\lambda$ gives $r' = (r-s)/2 - \lambda v_{d-1}/2 + sT + 2^{d-1}v_{d-1} +
\lambda v_{d-1}$, so $r' = (r + sq)/2$ demands
$\lambda/2 + 2^{d-1} = 0$, that is $\lambda = -2^{d}$.

### The top cell is a sign bit when the loop is open

The two regimes differ in exactly one weight: the top value cell $v_{d-1}$.

| | weight of $v_{d-1}$ in $r$ | range of $p$ |
|---|---|---|
| closed loop | $+2^{d-1}$ | non-positive |
| open loop | $-2^{d-1}$ | either sign |

Step 4 above is where that weight comes from: with the top cell holding its own
value, the clock reproduces $\sigma$ only under the flipped weight, and no other
constant weight works. The closed loop has no such term because its top cell is
overwritten by the output bit each clock instead of held.

There is a picture that agrees with the computation, and it is worth recording
with the step that makes it legitimate. Extend the register upward with value
cells $v_j = 1$ and carry cells $c_j = 0$ for every $j \geq d-1$, all untapped —
admissible, because an open loop has no taps at or above index $d-1$. The
extension starts at the register's own top cell, so this is the picture for a
state with $v_{d-1} = 1$; when that cell is $0$ there is no tail and the weight
contributes nothing either way. Each
extension position is then a half adder with
$v_j' = v_{j+1} \oplus c_j = 1 \oplus 0 = 1$ and
$c_j' = v_{j+1} \wedge c_j = 0$, so the extension is **invariant**, and it hands
a $1$ down to position $d-2$ at every clock — which is exactly what the
self-loop does when the top cell is set. The extension's value is

$$\sum_{j \geq d-1} 2^j = 2^{d-1}(1 + 2 + 4 + \cdots) = 2^{d-1}\cdot(-1) = -2^{d-1}$$

in $\mathbb{Z}_2$, using $\ldots 1111_2 = -1$ — the same weight step 4 forces.
This is the resemblance to two's complement, and that is all it is: a
fixed-width register standing in for an infinite tail of $1$s above it. The
invariance of the extension is what licenses the substitution; the integer
identity is what proves the weight.

This is why the open loop is where positive numerators live: $p = -r$, so a
negative reading is a positive numerator, and only the sign bit can make $r$
negative.

For `FCSR(2, 1)`, which is open ($q = 1 < 3$), the whole state space is eight
states — three cells, $\texttt{[}v_0, c_0, v_1\texttt{]}$ — and at $d = 2$ the
reading $r = a + 2c - 2^dv_{d-1}$ is $r = v_0 + 2c_0 - 2v_1$:

| state | $r$ | $p = -r$ | emits |
|---|---|---|---|
| `[0,0,0]` | $0$ | $0$ | $0,0,0,\ldots$ |
| `[1,0,0]` | $1$ | $-1$ | $1,1,1,\ldots$ |
| `[0,1,0]` | $2$ | $-2$ | $0,1,1,1,\ldots$ |
| `[1,1,0]` | $3$ | $-3$ | $1,0,1,1,1,\ldots$ |
| `[0,0,1]` | $-2$ | $+2$ | $0,1,0,0,\ldots$ |
| `[1,0,1]` | $-1$ | $+1$ | $1,0,0,0,\ldots$ |
| `[0,1,1]` | $0$ | $0$ | $0,0,0,\ldots$ |
| `[1,1,1]` | $1$ | $-1$ | $1,1,1,\ldots$ |

Two things are visible in the eight rows that the prose elsewhere states
abstractly. The map state $\mapsto r$ is many-to-one — `[0,0,0]` and `[0,1,1]`
both read $0$, `[1,0,0]` and `[1,1,1]` both read $1$ — which §The representable
range below turns into the reason a core numerator does not imply a cyclic
state. And the six attained numerators are exactly the integers of $[-3, +2]$,
each once or twice, which is the open-loop range that section computes at
$d = 2$.

By contrast `FCSR(3, 7)`, closed since $q = 7 = 2^3 - 1$, reads every cell
positively: state `[1,0,0,0,0]` has $r = 1$, so $p = -1$ and the output is the
expansion of $-1/7$, namely $1,0,0,1,0,0,\ldots$; state `[1,0,1,1,0]` has
$a = 3$, $c = 2$, so $r = 7$, $p = -7$, and the output is the expansion of
$-7/7 = -1$, all ones.

### The representable range

In the closed loop every weight is positive, so the bound follows from the
weighting by setting every cell:

$$\max r = (2^d - 1) + 2(2^{d-1} - 1) = 2^{d+1} - 3,$$

and the representable numerators are exactly

$$p \in [-(2^{d+1}-3),\; 0],$$

which is $2^{d+1}-2$ values, and does not depend on $q$. Every intermediate
value is attained: the value bits alone cover $[0, 2^d-1]$ in steps of one, and
the carries extend the run without gaps.

Against the periodic core $[-q, 0]$, of size $q+1$: the representable range is
wider whenever $q < 2^{d+1}-3$, and never narrower, since $q \leq 2^{d+1}-3$ is
the constructibility bound of §Closed and open loop. The numerators lying
between the two — those in $[-(2^{d+1}-3),\, -q-1]$ — are exactly the
representable numerators whose sequences are not purely periodic. That is the
pure-periodicity proposition of
[2-adic Integers and Rational Sequences](../theory/2-adic%20Integers%20and%20Rational%20Sequences.md)
(purely periodic iff $-q \leq p \leq 0$) restricted to the representable range.

It is a statement about numerators, and it transfers to *states* in one
direction only. A state whose numerator lies outside the core emits a sequence
with a transient, so that state cannot be on a cycle of the register: a state on
a cycle repeats, and so does everything read off it. The converse fails, because
state $\mapsto r$ is many-to-one. In `FCSR(2, 3)` — closed, since
$q = 3 = 2^2 - 1$ — the states `[0,1,0]` and `[0,0,1]` both read $r = 2$, the
first as $2c$ with $c = 1$, the second as $a = 2$ with $c = 0$. Both therefore
emit the expansion of $-2/3$, which is purely periodic. But both clock to
`[1,0,0]`: position $0$ is untapped, so $v_0' = v_1 \oplus c_0 = 1$ and
$c_0' = v_1 \wedge c_0 = 0$ in each case, and the closed top cell takes
$v_1' = v_0 = 0$. The clock map restricted to the states lying on cycles is a
bijection of that set — each cyclic state has its successor and its predecessor
inside its own cycle — so two distinct states with the same successor cannot
both be cyclic. At least one of `[0,1,0]`, `[0,0,1]` is therefore on a
transient while carrying a core numerator.

In the open loop the same count of states covers a range shifted by the sign
bit, running from $-(2^{d-1}-1) - 2\cdot(2^{d-1}-1)$ up to $+2^{d-1}$ — the
positive end being what a terminating expansion needs.

## Encoding a fraction: `state_from_frac`

Going the other way — given $p/q$, produce a state — inverts the weighting. The
two branches of `state_from_frac` split on the sign of $p$, which is not the
same division as the two regimes: the non-positive branch is regime-independent,
producing a state that reads identically whether the loop is open or closed,
while the positive branch is pinned to the open loop.

**Non-positive $p$.** The magnitude $|p|$ is split as $a + 2c$, distributed by
thirds: $\lfloor |p|/3 \rfloor$ to each of $a$ and $c$, with the remainder going
to $a$ (remainder 1) or $c$ (remainder 2). The size is

$$d = \max\big(\lceil \log_2(|p|/3 + 1)\rceil + 1,\ \lceil \log_2 q \rceil\big).$$

The second term is what makes the register exist at all. §Closed and open loop
requires $q \leq 2^{d+1} - 3$, and $d \geq \lceil \log_2 q \rceil$ gives
$q \leq 2^d$, which clears that bound: for $d \geq 2$, $2^d \leq 2^{d+1} - 3$
since $2^d \geq 3$; and $d = 1$ forces $q \leq 2$, hence $q = 1 = 2^{d+1} - 3$
for odd $q$. The first term is what makes the split fit, which is the next
paragraph.

The split is a *choice*, and what it has to satisfy is $a < 2^{d-1}$ and
$c < 2^{d-1}$: the value word must fit the $d$ value cells with the top one
clear, and the carry word must fit the $d-1$ carry cells. Not every split with
$a + 2c = |p|$ satisfies it. At $|p| = 3$, $q = 1$, $d = 2$ the thirds split
$a = c = 1$ gives `[1,1,0]`, which emits the expansion of $-3$; the alternative
$a = 3$, $c = 0$ gives `[1,0,1]`, which sets the top value bit — the sign bit of
§The top cell is a sign bit when the loop is open — and emits $1,0,0,0,\ldots$,
the expansion of $+1$ instead.

The thirds split does satisfy it, by residue. The sizing gives
$2^{d-1} \geq |p|/3 + 1$, and $2^{d-1}$ is an integer:

- $|p| = 3k$: $a = c = k$, and $2^{d-1} \geq k + 1 > k$.
- $|p| = 3k+1$: $a = k+1$, $c = k$, and $2^{d-1} \geq k + 4/3$; since
  $2^{d-1}$ is an integer and $k + 4/3$ is not, $2^{d-1} \geq k + 2 > k+1$.
- $|p| = 3k+2$: $a = k$, $c = k+1$, and $2^{d-1} \geq k + 5/3$, again a
  non-integer bound, so $2^{d-1} \geq k + 2 > k+1$.

Keeping both components small is a secondary benefit — it is what keeps $d$
itself small — but it is $a, c < 2^{d-1}$ that makes the state correct.

That bound is also what makes the branch regime-independent. With
$a < 2^{d-1}$ the top value bit $v_{d-1}$ is $0$, and at $v_{d-1} = 0$ the open
reading $a + 2c - 2^dv_{d-1}$ and the closed reading $a + 2c$ agree; either way
$r = a + 2c = |p|$, so $p = -r$ as required. The branch needs that
independence, because it lands in both regimes — the two fractions of the next
paragraph are one
of each: $-1/3$ sits on `FCSR(2, 3)`, closed at $q = 3 = 2^2 - 1$, and $-3/9$ on
`FCSR(4, 9)`, open at $q = 9 < 15$.

The fraction is deliberately **not** reduced first, and the reason is that
reducing it would silently change which register you are talking about.

A fraction and its reduction generate the *same sequence* — $-3/9$ and $-1/3$
both expand to $1,0,1,0,\ldots$ — but they name different registers. The
denominator is the connection integer, so reducing changes the tap word
$(q+1)/2$ and therefore the whole machine: $-3/9$ lives on `FCSR(4, 9)` in state
`[1,1,0,0,0,0,0]`, while $-1/3$ lives on `FCSR(2, 3)` in state `[1,0,0]`. Larger
$q$, larger numerator, larger register, same output.

Keeping the fraction as given is what preserves that distinction. A caller who
asks for $p/q$ gets a state for the register on $q$, not for whichever smaller
register happens to emit the same bits — which matters because `fromSeq` passes
the recovered denominator straight to the constructor, so the register and the
state have to agree on which $q$ they are built for.

**Positive $p$.** Here the encoding sets $a = 2^{d} - p$ with $c = 0$, and sizes

$$d = 1 + \lceil \log_2 \max(p, q) \rceil,$$

so that $2^{d} \geq 2\max(p,q)$ and hence $q < 2^{d} - 1$ — the open loop, which
is where a terminating expansion is possible.

The same sizing gives $2^{d-1} \geq p$, hence $2^{d-1} \leq a = 2^d - p < 2^d$:
the top value bit is set, and the open-loop reading returns
$r = a - 2^d = -p$, so the numerator is $-r = +p$. What is stored is the two's
complement pattern of $-p$ at width $d$, not the bits of $p$: for $p = 3$,
$q = 1$ the state is $a = 5 = 101_2$, i.e. `[1,0,0,0,1]`, whose value bits
low-first are $1,0,1$ while $p = 3$ is $11_2$.

What comes out is the expansion of $p/q$ — not of $p$. A terminating expansion
is a finite sum $\sum_{i < N} s_i 2^i$, and the values of those sums are exactly
the non-negative integers, so the output terminates if and only if $p/q$ is a
non-negative integer. Within this branch, where $p > 0$, that is the condition
$q \mid p$, and in particular it holds whenever $q = 1$. The hypothesis $p > 0$
is doing work: $-3/3 = -1$ also has $q \mid p$, and its expansion
$1,1,1,\ldots$ does not terminate. For $p = 3$,
$q = 1$ the output is $1,1,0,0,\ldots$, the binary digits of $3$ low bit first
— reconstructed by the carry chain, not the stored pattern shifting out, as the
stored value bits are $1,0,1$. For $p = 1$, $q = 3$ the state is `[1,0,1,0,1]`
and the output is $1,1,0,1,0,1,\ldots$, which never terminates: under $\sigma$
the numerators run $1 \mapsto -1 \mapsto -2 \mapsto -1$, and the last two
repeat.

**The two special cases.** Both are $q = 1$, and both are there because the
sizing above is tight at the very bottom of its range.

- $0/1$ returns a single cell holding $0$: the one fraction needing no register
  at all, emitting zeros forever.
- $1/1$ returns two value cells in state `[1,0,1]`. The sizing would give
  $d = 1$, where $q = 2^{d} - 1$ exactly — closed, and so emitting the expansion
  of $-1/1$, all ones, rather than of $1/1$. Two cells open the loop, and
  `[1,0,1]` is the state whose sign bit supplies the $-2$ making $r = -1$ and
  therefore $p = +1$. It emits the leading $1$ and then falls into `[0,1,1]`,
  which is absorbing and emits zeros.

$p = q = 1$ is the only input where the logarithmic sizing fails to clear
$q < 2^d - 1$; every other positive numerator has $\max(p,q) \geq 2$ and clears
it with room.

## API

```python
FCSR(diadic_complexity, q)          # d value cells, d-1 carry cells, width 2d-1
FCSR.fromSeq(seq)                   # -> (state, register)
FCSR.fromReg(F, bit, numIters)      # -> (state, register), observing one bit
FCSR.state_from_frac(num, den)      # -> (d, state)
fcsr.carries                        # the d-1 carry functions
fcsr.values                         # the paired value functions
```

`state_from_frac` takes the fraction as it is: `num` carries its own sign and
`den` is the positive connection integer. `fromSeq` runs `BM_FCSR` to recover
$p/q$ from the observation, builds the state from it, then constructs the
register from $q$.

## Source

| | |
|---|---|
| `FeedbackFunctions/FCSR.py` | `__init__`, `state_from_frac`, `fromSeq`, `fromReg`, `carries`, `values` |
| `Tools/RegisterSynthesis/fcsrSynthesis.py` | `BM_FCSR`, `FCSR_size` |
