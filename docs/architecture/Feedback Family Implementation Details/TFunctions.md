# T-Functions

This page describes `TFunction` (`FeedbackFunctions/TFunction.py`), the feedback
function PyPR provides for triangular state-update maps: what a T-function is,
which one the constructor builds, the bit order it uses, and how the class sits
inside the `CMPR` interface. The period and linear-complexity theory of
T-functions is developed in
[Root Multiplicities and Jordan Decomposition](../../theory/Root%20Multiplicities%20and%20Jordan%20Decomposition.md#t-functions-and-binomial-sequences);
this page covers only the implementation.

## Definition

A **T-function**, as introduced by Klimov and Shamir (`KlimovShamir2002`), is a
state-update map $F : \mathbb{F}_2^n \to \mathbb{F}_2^n$ that is triangular in
the fixed order of the bit indices: output bit $i$ depends only on input bits
$0, \ldots, i$. This library uses the term slightly more generally, allowing the
order to be chosen: the bits are placed in a **dependency order** such that the
update of each bit depends only on that bit and on the bits before it in the
order. In the form
this library works with, a bit keeps its own value and adds a function of the
earlier bits, so any nonlinear dependence is on earlier bits only:

$$b_k[t+1] = b_k[t] \oplus f_k(b_0[t], \ldots, b_{k-1}[t]),$$

where $k$ counts positions in the dependency order.

T-functions are a general class. The binary counter described next is one
member of it, and it is the one `TFunction` constructs.

## The default construction

`TFunction(n)` builds an $n$-bit binary counter: each clock adds $1$ modulo
$2^n$. Writing $x_0, \ldots, x_{n-1}$ for the state bits and $x_i'$ for their
values after one clock, the update is

$$x_{n-1}' = x_{n-1} \oplus 1, \qquad x_i' = x_i \oplus \prod_{j=i+1}^{n-1} x_j \quad (0 \leq i < n-1).$$

Bit $n-1$ is the least significant bit. It flips at every clock, and bit $i$
flips exactly when every bit of higher index is $1$, which is the carry rule of
binary addition: a carry reaches a position when all less significant positions
are set. The integer the state represents is therefore

$$\sum_{i=0}^{n-1} x_{n-1-i}\, 2^i,$$

and the state list $[x_0, \ldots, x_{n-1}]$ reads as that integer's binary
numeral with the most significant digit first. For $n = 3$, starting from the
zero state, the register runs

```
[0,0,0] -> [0,0,1] -> [0,1,0] -> [0,1,1] -> [1,0,0] -> ...
```

The step from `[0,1,1]` shows all three rules at once: $x_2' = 1 \oplus 1 = 0$,
$x_1' = 1 \oplus x_2 = 0$, and $x_0' = 0 \oplus x_1 x_2 = 1$.

Because the counter passes through every residue modulo $2^n$ before returning,
the state space is a single cycle of length $2^n$. The tests in
`test/FeedbackFunctions/TFunction/tfunction_test.py` check all three statements
at $n = 4$: the represented integer increases by one per clock, sixteen distinct
states are visited, and the period is $16$.

## Bit order

Placing the least significant bit at index $n-1$ is the reverse of the usual
counter layout, and it is deliberate. In a `CMPR`, chaining propagates from
higher bit indices toward lower ones: block $0$, the source block with no
chaining inputs, occupies the highest bit indices
([Notation and Terminology](../../conventions/Notation%20and%20Terminology.md),
§CMPR Block Ordering and T-Functions). With the least significant bit at the top
index, each bit of the counter depends on the bits of higher index, so the
T-function's dependency order runs in the same direction as CMPR chaining.

The attribute `induction_order` records that order as a list of bit indices,
`[n-1, n-2, ..., 0]`: bit $i$ appears only after every bit $j > i$. It is the
order in which an inductive argument over the triangular structure visits the
bits, starting from the bit with no dependencies.

The theory document indexes bits by position in the dependency order instead,
so its bit $k$ is this class's bit $n-1-k$. Its bit $0$, the one with no
chaining inputs and the constant $\oplus\,1$, is bit $n-1$ here.

## Relation to CMPR

`TFunction` is a subclass of `CMPR`. It is assembled from $n$ single-bit blocks,
each an `MPR` of size $1$ with polynomial $1 + x$, whose update leaves its bit
unchanged; the products and the constant in the update formula above are then
added to those blocks as chaining terms. The consequences for a user are the
following.

- A `TFunction` can be passed to `FeedbackRegister` directly, or used as a
  component of a larger register. The tests build `CMPR([T, M7, M5, M3])` from a
  `TFunction(5)` and three `MPR`s, install chaining with `generate_chaining`,
  run the result, and compute its monomial profiles.
- `root_expressions` is overridden. For the bit at position $k$ of the
  dependency order it returns a single `JordanSet` with no field roots and
  Jordan lengths $\{1, \ldots, 2^k + 1\}$. The maximum $2^k + 1$ is the value
  derived in
  [Root Multiplicities and Jordan Decomposition](../../theory/Root%20Multiplicities%20and%20Jordan%20Decomposition.md#period-doubling-through-chaining).
- `monomial_profiles` is overridden as well, with a construction specific to the
  counter in place of the general `CMPR` computation.
- The mesh optimization does not apply, because every block has the same
  polynomial and size $1$
  ([Mesh Optimization](../Mesh%20Optimization.md), §Why Root Expressions Need
  Pairwise Coprime Sizes).
- `estimate_LC` runs on a T-function but was not designed for one; see
  [Root Expressions and LC Estimation](../../theory/Root%20Expressions%20and%20LC%20Estimation.md#t-functions)
  for its limits there.

## Usage

```python
from PyPR.FeedbackFunctions import TFunction
from PyPR.FeedbackRegister import FeedbackRegister

T = TFunction(4)                    # 4-bit counter, bit 3 least significant
reg = FeedbackRegister(0, T)        # seeded with the zero state
for state in reg.run(compiled=False, limit=16):
    ...                             # visits all 16 states
```

## Source

| | |
|---|---|
| `FeedbackFunctions/TFunction.py` | `TFunction`, `induction_order`, `root_expressions`, `monomial_profiles` |
| `test/FeedbackFunctions/TFunction/tfunction_test.py` | counter behaviour, period, use as a `CMPR` component |

## References

`KlimovShamir2002` — "A New Class of Invertible Mappings", the paper that
introduces T-functions.
