# Logic Optimization

`PyPR.BooleanLogic.Optimization.optimize` rewrites a list of
`BooleanFunction`s as one DAG that computes the same outputs in fewer two-input
operations. It is reached through the `optimize` flag of `generate_python`, on
both `BooleanFunction` and `FeedbackFunction`. The flag is off by default, and
`compile` turns it on, since compiled code evaluates exactly the statements
`generate_python` returns. This document states what the rewrite assumes, how
the network is represented, the rule and source of each pass, what is and is
not guaranteed, how the flag is wired, and the measurements behind those
statements.

```python
fn.generate_python()                # each bit as written (the default)
fn.generate_python(optimize=True)   # the optimized update, where that saves operations

fn.compile()                        # compiles generate_python(optimize=True)
fn.compile(optimize=False)          # compiles each bit as written

smaller = fn.optimize()             # the optimized update as a FeedbackFunction, to inspect

from PyPR.BooleanLogic.Optimization import operation_count, optimize
operation_count(fn.fn_list), operation_count(optimize(fn.fn_list))   # before, after
```

Source: `src/PyPR/BooleanLogic/Optimization.py`,
`FeedbackFunction.optimize` / `generate_python` / `compile` /
`_straight_line_program` in
`src/PyPR/FeedbackFunctions/FeedbackFunction.py`. Tests:
`test/BooleanLogic/test_optimization.py`.

## 1. What is assumed

The functions are read as Boolean: every `VAR` takes the value 0 or 1, and a
`CONST` whose value is the integer 0 or 1 is that bit. Every rewrite below is an
identity of Boolean algebra, so under this reading each output computes the same
function of the bits before and after.

The reading is narrower than what a `BooleanFunction` can be evaluated on.
`eval_ANF` and `remap_constants` evaluate the same DAG over other objects
(ANFs, monomial profiles), where identities such as $x \wedge x = x$ need not
hold. `optimize` is therefore not applied to `fn_list` itself: `compile` uses
the optimized form only for the code it emits, and `FeedbackFunction.optimize`
returns a new object.

Two kinds of node are **atoms**, carried through unchanged over optimized
arguments: a node whose class is not exactly one of `VAR`, `CONST`, `NOT`,
`XOR`, `AND`, `OR`, `XNOR`, `NAND`, `NOR` (a subclass may mean something else),
and a `CONST` whose value is not the integer 0 or 1. An atom is never merged
with another atom, even an identical one.

## 2. The network

The functions are first translated into a network held in `_Network`.

- **Literals.** A reference to a node is a *literal* $2n + c$, where $n$ is the
  node's number and $c \in \{0, 1\}$ says whether the reference is
  complemented. Node 0 is the constant 0, so literal 0 is false and literal 1
  is true. A `NOT` is not a node: it is the complement bit on a reference.
  This is the and-inverter-graph convention of complemented edges
  (`KuehlmannParuthiKrohmGanai2002`).
- **Nodes.** Variables, atoms, and three kinds of gate: XOR, AND and OR, each
  over a set of two or more literals, stored as a sorted tuple.
- **Structural hashing.** A gate is looked up by its kind and its sorted terms
  before it is created, so two gates with the same kind and the same inputs,
  in any order, are one node (`KuehlmannParuthiKrohmGanai2002`). This is what
  shares logic between outputs: the same product written in two bits'
  functions is built once.

Every gate is simplified as it is built, by these identities:

| Kind | Identities applied |
|------|--------------------|
| XOR | $x \oplus \bar y = \overline{x \oplus y}$, so complements move to the gate's output and XOR terms are always uncomplemented; $x \oplus x = 0$; $x \oplus 0 = x$ |
| AND | $x \wedge 0 = 0$; $x \wedge 1 = x$; $x \wedge x = x$; $x \wedge \bar x = 0$ |
| OR | $x \vee 1 = 1$; $x \vee 0 = x$; $x \vee x = x$; $x \vee \bar x = 1$ |
| AND, OR | De Morgan when every input is complemented: $\bar a \wedge \bar b = \overline{a \vee b}$ and $\bar a \vee \bar b = \overline{a \wedge b}$, one complement in place of one per input |

A gate left with no terms is a constant, and with one term is that term.
`XNOR`, `NAND` and `NOR` are the complemented literals of XOR, AND and OR gates.

## 3. The cost model

The cost of a network (`_cost`) is the number of two-input operations needed to
evaluate its outputs with each node computed once:

- a gate with $n$ terms costs $n - 1$;
- a complement costs 1 per node that is referenced complemented anywhere, since
  the complemented value is computed once and shared;
- except that an XOR gate referenced *only* complemented costs nothing extra
  when one of its terms $a$ is itself referenced complemented elsewhere:
  $\overline{a \oplus b} = \bar a \oplus b$, so the gate is built directly over
  $\bar a$ (`_free_complements`).

`operation_count` is the same count taken on `BooleanFunction` objects: each
node object once, $n - 1$ per $n$-argument gate, 1 per `NOT`, and one more for
`XNOR`, `NAND`, `NOR`. It counts the DAG as built, so equal logic held in two
separate objects is counted twice. That is the baseline `optimize` is measured
against.

## 4. The passes

`optimize` builds the network, then runs three passes. A pass's result is kept
only if it lowers the cost, so the cost never rises from one step to the next.

### Absorbing (`_rebuild`)

A gate used exactly once, by a gate of the same kind, is merged into its user:
$(a \oplus b) \oplus c$ becomes $a \oplus b \oplus c$. For AND and OR the
reference must be uncomplemented. A gate with any other user is kept whole,
because that user shares it. Every pass ends with a rebuild, which also
re-simplifies and re-hashes each gate.

### Pair extraction (`_extract_pairs`)

Run separately on the XOR gates, the AND gates and the OR gates. If two terms
$a, b$ occur together in $k \ge 2$ gates of one kind, a gate $p = a \circ b$ is
made once and each of the $k$ gates uses $p$ in place of the pair. Each gate
loses one operation and $p$ costs one, so the saving is $k - 1$; when a gate
that is exactly $\{a, b\}$ already exists it serves as $p$ and the saving is
the number of other gates. The pair with the largest $k$ is taken first, ties
broken by literal number, and the counts are updated as gates change, including
for pairs involving the new term $p$.

For XOR this is Paar's greedy algorithm (`Paar1997`). For AND, a shared pair is
a common two-literal cube, the single-cube divisor of
`RajskiVasudevamurthy1992`; the same procedure is applied to OR by duality.

Counting pairs is quadratic in a gate's size, so the pass counts at most
`PAIR_BUDGET` pairs per kind. Gates are counted smallest first, only among
terms that occur in two or more gates, and a gate that does not fit the
remaining budget takes no part. The budget bounds running time, not
correctness: a skipped gate is simply left as it is.

### Factoring (`_factor`)

Within an XOR gate, if a literal $x$ occurs in $k \ge 2$ of the AND gates that
are its terms, those terms are replaced by one:

$$x a_1 \oplus x a_2 \oplus \dots \oplus x a_k = x \,(a_1 \oplus a_2 \oplus \dots \oplus a_k)$$

where $a_i$ is the product with $x$ removed. The number of XORs is unchanged
($k - 1$ move inside the bracket) and $k$ ANDs by $x$ become one, a saving of
$k - 1$. The literal in the most products is taken first, ties broken by
literal number, and the cofactor sum $a_1 \oplus \dots \oplus a_k$ is factored
in turn. The same identity is applied to an OR of ANDs,
$x a \vee x b = x (a \vee b)$, and an AND of ORs,
$(x \vee a)(x \vee b) = x \vee a b$.

Only an inner gate used once, by the gate being factored, takes part. One with
other users is shared logic, and factoring it out would leave those users
computing it again.

This is algebraic division by a single-literal divisor, the simplest case of
the kernel-based factoring of `BraytonMcMullen1982`;
`HosangadiFallahKastner2006` apply the same methods to polynomials, of which a
sum (XOR) of products (AND) over GF(2) is an instance.

### Order

Extraction, factoring, then extraction again. Extraction runs first so that a
product shared between outputs is shared before factoring could consume it, and
again afterwards because the cofactor sums factoring creates can have pairs in
common.

## 5. What is guaranteed

- **Same outputs.** Each returned function equals its input on every
  assignment of bits (section 1).
- **Never more operations.** `operation_count(optimize(fs)) <=
  operation_count(fs)`. Each pass is accepted only on a lower cost, and
  `optimize` finally compares `operation_count` on the result against the
  input and returns a copy of the input if the result is larger. The final
  check exists because the network's cost model and `operation_count` can
  differ where the input happened to place a complement better.
- **Deterministic.** Every choice is broken by node number, so the same input
  gives the same output, and `compile` the same code.
- **Not minimal.** Finding a smallest program is NP-hard already for XOR gates
  alone: that is the shortest linear program problem
  (`BoyarMatthewsPeralta2013`). Every pass here is a greedy heuristic.

Known limits, each a place a better result exists that these passes do not
find:

- XOR sharing is cancellation-free: a pair is shared only where both terms
  literally occur. `BoyarMatthewsPeralta2013` give a heuristic that allows
  cancellation.
- A gate with several users is never flattened into them, so a cancellation
  that only appears after expanding a shared sum is not found.
- Factoring divides by single literals only, not by multi-term kernels.
- AND and OR gates are matched as written, apart from the all-complemented De
  Morgan case, so the same function written once as an AND and once as an OR of
  complements is not recognised as shared.

## 6. The `optimize` flag: `generate_python` and `compile`

The flag belongs to `generate_python`, and defaults to False there.

- **`BooleanFunction.generate_python(optimize=False)`.** With False, the
  function is written as it stands, each of its subfunctions assigned to a
  variable. With True, `optimize([self])` is written the same way instead,
  provided it has strictly fewer operations. `overrides` cannot be combined
  with True and raises `ValueError`: an override names a node of the original
  function, and the optimized form is made of different nodes.
- **`FeedbackFunction.generate_python(optimize=False)`.** With False, each
  bit's function is written by `BooleanFunction.generate_python`, highest bit
  first, and a subfunction already named by a higher bit is reused. With True,
  all the bits are optimized together, and if that has strictly fewer
  operations the result is written by `_straight_line_program`, which names
  every gate shared between *any* of the bits.
- **`compile`.** `BooleanFunction.compile()` compiles
  `generate_python(optimize=True)`. `FeedbackFunction.compile(optimize=True)`
  passes its argument through to `generate_python`, so by default it compiles
  the optimized update, and `compile(optimize=False)` compiles each bit as
  written. `FeedbackRegister(seed, fn, compile=True)` calls `fn.compile()`.

Two properties make True a safe default for `compile` while False is the
default for `generate_python`:

- With True, code that differs from the function as written is produced only
  when the rewrite saves at least one operation. Otherwise the statements are
  those of False, so an update the optimizer cannot improve (an LFSR, an MPR)
  compiles to exactly the code `compile(optimize=False)` gives.
- Compiled code only ever runs on arrays of bits, which is the reading the
  rewrite assumes (section 1). `generate_python` has no such guarantee about
  what its caller does with the text, and with False its output corresponds
  statement for statement to the function's own nodes, which is what a caller
  passing `overrides`, or reading the code against `fn_list`, relies on.

`fn_list` is never modified by any of these.

**Where it pays.** Section 7 shows the pattern: an update built by composing or
unrolling other functions (`iterator`), a T-function, or an update with many
products drawn from a common pool. `operation_count` before and after `optimize`
is a quick test: if the count barely moves, the clock speed will not either.
An LFSR, an MPR, or a chained CMPR as the library generates them gains nothing,
because its clock time is mostly the call into the compiled update (section 7).

`_straight_line_program` emits one Python statement per shared gate and one per
output bit: a gate object with more than one user is assigned to a local, and
every other gate is written inline in its user's expression. Three layout
choices were measured (section 7) and follow the unoptimized emitter:

- output bits are written from the highest index down, each as soon as its
  value is ready;
- a state bit is read from the array at each use, not copied to a local first;
- consequently reads and writes interleave, and the state array and the output
  array must be distinct, which `FeedbackRegister` guarantees by swapping two
  buffers.

## 7. Measurements

Measured 2026-10-01 on one Windows 11 development machine, numba-compiled, as
nanoseconds per clock in a compiled loop calling `_compiled_inplace` through a
function pointer, the way `period` and `run` do
([Compiled Evaluation](Compiled%20Evaluation.md) §3). Each row compares
`compile(optimize=False)` (left of the arrow) with `compile()` (right). Each figure is the minimum
over interleaved rounds. Run-to-run noise on this machine is about ±7% at these
durations: Fibonacci and Galois compile to byte-identical code in both columns
and still read 0.93 and 0.94.

| Update | Operations | ns per clock | Speed-up | numba compile |
|--------|-----------:|-------------:|---------:|--------------:|
| MPR 31 | 5 → 5 | 26.5 → 26.7 | 0.99× | 0.45 → 0.34 s |
| CMPR 7+5+3, `fast_template` | 37 → 30 | 24.3 → 24.2 | 1.01× | 0.46 → 0.35 s |
| CMPR 31+17+13, `arman_template(3)` | 162 → 154 | 45.4 → 44.8 | 1.01× | 1.68 → 1.53 s |
| Fibonacci 31 | 5 → 5 | 25.0 → 26.8 | 0.93× | 0.34 → 0.40 s |
| Galois 31 | 5 → 5 | 25.0 → 26.5 | 0.94× | 0.40 → 0.35 s |
| CrossJoin 31 | 69 → 58 | 30.2 → 29.0 | 1.04× | 0.73 → 0.70 s |
| FCSR(32) | 212 → 152 | 39.8 → 39.6 | 1.01× | 2.18 → 1.58 s |
| TFunction(32) | 497 → 62 | 76.3 → 25.4 | 3.01× | 4.57 → 0.68 s |
| CMPR 61 unrolled 4 clocks | 810 → 714 | 140.3 → 105.1 | 1.34× | 11.3 → 4.6 s |
| CMPR 61 unrolled 8 clocks | 1458 → 1264 | 187.2 → 161.3 | 1.16× | 16.5 → 7.9 s |
| CMPR 61 unrolled 16 clocks | 2754 → 2351 | 284.5 → 254.2 | 1.12× | 30.4 → 15.6 s |
| Dense ANF, 64 bits, 12 products per bit from a pool of 120 | 1826 → 758 | 256.2 → 112.6 | 2.28× | 116.4 → 4.2 s |

`optimize` itself took 0.5 to 50 ms on these inputs.

What the table shows:

- **Registers that are mostly shifts gain nothing.** Their per-clock time is
  dominated by a cost the logic does not control. Measured together in a
  separate run of the same loop, an update that does nothing at all took
  18.9 ns per clock and the 31-bit MPR took 23.8 ns, so about four fifths of
  that register's clock is the call and the buffer swap around it. For these
  the optimizer either finds no saving (and `compile` emits the original code)
  or saves a few operations with no measurable effect.
- **Updates with shared or redundant logic gain 1.1× to 3×.** For the two rows
  where the optimizer removes the most, the speed-up is below the operation
  ratio: TFunction(32) has 8.0× fewer operations and is 3.01× faster, and the
  dense ANF has 2.41× fewer and is 2.28× faster. A fixed per-clock call cost of
  roughly 19 ns, measured in a separate run
  ([Compiled Evaluation](Compiled%20Evaluation.md) §3), accounts for the
  direction and rough size of the difference.
- **Compile time falls**, since numba is given less code. The optimizer's own
  time is small beside it.

The layout choices in section 6 come from the same harness. Reading every state
bit into a local before writing anything measured 5 to 20% slower than the
unoptimized code on the shift-dominated registers (Fibonacci 0.80 to 0.87×,
FCSR 0.79×), and reading at each use measured level with it (0.95 to 1.06×).
The order of the output bits made no measurable difference.
