# Getting Started

This page gives worked examples of the main things PyPR does: clocking a
register, building Boolean functions, synthesizing a register from a sequence,
running an algebraic attack, controlling output, and saving objects. Every code
block here was run as written. Installation is in the [README](../README.md#setup).

## A first register

The first example builds an MPR (Mersenne product register) of size 5 with
feedback polynomial `0x12`, runs it, and finds its period. `compile=True` JIT-compiles
the update up front, so later clocking runs at numba speed.

```python
from PyPR import FeedbackRegister
from PyPR.FeedbackFunctions import MPR

m = MPR(5, "12")
reg = FeedbackRegister(1, m, compile=True)    # seed 1: bit 0 set

print([int(r[0]) for r in reg.run(12)])       # bit 0 over 12 clocks
# [1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1]

print(FeedbackRegister(1, m, compile=True).period())
# (31, 0)  -- (period, preperiod): a maximal-length 5-bit register
```

`run` yields the register itself at each clock, so `r[i]` reads bit `i` of the
current state. Pass `compiled=False` to `run`, `clock` or `period` to use the
pure-Python path instead.

Compiling optimizes the update first: logic that several bits share is computed
once, and redundant terms are simplified away. A plain shift register has
nothing to remove and compiles exactly as written; updates with a lot of shared
logic clock up to three times faster. `fn.compile(optimize=False)` skips the
optimization, and
`fn.generate_python(optimize=True)` shows the optimized code. See
[`architecture/Logic Optimization.md`](architecture/Logic%20Optimization.md).

## Boolean functions

Functions are built from gates and leaves, and every operation works on the DAG,
so shared subexpressions are evaluated and emitted once.

```python
from PyPR.BooleanLogic import AND, NOT, VAR, XOR

f = XOR(AND(VAR(0), VAR(1)), VAR(2), NOT(VAR(3)))

f.eval([1, 1, 0, 1])               # 1
f.translate_ANF().degree()         # 2

print("\n".join(f.generate_python()))
# output = ((array[0] & array[1]) ^ array[2] ^ (1-(array[3])))
print("\n".join(f.generate_c()))
# output = ((array[0] & array[1]) ^ array[2] ^ (!(array[3])));
print("\n".join(f.generate_latex()))
# f &= x_{0} x_{1} \oplus x_{2} \oplus \overline{x_{3}}
```

LaTeX output is configured through `LatexStyle`: operator symbols,
precedence, negation, variable and constant templates, and line layout. Whole
registers render as a system of equations with
`feedback_fn.generate_latex()`, and as hardware with `feedback_fn.write_VHDL(path)`.

## Register synthesis

A register can be recovered from its output. `Fibonacci.from_seq` runs
Berlekamp–Massey and returns a seed and a Fibonacci LFSR that reproduce the
sequence:

```python
from PyPR.FeedbackFunctions import Fibonacci

stream = [int(r[0]) for r in FeedbackRegister(1, m, compile=True).run(40)]
seed, lfsr = Fibonacci.from_seq(stream)

print(lfsr.size)                                                       # 5
print([int(r[0]) for r in FeedbackRegister(seed, lfsr).run(12, compiled=False)])
# [1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1]
```

`from_seq(stream, nonlinear=True)` fits an NLFSR instead, and `FCSR.from_seq`
recovers an FCSR with the 2-adic analogue of Berlekamp–Massey (rational
approximation; `BM_FCSR`).

## An algebraic attack

This example recovers the secret initial state of a small chained CMPR
(composite Mersenne product register) from its keystream, using the naive
algebraic attack (NAA). Attacks split into an offline phase, which depends only on the register and the
output function, and an online phase, which solves for one keystream.

```python
import logging
import random

import numpy as np

import PyPR
from PyPR import FeedbackRegister
from PyPR.BooleanLogic import AND, VAR, XOR
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template
from PyPR.FeedbackFunctions import CMPR, MPR
from PyPR.Cryptanalysis.Attacks.naive_algebraic_attack import NAA_offline, NAA_online

random.seed(1); np.random.seed(1)             # chaining is randomized
cmpr = CMPR([MPR(5, "12"), MPR(3, "5")])
cmpr.generate_chaining(template=arman_template(max_and=2))
output_fn = XOR(AND(VAR(0), VAR(1)), VAR(2))

PyPR.logging.level = logging.INFO             # watch the run

attack_data = NAA_offline(cmpr, output_fn, 0, 60,
                          monomial_profiles=cmpr.monomial_profiles(),
                          variable_blocks=cmpr.blocks)

secret = FeedbackRegister(42, cmpr)
keystream = np.array([output_fn.eval(r) for r in secret.copy().run(
    attack_data.keystream_needed, compiled=False)], dtype=np.uint8)

recovered = NAA_online(cmpr, output_fn, keystream, attack_data)
print([int(b) for b in recovered] == [int(b) for b in secret])   # True
```

At `INFO`, the run prints its phases nested by depth, with a timing line when
each one finishes. In a terminal, progress meters are drawn as live bars while
they run.

```
Offline phase (Naive Algebraic Attack)
|   Monomial profile optimization: on
|   Monomial profile
|   Monomial profile finished: 0.00 s
|   Variable map
|   Variable map finished: 0.00 s
|   Generating equations
|   |   Equations found: 82/82 in 0.33 s (246 eq/s)
|   |   Linear complexity reached at t=82
|   |   Keystream needed: 82 bits; unsolved monomials: 0
|   Generating equations finished: 0.33 s
Offline phase (Naive Algebraic Attack) finished: 0.48 s
Online phase (Naive Algebraic Attack)
|   Solving
|   |   LU solve
|   |   |   Collecting guess effect vectors
|   |   |   Collecting guess effect vectors finished: 0.00 s
|   |   |   Guessing
|   |   |   |   Guessing remaining state
|   |   |   |   |   Base solution is already correct
|   |   |   |   Guessing remaining state finished: 0.16 s
|   |   |   Guessing finished: 0.16 s
|   |   LU solve finished: 0.17 s
|   Solving finished: 0.17 s
Online phase (Naive Algebraic Attack) finished: 0.17 s
```

The reduced (`RAA_*`) and fast (`FAA_*`) attacks take an annihilator and its
multiple, computed with
`PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator.annihilators`.
Offline results are dataclasses (`NAAOfflineData`, `RAAOfflineData`,
`FAAOfflineData`, `CubeAttackData`). Which equation stores and solvers each
attack accepts is listed in
[`architecture/Attack_Compatibilities.md`](architecture/Attack_Compatibilities.md).

## Output

PyPR prints nothing but warnings and errors unless asked. One setting controls
how much it prints:

```python
import logging
import PyPR

PyPR.logging.level = logging.INFO     # phases, results, progress
PyPR.logging.level = logging.DEBUG    # plus per-step solver and store detail
PyPR.logging.level = "OFF"            # nothing, not even warnings
```

Output goes to standard output. Piped output gets a plain progress line every
ten seconds in place of live bars, so `python run.py > run.log` keeps a usable
record.

## Saving and loading

Registers, feedback functions, Boolean functions, ANFs and rational generating
functions all serialize to JSON:

```python
cmpr.to_file("register.json")
same = CMPR.from_file("register.json")
```

`PyPR.JSON_Serialization.generate_JSON(a, b, ...)` stores several objects in one
document and keeps any structure they share shared when they are loaded back.

## Layout

| Package | Contents |
|---------|----------|
| `PyPR.FeedbackRegister` | `FeedbackRegister`: seeding, clocking, running, period finding |
| `PyPR.FeedbackFunctions` | `MPR`, `CMPR`, `Fibonacci`, `Galois`, `FCSR`, `CrossJoin`, `TFunction` |
| `PyPR.BooleanLogic` | Boolean function DAGs, gates, `BooleanANF`, `BooleanGF`, SAT encoding, LaTeX styles, chaining templates |
| `PyPR.Tools` | Linear-complexity estimation (root expressions, monomial profiles), register synthesis, resolvents, number theory |
| `PyPR.Cryptanalysis` | Algebraic and cube attacks, and the generators, stores, annihilators and solvers they are built from |
