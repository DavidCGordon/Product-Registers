---
description: Guide for writing and running PyPR experiment scripts to test hypotheses instead of reasoning through complex math.
---

# PyPR Probe

Use this skill whenever you are **uncertain about how a PyPR mathematical object behaves**. Do NOT try to mentally trace through annihilator theory, ANF manipulations, period computations, or monomial profile arithmetic. Run an experiment instead.

## When to invoke

- Unsure what period/ANF/output stream a register configuration produces
- Checking whether a property holds (annihilator correctness, degree, etc.)
- Verifying an algebraic claim before writing it into a docstring
- Exploring how an API argument changes mathematical output
- Sanity-checking that a code change preserved mathematical behavior

## Workflow

1. State the specific claim or question.
2. Choose the **smallest register** that exercises the property (default: M3, period 7).
3. Write a minimal script that **prints** the value — prefer `print` over `assert` so you see actuals.
4. Write to `experiments/scratch/_probe.py` (scratch file, overwrite freely).
5. Run: `python experiments/scratch/_probe.py`
6. Read the output and update your reasoning before continuing.

If it takes >30s, switch to a smaller register or add a `limit=` parameter.

---

## Standard Imports

```python
from PyPR.BooleanLogic import AND, XOR, OR, NOT, VAR, CONST
from PyPR.BooleanLogic.BooleanFunction import BooleanFunction
from PyPR.BooleanLogic.BooleanANF import BooleanANF
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, CMPR, Fibonacci, Galois, CrossJoin, TFunction, FCSR

from PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator import annihilators
from PyPR.Cryptanalysis.Components.Annihilators.GaussianAnnihilator import annihilators as annihilators_gaussian
from PyPR.Cryptanalysis.Components.EquationGenerators.SymbolicEqGenerator import SymbolicEqGenerator
from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import CubeEqGenerator
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import SubstitutionEqGenerator
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
from PyPR.Cryptanalysis.Attacks.fast_algebraic_attack import FAA_offline, FAA_online
from PyPR.Cryptanalysis.Attacks.naive_algebraic_attack import NAA_offline, NAA_online

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey
from PyPR.Tools.RootCounting.RootExpression import RootExpression
from PyPR.Tools.CostEstimation import estimate_cost_comp, estimate_cost_cube

import numpy as np
```

---

## Canonical Small Objects

Always start with the smallest object that still exercises the property.

```python
# Fastest: M3 — period 7, good for ANF/annihilator/stream checks
M3 = MPR(3, [1,1,0,1], [1,0,1])

# Medium: M5 — period 31
M5 = MPR(5, [1,0,1,0,0,1], [1,1,0,0,1])

# Larger: M7 — period 127
M7 = MPR(7, [1,1,0,0,0,0,0,1], [1,0,0,0,0,1,0])

# Small CMPR: M5+M3, period = lcm(31,7) = 217
C = CMPR([M5, M3])
C.generateChaining(template=arman_template(max_and=2))

# Register — always pass compiled=False in probe scripts
F = FeedbackRegister(1, M3)
```

---

## run() — Critical Behavior

`run()` **yields `self` (the register object) on each cycle** — not a copy. The state array is shared:

```python
F = FeedbackRegister(1, M3)

# run(n) yields n states: t=0 through t=n-1
# After the loop, the register is at state t=n (not yielded)
for reg in F.run(7, compiled=False):
    print(reg._state.copy())   # copy if you need to keep it; reg._state is overwritten next cycle

# Sequential run() calls chain seamlessly — no states are skipped or repeated
for reg in F.run(3, compiled=False):   # yields t=0,1,2; register lands at t=3
    pass
for reg in F.run(4, compiled=False):   # yields t=3,4,5,6 from where we left off
    pass

# Extracting a stream (most common pattern):
F = FeedbackRegister(1, M3)
output_fn = VAR(0)
stream = [int(output_fn.eval(reg._state)) for reg in F.run(20, compiled=False)]
print(stream)

# Collecting states:
F = FeedbackRegister(1, M3)
states = [reg._state.copy() for reg in F.run(10, compiled=False)]
```

---

## Composition Patterns

### BooleanFunction.compose() — substitute variables

```python
f = XOR(AND(VAR(0), VAR(1)), VAR(2))
g = AND(VAR(3), VAR(4))

# Replace VAR(0) → g, VAR(1) → VAR(1), VAR(2) → VAR(2)
composed = f.compose([g, VAR(1), VAR(2)])   # input_map[i] replaces VAR(i)
print(composed.dense_str())
```

### BooleanFunction.eval_ANF() — symbolic substitution

```python
# Replace VAR(i) with arbitrary BooleanANF objects and evaluate symbolically
# This is how SymbolicEqGenerator propagates state expressions
fns = [BooleanANF([[b]]) for b in range(M3.size)]   # symbolic vars x_0, x_1, x_2

output_fn = XOR(AND(VAR(0), VAR(1)), VAR(2))
result_anf = output_fn.eval_ANF(fns)   # returns BooleanANF
print(result_anf)
```

### SymbolicEqGenerator — propagated equations

```python
# Yields limit+1 BooleanFunction objects: equations at t=0 through t=limit
# Each equation is in terms of the initial state variables
for t, eq in enumerate(SymbolicEqGenerator(M3, VAR(0), limit=6)):
    print(f"t={t}: {eq.dense_str()}")
    print(f"  evaluates to: {int(eq.eval([1,0,0]))}")   # check against actual output
```

### FeedbackFunction.iterator() — symbolic unrolling

```python
# Yields list of BooleanFunction objects (one per bit) at each time step
# Each list represents the full state at time t in terms of t=0
for t, fns in enumerate(M3.iterator(5)):
    vals = [int(f.eval([1,0,0])) for f in fns]   # evaluate on a specific state
    print(f"t={t}: {vals}")
```

---

## Equation Pipeline: Generator → Store → Inspect

```python
# Build a store from symbolic equations
M3 = MPR(3, [1,1,0,1], [1,0,1])
output_fn = VAR(0)
store = LUEqStore()   # dynamic (grows as new monomials appear)

for eq in SymbolicEqGenerator(M3, output_fn, limit=10):
    was_new = store.insert_equation(eq, translate_ANF=True)
    print(f"rank={store.rank}, num_vars={store.num_vars}, new={was_new}")
```

### Inspect the equation store

```python
print("rank:", store.rank)           # number of linearly independent equations
print("num_vars:", store.num_vars)   # number of monomials tracked
print("idx_to_comb:", store.idx_to_comb)   # column index → monomial tuple
# store.upper_matrix[:store.rank, :store.num_vars]  — the pivot rows
```

---

## CMPR Structure

```python
C = CMPR([M7, M5, M3])
C.generateChaining(template=arman_template(max_and=2))

print("total size:", C.size)                   # 7+5+3 = 15
print("num_components:", C.num_components)     # 3
print("blocks:", C.blocks)                     # [[14..8], [7..3], [2..0]] — bit indices per block
print("divisions:", C.divisions)               # boundary indices

# Monomial profiles (what monomials appear in each bit's output ANF)
profiles = C.monomial_profiles()
for i, p in enumerate(profiles):
    print(f"bit {i}: {p}")

# Root expressions (linear complexity bounds per bit)
res = C.root_expressions()
for i, re in enumerate(res):
    print(f"bit {i}: lower={re.lower()}, upper={re.upper()}")

# Max achievable period
print("max_period:", C.max_period)
print("expected_period:", C.expected_period)

# Cost estimation (compare equation generation strategies)
output_fn = XOR(AND(VAR(0), VAR(1)), VAR(2))
profiles = C.monomial_profiles()
comp_cost = estimate_cost_comp(C, output_fn, profiles)
cube_cost = estimate_cost_cube(C, output_fn, profiles)
print(f"composition cost: {comp_cost}, cube cost: {cube_cost}, ratio: {comp_cost/cube_cost:.2f}")
```

---

## Berlekamp-Massey and Register Recovery

```python
import numpy as np
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

F = FeedbackRegister(1, M3)
stream = np.array([int(VAR(0).eval(reg._state)) for reg in F.run(30, compiled=False)])
lc, poly = berlekamp_massey(stream)
print(f"linear complexity: {lc}, poly: {poly}")

# Recover as Fibonacci LFSR from stream
seed, fib = Fibonacci.fromSeq(stream.tolist())
print(f"recovered seed: {seed}")
print(f"recovered poly: {fib.primitive_polynomial}")
```

---

## Truth Table Checks

```python
# Exhaustive check over 2^n inputs for small n
f = XOR(AND(VAR(0), VAR(1)), VAR(2))

for b0 in [0,1]:
    for b1 in [0,1]:
        for b2 in [0,1]:
            out = int(f.eval([b0, b1, b2]))
            print(f"[{b0},{b1},{b2}] → {out}")
```

---

## Annihilator Checks

```python
f = XOR(AND(VAR(0), VAR(1)), VAR(2))
(d_ann, d_mult), basis = annihilators(f, verbose=False)
print(f"ann degree: {d_ann}, mult degree: {d_mult}, basis size: {len(basis)}")

for ann in basis:
    anf = ann.translate_ANF()
    print("  ANF:", anf.dense_str())
    # Verify: f * ann = 0 on all inputs
    product = AND(f, ann)
    ok = all(int(product.eval([b0,b1,b2])) == 0
             for b0 in [0,1] for b1 in [0,1] for b2 in [0,1])
    print("  annihilates correctly:", ok)
```

---

## Serialization Round-Trip

```python
import json
C = CMPR([M5, M3])
C.generateChaining(template=arman_template(max_and=2))

json_data = C.to_JSON()
C2 = CMPR.from_JSON(json_data)
print("round-trip identical:", C.dense_str() == C2.dense_str())
```

---

## Gotchas

- **`run()` / `clock()` default `compiled=True`** — raises `ValueError` unless `fn.compile()` was called first. Always pass `compiled=False` in probe scripts.
- **State is mutable and shared** — `reg._state` in a `run()` loop points to the register's live array. It is overwritten each cycle. Use `reg._state.copy()` to persist a snapshot.
- **`run(n)` convention** — yields n states (t=0 through t=n-1). Register is at t=n after the loop. Two successive `run(a)` + `run(b)` calls cover t=0 to t=a+b-1 without overlap.
- **`seed()` with ndarray** — stores the reference without copying. Pass a copy (`arr.copy()`) if the array might change.
- **`eval()` returns `np.uint8`** — use `int(val)` when printing or comparing if you want clean output.
- **`SymbolicEqGenerator(fn, output_fn, limit)` yields `limit+1` equations** — t=0 through t=limit inclusive.
- **`translate_ANF=True`** in `insert_equation` — always keep this unless you know the equation is already in ANF form.
- **ANF display quirk** — single-variable monomials display as `AND(VAR(i))` (one-arg AND = identity). This is correct.
- **`period()` does not modify state** — it uses internal copies. Safe to call on a running register.
- **Dynamic vs static `LUEqStore`** — `LUEqStore()` with no args is dynamic (variable set grows). `LUEqStore(comb_to_idx)` is static and won't accept new monomials.
- **Import path for new modules** — when unsure, grep `src/PyPR/` for the class name; the file path maps to the module path.
