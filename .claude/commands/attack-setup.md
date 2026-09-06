---
description: Guided workflow for constructing end-to-end algebraic attack experiments with correct component wiring.
---

# Attack Setup

Use this skill when setting up a new algebraic attack experiment (NAA, RAA, or FAA) to ensure components are wired correctly.

## Step 1 — Choose the Attack Type

| Attack | When to Use | Key Advantage |
|--------|------------|---------------|
| **NAA** | Output function degree is already low, or as a baseline | Simplest; no annihilator computation needed |
| **RAA** | A true annihilator exists (or a low-degree pair with h=0) | Reduces system degree via annihilator |
| **FAA** | A low-degree pair (g, h) exists where h has low linear complexity | Exploits linear recurrence of the multiple for fewer equations |

## Step 2 — Select Components

Read `docs/architecture/Attack_Compatibilities.md` for the full compatibility reference. Summary:

### NAA Components
- **Generator:** `CubeEqGenerator` (with monomial profiles) or `SubstitutionEqGenerator` (dynamic)
- **Offline store:** `LUEqStore` (always — no other option in current code)
- **Solver:** `LUSolver` (default, recommended) or `GaussElimSolver`
- **NOT compatible:** `GrobnerSolver` — NAA defers keystream constants, which Grobner cannot accept separately

### RAA Components
- **Annihilator:** Compute via `SparseAnnihilator.annihilators(output_fn)` or `GaussianAnnihilator`
- **Generator:** `CubeEqGenerator` (with profiles) or `SubstitutionEqGenerator` (dynamic)
- **Offline stores:** `EqStore` for equations, linked `LUEqStore` for rank tracking
- **Online store:** `LUEqStore(comb_to_idx, consistent=True)` (default)
- **Solver:** `LUSolver` (default) — `GrobnerSolver` is untested but plausible

### FAA Components
- **Annihilator + multiple:** Same computation as RAA, but keep both g and h = f*g
- **Generator:** `CubeEqGenerator` (with profiles) or `SubstitutionEqGenerator` (dynamic)
- **Offline store:** `EqStore` for reducer equations
- **Online store:** `LUEqStore(comb_to_idx, consistent=True)` (default)
- **Solver:** `LUSolver` (default)
- **Linear relation:** Computed automatically via Berlekamp-Massey in the offline phase

## Step 3 — Construct the Register

```python
from PyPR.FeedbackFunctions import MPR, CMPR
from PyPR.BooleanLogic import AND, XOR, VAR
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template

# Choose component MPRs (sizes and primitive polynomials)
M5 = MPR(5, [1,0,1,0,0,1], [1,1,0,0,1])
M3 = MPR(3, [1,1,0,1], [1,0,1])
C = CMPR([M5, M3])
C.generateChaining(template=arman_template(max_and=2))

# Choose an output function
output_fn = XOR(AND(VAR(0), VAR(1)), VAR(2))
```

## Step 4 — Compute Annihilators (RAA/FAA only)

```python
from PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator import annihilators

(d_ann, d_mult), basis = annihilators(output_fn, verbose=True)
ann = basis[0]

# For FAA: also compute the multiple
from PyPR.BooleanLogic import AND
mult_anf = AND(output_fn, ann).translate_ANF()
mult_fn = mult_anf.to_BooleanFunction()
```

## Step 5 — Run the Attack

### NAA
```python
from PyPR.Cryptanalysis.Attacks.naive_algebraic_attack import NAA_offline, NAA_online
from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.FeedbackRegister import FeedbackRegister
import numpy as np

# Offline
attack_data = NAA_offline(C, output_fn, 0, 500,
    time_limit=30, verbose=True,
    monomial_profiles=C.monomial_profiles(), variable_blocks=C.blocks)

# Generate keystream from a secret state
secret = 42
F = FeedbackRegister(secret, C)
ks = np.array([int(output_fn.eval(reg._state))
               for reg in F.run(attack_data['keystream needed'], compiled=False)])

# Online
result = NAA_online(C, output_fn, ks, attack_data,
    verbose=True, solver=LUSolver())
```

### RAA
```python
from PyPR.Cryptanalysis.Attacks.reduced_algebraic_attack import RAA_offline, RAA_online

attack_data = RAA_offline(C, ann, mult_fn, 0, 500,
    time_limit=30, verbose=True,
    monomial_profiles=C.monomial_profiles(), variable_blocks=C.blocks)

result = RAA_online(C, output_fn, ks, attack_data, verbose=True)
```

### FAA
```python
from PyPR.Cryptanalysis.Attacks.fast_algebraic_attack import FAA_offline, FAA_online

attack_data = FAA_offline(C, ann, mult_fn, 0, 500,
    time_limit=30, verbose=True,
    monomial_profiles=C.monomial_profiles(), variable_blocks=C.blocks)

result = FAA_online(C, output_fn, ks, attack_data, verbose=True)
```

## Step 6 — Verify the Result

```python
# Check if the recovered state matches the secret
F_check = FeedbackRegister(secret, C)
print("Secret state:", F_check._state)
print("Recovered:   ", result)  # compare
```

## Common Mistakes

- **Forgetting `compiled=False`** when generating keystream in probe scripts.
- **Using `GrobnerSolver` with NAA** — this will raise `ValueError`.
- **Wrong keystream length** — always use `attack_data['keystream needed']`, not a hardcoded value.
- **Mismatched output function** — the output_fn passed to the online phase must be the same one used in the offline phase.
- **Skipping monomial profiles** — without `monomial_profiles` and `variable_blocks`, the attack falls back to the slower dynamic path. For experiments where speed matters, always compute and pass these.
