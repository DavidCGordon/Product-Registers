---
description: Comprehensive API reference for all PyPR classes and functions — method signatures, return types, what generators yield.
---

# PyPR API Reference

Quick lookup for method signatures, constructor arguments, and generator behavior across all major PyPR modules.

---

## FeedbackRegister

```python
from PyPR.FeedbackRegister import FeedbackRegister

FeedbackRegister(
    seed: int | list[int] | np.ndarray,  # int: binary, LSB=index 0
    fn: FeedbackFunction,
    compile: bool = False   # if True, calls fn.compile() eagerly
)
```

**Key attributes:** `size: int`, `fn: FeedbackFunction`, `_state: np.ndarray[uint8]`, `_seed: np.ndarray[uint8]`

| Method | Signature | Notes |
|--------|-----------|-------|
| `seed(s)` | `int\|list\|ndarray → None` | Resets the stored seed; does not change current state |
| `reset()` | `→ None` | Copies `_seed` back into `_state` |
| `set_state(s)` | `int\|list\|ndarray → None` | Overwrites current state (copies ndarray in) |
| `clock(compiled=True)` | `→ None` | One cycle. Raises if compiled and fn not compiled |
| `run(limit=None, compiled=True)` | `→ Iterator[self]` | Yields self n times (t=0..n-1); register lands at t=n after |
| `period(compiled=True, safe=True, limit=None)` | `→ tuple[int,int]\|None` | `(period, preperiod)`. Does not modify state. Brent (safe) or naive |
| `copy()` | `→ FeedbackRegister` | Deep copy |
| `to_JSON()` / `from_JSON(d)` | | Round-trip through dict |
| `to_file(path)` / `from_file(path)` | | JSON file I/O (`.json` extension required) |

**`run()` semantics:**
- Yields `self` — the register object, not a copy. `_state` changes every cycle.
- `run(n)` yields states t=0 through t=n-1. After the loop, register is at t=n.
- Two sequential `run(a)` / `run(b)` cover t=0..a+b-1 with no gaps or repeats.
- `run()` with no limit is infinite — break manually.

**Period defaults (limit):** compiled+safe=2²⁵, compiled+unsafe=2²⁶, uncompiled+safe=2¹⁴, uncompiled+unsafe=2¹⁶

---

## FeedbackFunction (base + subclasses)

All subclasses share the base interface.

```python
from PyPR.FeedbackFunctions import MPR, CMPR, Fibonacci, Galois, CrossJoin, TFunction, FCSR
```

**Base attributes:** `fn_list: list[BooleanFunction]`, `size: int`

| Method | Notes |
|--------|-------|
| `compile()` | Numba-JIT the function. Must call before `compiled=True` in register |
| `copy()` | Deep copy |
| `iterator(n)` | Yields `list[BooleanFunction]` — one per bit, symbolically unrolled for n steps |
| `gateSummary()` | `→ dict[str,int]` — counts of AND/XOR/NOT/VAR/CONST nodes |
| `isLinear(allowAffine=False)` | `→ bool` |
| `dense_str()` / `pretty_str()` / `anf_str()` | String representations |
| `to_JSON()` / `from_JSON(d)` / `to_file(p)` / `from_file(p)` | Serialization |

### MPR

```python
MPR(
    size: int,
    primitive_poly: str | list[int],  # list: [c0,c1,...,cn], coeff of x^i at index i
    update_poly: list[int] | None = None   # length n; default multiplies by x
)
```

`primitive_poly` as Koopman hex string is also accepted.

**Attributes:** `primitive_polynomial: list[int]`, `update_polynomial: list[int]`  
**Cached:** `minimal_polynomial: list[int]`

### CMPR

```python
CMPR(components: list[MPR | CMPR | FeedbackFunction])
# components ordered high-to-low block; nested CMPRs are flattened
```

**Key attributes:** `num_components`, `blocks`, `divisions`, `primitive_polynomials`, `update_polynomials`

| Method / Property | Notes |
|-------------------|-------|
| `generateChaining(template)` | Apply chaining template (arman_template, fast_template, etc.) |
| `update_MPR(mpr_index, new_update_poly)` | Replace one block's update poly (invalidates cached props) |
| `blocks` | `@cached_property → list[list[int]]` — bit indices per block, high-to-low |
| `divisions` | Bit-index boundaries between blocks |
| `max_period` | `@cached_property → int` — LCM of all block periods |
| `expected_period` | `@cached_property → float` |
| `cycle_lengths` | `@cached_property → list[tuple[int,int]]` — (period, multiplicity) |
| `monomial_profiles(verbose=False, force_default=False)` | `→ list[MonomialProfile]` — one per bit |
| `root_expressions(locked_list=None, verbose=False, force_default=False)` | `→ list[RootExpression]` — one per bit |
| `estimate_LC(output_bit, locked_list=None, verbose=False)` | `→ tuple[int,int]` — (lower, upper) bounds |
| `fixpoint` | `@property → list[int]` — state where F(s)=s |
| `reverse_clock(state)` | `→ list[int]` — predecessor state |
| `update_matrices` | `@cached_property → list[ndarray]` — GF(2) matrix per block |
| `resolvent_matrices` | `@cached_property → list[ndarray]` — (I+UD)^{-1} per block |
| `has_chaining` | `@property → list[int]` — number of chaining terms per bit |
| `component_feedback` | `@property → list[BooleanFunction]` — linear (MPR) part per bit |
| `chaining_feedback` | `@property → list[BooleanFunction]` — nonlinear part per bit |

### Fibonacci / Galois

```python
Fibonacci(size: int, primitive_polynomial: str | list[int])
Galois(size: int, primitive_polynomial: str | list[int])
```

| Class method | Notes |
|-------------|-------|
| `fromSeq(seq, nonlinear=False)` | `→ (seed, register)` from output sequence (BM internally) |
| `fromReg(F, bit=0, numIters=None, nonlinear=False)` | `→ (seed, register)` from running register |
| `invert()` | Toggle time-reversed mode |

**Attribute:** `update_matrix: @cached_property ndarray`

### CrossJoin

```python
CrossJoin(size: int, primitive_poly: str | list[int])
```

| Method | Notes |
|--------|-------|
| `generateNonlinearity(maxAnds=4, tapDensity=0.75)` | Add random AND terms |
| `addNonLinearTerm(maxAnds)` | Add one random AND term |
| `blocks` | `@property → list[list[int]]` — single block |
| `monomial_profiles()` | `→ list[MonomialProfile]` |
| `root_expressions()` | `→ list[RootExpression]` |
| `filter_generator()` | `→ (Fibonacci, tuple[BooleanFunction,...])` — base LFSR + filters |
| `convert_state(state)` | Convert CrossJoin state to base LFSR state |
| `compensation_list()` | Filter functions for base LFSR |

### TFunction

```python
TFunction(size: int)   # n-bit binary counter (triangular dependencies)
```

Inherits CMPR. Extra attribute: `induction_order: list[int]`

### FCSR

```python
FCSR(diadic_complexity: int, q: int)   # q must be odd positive
# size = 2*diadic_complexity - 1 (interleaved values + carries)
```

| Class method | Notes |
|-------------|-------|
| `fromSeq(seq)` | `→ (seed, register)` from sequence (2-adic BM) |
| `fromReg(F, bit=0, numIters=None)` | `→ (seed, register)` from running register |
| `state_from_frac(num, den)` | `→ (size, state)` from 2-adic fraction |

**Properties:** `carries → list[BooleanFunction]`, `values → list[BooleanFunction]`

---

## BooleanFunction and Gates

```python
from PyPR.BooleanLogic import AND, XOR, OR, NOT, VAR, CONST
from PyPR.BooleanLogic.BooleanFunction import BooleanFunction
```

**Constructors:**
```python
VAR(index: int)          # leaf: reads state[index]
CONST(value: Any)        # leaf: returns value (typically 0 or 1)
AND(*args)               # n-ary AND
XOR(*args)               # n-ary XOR
OR(*args)                # n-ary OR
NOT(arg)                 # unary NOT (arg_limit=1)
NAND(*args)
NOR(*args)
XNOR(*args)
```

**Core evaluation:**

| Method | Notes |
|--------|-------|
| `eval(state, cache=None)` | `→ np.uint8`. state is any indexable container. Cache is dict for subexpression reuse |
| `eval_ANF(array, cache=None)` | `→ BooleanANF`. array is list of BooleanANF; evaluates symbolically |
| `translate_ANF()` | `→ BooleanANF`. Convert to Algebraic Normal Form |

**Composition:**

| Method | Notes |
|--------|-------|
| `compose(input_map, in_place=False)` | Replace VAR(i) with input_map[i] (any indexable container of BooleanFunction) |
| `remap_indices(index_map, in_place=False)` | Remap VAR indices |
| `shift_indices(shift_amount, in_place=False)` | Add shift to all VAR indices |
| `remap_constants(const_map, in_place=False)` | Replace constant values; const_map is `list[tuple[old_val, new_val]]` |

**Structural:**

| Method | Notes |
|--------|-------|
| `copy()` | Deep copy |
| `subfunctions()` | `→ list[BooleanFunction]` — reused subexpressions, topologically sorted |
| `inputs()` | `→ list[BooleanFunction]` — all VAR/CONST leaves |
| `idxs_used()` | `→ set[int]` — VAR indices that appear |
| `max_idx()` | `→ int` — highest VAR index |
| `is_leaf()` | `→ bool` |
| `gateSummary()` | `→ dict[str,int]` — node counts by type |
| `num_nodes()` | `→ int` |
| `merge_redundant(in_place=False)` | Remove duplicate subexpressions |
| `binarize(in_place=False)` | Convert to binary (2-arg) gates |

**String / code generation:**

| Method | Notes |
|--------|-------|
| `dense_str()` | Single-line nested form |
| `pretty_str()` | Multiline with subfunctions listed |
| `generate_python(output_name, subfunction_prefix, array_name, overrides)` | `→ list[str]` |
| `generate_c(...)` | `→ list[str]` |
| `generate_VHDL(...)` | `→ list[str]` |

**SAT:**

| Method | Notes |
|--------|-------|
| `tseytin_encode()` | `→ (clauses, node_labels, var_labels)` — CNF encoding |
| `solve_SAT(assumptions=None, ...)` | `→ (bool, dict) \| None` — satisfying assignment |

---

## BooleanANF

```python
from PyPR.BooleanLogic.BooleanANF import BooleanANF

BooleanANF(nested_iterable=None, fast_init=False)
# BooleanANF()          → zero function
# BooleanANF([True])    → constant 1
# BooleanANF([[0,1]])   → x_0 * x_1
# BooleanANF([[0],[1]]) → x_0 XOR x_1
```

Terms are `frozenset[frozenset]`. Empty frozenset = constant 1 term.

| Operator / Method | Notes |
|-------------------|-------|
| `__xor__(other)` / `__add__(other)` | XOR: symmetric difference of term sets |
| `__and__(other)` / `__mul__(other)` | AND: Cartesian product of terms |
| `__invert__()` | NOT: XOR with constant 1 |
| `degree()` | `→ int` — max term length |
| `__len__()` | Number of terms |
| `__iter__()` | Iterate over terms (each is a frozenset of variable indices) |
| `__eq__(other)` | Term-set equality |
| `from_BooleanFunction(bf)` | `@classmethod → BooleanANF` |
| `to_BooleanFunction()` | `→ BooleanFunction` |
| `dense_str()` | String representation |

**Class methods:**
```python
BooleanANF.from_BooleanFunction(bf)   # convert BooleanFunction → BooleanANF
```

---

## Chaining Templates

```python
from PyPR.BooleanLogic.ChainingGeneration.Templates import arman_template, fast_template

# arman_template: structured random AND terms between blocks
template_fn = arman_template(max_and: int = 4)
C.generateChaining(template=template_fn)

# fast_template: simpler random AND terms
template_fn = fast_template(max_and: int = 4)
```

Templates are callables: `template(cmpr: CMPR) → dict[int, BooleanFunction]`

---

## Annihilators

```python
from PyPR.Cryptanalysis.Components.Annihilators.SparseAnnihilator import annihilators
from PyPR.Cryptanalysis.Components.Annihilators.GaussianAnnihilator import annihilators as annihilators_gaussian

(d_ann, d_mult), basis = annihilators(output_fn, verbose=False)
# d_ann: optimal annihilator degree
# d_mult: optimal (f * annihilator) degree
# basis: list[BooleanFunction] — basis for the annihilator space

(d_ann, d_mult), basis = annihilators_gaussian(output_fn)
# Same return format; uses Gaussian elimination instead of sparse algorithm
```

Both algorithms find the minimum-degree annihilator. Sparse is generally faster for low-degree functions; Gaussian is more straightforward.

---

## Equation Generators

All generators are **functions** (not classes) returning iterators.

### SymbolicEqGenerator

```python
from PyPR.Cryptanalysis.Components.EquationGenerators.SymbolicEqGenerator import SymbolicEqGenerator

SymbolicEqGenerator(
    feedback_fn: FeedbackFunction,
    output_fn: BooleanFunction | list[BooleanFunction],
    limit: int,
    initialization=None   # list[BooleanANF] for partial state knowledge; None = fresh symbolics
) → Iterator[BooleanFunction | list[BooleanFunction]]
```

- Yields `limit+1` equations: t=0 through t=limit (inclusive).
- With single `output_fn`: yields `BooleanFunction` each step.
- With `list[BooleanFunction]`: yields `list[BooleanFunction]` each step.
- `initialization`: substitute known variables (e.g. replace VAR(i) with a constant).
- Slowest generator but most flexible; supports partial initialization.

### CubeEqGenerator

```python
from PyPR.Cryptanalysis.Components.EquationGenerators.CubeEqGenerator import CubeEqGenerator

CubeEqGenerator(
    feedback_fn: FeedbackFunction,
    output_fn: BooleanFunction | list[BooleanFunction],
    limit: int,
    monomial_profiles: list[MonomialProfile] | None = None,
    variable_blocks: list[list[int]] | None = None,
    include_variables: bool = True,
    complete_subsets: bool = False,
    include_constant: bool = True,
    time_limit: float | None = None,
    verbose: bool = False,
    print_depth: int = 0
) → Iterator[...]
```

- Fastest generator (~1000x faster than Symbolic).
- Requires `monomial_profiles` and `variable_blocks` for full functionality.
- Use `cmpr.monomial_profiles()` and `cmpr.blocks` to get these.

### SubstitutionEqGenerator

```python
from PyPR.Cryptanalysis.Components.EquationGenerators.SubstitutionEqGenerator import SubstitutionEqGenerator

SubstitutionEqGenerator(
    feedback_fn: FeedbackFunction,
    output_fn: BooleanFunction | list[BooleanFunction],
    limit: int,
    initialization=None
) → Iterator[BooleanFunction | list[BooleanFunction]]
```

- Similar to Symbolic but composes in opposite direction (faster when few bits feed the output).
- Cannot accept partial initialization.

---

## Equation Stores

### LUEqStore (most common)

```python
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

LUEqStore(
    comb_to_idx: dict[tuple[int,...], int] | None = None,  # None → dynamic (grows)
    consistent: bool = False   # True → check for contradictions (needs CONST column)
)
```

**Key attributes:** `rank` (= `num_eqs`), `num_vars`, `comb_to_idx`, `idx_to_comb`, `upper_matrix`, `lower_matrix`, `solved_for`

| Method | Notes |
|--------|-------|
| `insert_equation(eq, identifier=None, translate_ANF=True)` | `→ bool` (True if linearly independent / new info) |

- Dynamic store: grows `comb_to_idx` as new monomials appear. Accepts `BooleanFunction`.
- Static store: fixed index map. Can also accept `np.ndarray` coefficient vectors.
- `identifier`: metadata (clock time, cube) stored for debugging; does not affect solving.

### Other Stores

```python
from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore         # dense matrix, no reduction
from PyPR.Cryptanalysis.Components.EquationStores.SymbolicEqStore import SymbolicEqStore  # Groebner-based
from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GrobnerEqStore
from PyPR.Cryptanalysis.Components.EquationStores.IndexedEqStore import IndexedEqStore    # base class
```

All share `insert_equation(eq, identifier=None, translate_ANF=True)` interface.

---

## Equation Solving

```python
from PyPR.Cryptanalysis.Components.EquationSolving.LU_Solver import LUSolver
from PyPR.Cryptanalysis.Components.EquationSolving.GaussElim import GaussElimSolver
from PyPR.Cryptanalysis.Components.EquationSolving.Grob_Solver import GrobnerSolver
from PyPR.Cryptanalysis.Components.EquationSolving.GuessSolver import guess_and_solve
```

Each solver class has two methods:
- `reduce(equation_store, ...)` — applies the solver's reduction (LU back-substitution, RREF, or Gröbner basis) without guessing free variables.
- `solve(equation_store, feedback_fn, output_fn, keystream, *, ...)` — full pipeline: `reduce` + exhaustive guess-and-prune via `GuessSolver`. Derives `guess_bits` and `variable_indices` internally from the store. `LUSolver` is the most efficient (reuses cached decomposition); `GaussElimSolver` re-derives effect vectors from the RREF.

Solver configuration is set at construction:
- `LUSolver(additional_constants=None)` / `GaussElimSolver(additional_constants=None)` — `additional_constants` is set by NAA internally
- `GrobnerSolver(simplify_mode=None)` — Gröbner simplification strategy

Module-level `reduce()` / `solve()` functions remain available alongside the classes.

---

## Attacks

### Naive Algebraic Attack (NAA)

```python
from PyPR.Cryptanalysis.Attacks.naive_algebraic_attack import NAA_offline, NAA_online

attack_data = NAA_offline(
    feedback_fn, output_fn, init_rounds,
    time_limit,             # seconds for offline phase
    verbose=False,
    print_depth=0,
    monomial_profiles=None,  # from cmpr.monomial_profiles()
    variable_blocks=None     # from cmpr.blocks
)
# attack_data keys: 'keystream needed', 'guess vars', 'upper matrix', 'lower matrix',
#                   'idx_to_comb', 'comb_to_idx', 'equation times'

result = NAA_online(
    feedback_fn, output_fn,
    keystream,      # np.ndarray of output bits
    attack_data,
    test_length=1000,
    verbose=False,
    print_depth=0,
    solver=LUSolver(),         # or GaussElimSolver() (GrobnerSolver not supported for NAA)
)
# Solver classes: LUSolver, GaussElimSolver, GrobnerSolver
# from PyPR.Cryptanalysis.Components.EquationSolving.{LU_Solver,GaussElim,Grob_Solver}
```

### Fast Algebraic Attack (FAA)

```python
from PyPR.Cryptanalysis.Attacks.fast_algebraic_attack import FAA_offline, FAA_online

# First compute annihilators:
(d_ann, d_mult), basis = annihilators(output_fn)
annihilator = basis[0]
multiple = AND(output_fn, annihilator).translate_ANF()   # f * g → BooleanANF → back to BooleanFunction
multiple_fn = multiple.to_BooleanFunction()

attack_data = FAA_offline(
    feedback_fn, annihilator, multiple_fn,
    init_rounds, max_time,
    time_limit=120,
    verbose=False,
    monomial_profiles=None,
    variable_blocks=None
)
# attack_data keys: 'annihilator equations', 'idx to comb map', 'keystream needed', ...

comb_to_idx = attack_data['comb to idx map']
result = FAA_online(
    feedback_fn, output_fn,
    keystream, attack_data,
    verbose=False,
    solver=LUSolver(),                              # or GaussElimSolver(), GrobnerSolver()
    online_store=LUEqStore(comb_to_idx, consistent=True),  # or EqStore(comb_to_idx), GroebnerEqStore()
)
```

### Reduced Algebraic Attack (RAA)

```python
from PyPR.Cryptanalysis.Attacks.reduced_algebraic_attack import RAA_offline, RAA_online
# Same online interface as FAA (solver, online_store)
```

---

## Tools

### Berlekamp-Massey

```python
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

lc, poly = berlekamp_massey(seq)  # seq: list[int] or np.ndarray
# lc: linear complexity (int)
# poly: coefficient array of minimal polynomial
```

### RootExpression

```python
from PyPR.Tools.RootCounting.RootExpression import RootExpression

re = cmpr.root_expressions()[i]   # get for bit i
re.upper()   # → int: upper bound on linear complexity
re.lower()   # → int: lower bound

# Arithmetic:
re1 + re2    # union (XOR of sequences)
re1 * re2    # product (AND/convolution of sequences)
re.extend(jordan_set)  # add roots with multiplicity
```

### Cost Estimation

```python
from PyPR.Tools.CostEstimation import estimate_cost_comp, estimate_cost_cube

profiles = cmpr.monomial_profiles()
comp_cost = estimate_cost_comp(feedback_fn, output_fn, profiles)
cube_cost = estimate_cost_cube(feedback_fn, output_fn, profiles)
# Both return int (coefficient flip count)
```

### BooleanGF

```python
from PyPR.BooleanLogic.BooleanGF import BooleanGF

D   = BooleanGF.delay()       # delay operator
one = BooleanGF.one()
zero = BooleanGF.zero()
x   = BooleanGF.from_int(3)   # integer value as polynomial

# Arithmetic: +, *, /, ** all work; returns BooleanGF
bg = D / (D + one)
bg.simplify()

# From sequence:
bg = BooleanGF.from_seq(seq)   # BM on integer sequence → rational poly
```

---

## Common Patterns

### Run + collect stream

```python
F = FeedbackRegister(1, M3)
output_fn = VAR(0)
stream = np.array([int(output_fn.eval(reg._state)) for reg in F.run(100, compiled=False)])
```

### Chained sequential runs

```python
F = FeedbackRegister(1, M3)
# Process first 50 states
for reg in F.run(50, compiled=False):
    process(reg._state.copy())
# Continue from state 50
for reg in F.run(50, compiled=False):
    process(reg._state.copy())
```

### Symbolic equation pipeline

```python
M3 = MPR(3, [1,1,0,1], [1,0,1])
output_fn = VAR(0)
store = LUEqStore()

for eq in SymbolicEqGenerator(M3, output_fn, limit=15):
    store.insert_equation(eq)

print(f"rank {store.rank} / {store.num_vars} monomials")
```

### Full FAA experiment (small scale)

```python
M5 = MPR(5, [1,0,1,0,0,1], [1,1,0,0,1])
M3 = MPR(3, [1,1,0,1], [1,0,1])
C = CMPR([M5, M3])
C.generateChaining(template=arman_template(max_and=2))
output_fn = XOR(AND(VAR(0), VAR(1)), VAR(2))

(d_ann, d_mult), basis = annihilators(output_fn)
ann = basis[0]
mult_anf = AND(output_fn, ann).translate_ANF()
mult_fn = mult_anf.to_BooleanFunction()

attack_data = FAA_offline(C, ann, mult_fn, 0, 500,
    time_limit=30, verbose=True,
    monomial_profiles=C.monomial_profiles(), variable_blocks=C.blocks)

secret = 42
F = FeedbackRegister(secret, C)
ks = np.array([int(output_fn.eval(reg._state)) for reg in F.run(attack_data['keystream needed'], compiled=False)])

result = FAA_online(C, output_fn, ks, attack_data, verbose=True)
```

### Enumerate CMPR bit structure

```python
C = CMPR([M7, M5, M3])
for b_idx, block in enumerate(C.blocks):
    print(f"Block {b_idx} (n={len(block)}): bits {block}")
```

### CrossJoin with nonlinearity

```python
from PyPR.FeedbackFunctions import CrossJoin
cj = CrossJoin(7, [1,0,0,0,0,1,1,1])   # primitive degree-7 poly
cj.generateNonlinearity(maxAnds=3, tapDensity=0.5)
F = FeedbackRegister(1, cj)
print(F.period(compiled=False))
```

### BooleanANF arithmetic

```python
from PyPR.BooleanLogic.BooleanANF import BooleanANF

x0 = BooleanANF([[0]])   # x_0
x1 = BooleanANF([[1]])   # x_1
one = BooleanANF([True]) # constant 1

f = x0 * x1 + x0 + one   # x_0*x_1 XOR x_0 XOR 1
print(f.degree(), len(f))
```

---

## Serialization

All major classes support the same protocol:

```python
obj.to_JSON()            # → dict  
Class.from_JSON(d)       # → instance
obj.to_file("path.json") # writes JSON (must end .json)
Class.from_file("path.json")

# Shared object graph (multiple objects sharing subgraphs):
ids = obj1.generate_ids()
ids = obj2.generate_ids(ids)   # reuses shared nodes
```
