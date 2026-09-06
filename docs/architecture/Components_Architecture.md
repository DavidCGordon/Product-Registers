# Cryptanalysis/Components Architecture

`src/PyPR/Cryptanalysis/Components/` has five subdirectories: `EquationStores/`,
`Adapters/`, `Annihilators/`, `EquationGenerators/`, `EquationSolving/`. This
doc explains what each one is and — more importantly — how they interoperate:
which representations feed into which, what's native/fast vs. converted/slow,
and the one hard compatibility restriction that exists between an attack and
a solver.

For attack-level usage (offline/online phase structure, which generator/store
each attack picks), see [Algebraic Attacks](../theory/Algebraic%20Attacks.md) and
[Attack Compatibilities](Attack_Compatibilities.md) (attack-by-attack store/solver
compatibility, kept current) — this doc covers the Components directory's
internal shape, not attack-by-attack usage.

## 1. EquationStores — the equation-holding data structures

All stores share a minimal contract defined in `BaseEqStore` (`EquationStores/BaseEqStore.py:14`):
- `consistent` — whether contradiction-checking is enabled
- `eager` — reduce per-insertion (`True`) vs. defer for batch processing (`False`)
- `filtering` — discard redundant info as equations arrive (progressively "more solved") vs. accumulate everything as a bag
- `insert_equation` / `queue_equation` / `process_pending` / `is_determined`

A second shared layer, `IndexedEqStore` (`EquationStores/IndexedEqStore.py:17`),
factors out monomial-indexing logic shared by `EqStore`, `LUEqStore`, and
`SymbolicEqStore` (but **not** `GroebnerEqStore` — see below):
- **Dynamic mode** (`comb_to_idx=None`): monomials discovered on insertion, backing arrays grow as needed.
- **Static mode** (`comb_to_idx` given): fixed monomial set known ahead of time; accepts raw `ndarray` rows directly; raises `ValueError` on an unseen monomial.
- **Linking** (`link()`, `IndexedEqStore.py:62-77`): one store propagates newly-discovered monomials to others — how RAA's offline phase keeps parallel `g`/`h` stores and their LU rank-trackers in sync (`reduced_algebraic_attack.py:88-93`).

| Store | Base | Representation | Reduces? | Notes |
|---|---|---|---|---|
| `EqStore` | `IndexedEqStore` | dense `ndarray[uint8]` coefficient matrix, one row per equation | No — pure accumulator (`eager=True`, `filtering=False`) | Fed to a batch solver later (Gauss-elim, Gröbner via conversion, etc.) |
| `SymbolicEqStore` | `IndexedEqStore` | `list[BooleanANF]` | No — "the symbolic analog of EqStore" (own docstring) | Same bag semantics as `EqStore`, but keeps `comb_to_idx`/`idx_to_comb` like `EqStore`, so it interoperates with matrix stores via linking/Adapters |
| `LUEqStore` | `IndexedEqStore` | incremental LU decomposition (`upper_matrix`, `lower_matrix`, `solved_for` bit-vector = rank/pivot columns, **not** variable values) | Yes — rejects linearly dependent rows (`filtering=True`, `eager=True`) | Real variable values require back-substitution via `LU_Solver`; `rank`/`is_determined` exposed directly |
| `GroebnerEqStore` | `BaseEqStore` only (**not** `IndexedEqStore`) | `list[BooleanANF]`, no index map — monomials are `frozenset[int]`s directly | Yes — incremental Buchberger algorithm (`eager=False`, `filtering=True`); insertions are queued and only processed via `process_pending`/`_consume_queue` | The genuine odd one out: a live Gröbner basis, not a coefficient-matrix accumulator. This is exactly why Adapters needs a special case for it (§2) |

**Mental model**: `EqStore`/`SymbolicEqStore` are dual representations (matrix vs.
symbolic-list) of the same "no-reduction bag" concept. `LUEqStore` is a
filtering/reducing alternative to `EqStore` for the same indexed-matrix
representation. `GroebnerEqStore` is structurally different — self-indexing,
symbolic-only — and requires monomial enumeration to interoperate with the
other three.

## 2. Adapters — converting between representations

Three files, re-exported via `Adapters/__init__.py`:

- **`equation_repr.py`** — single-equation converters between coefficient
  vector (`ndarray[uint8]`, used by `EqStore`/`LUEqStore`), `BooleanANF` (used
  by `GroebnerEqStore`/`SymbolicEqStore`), and `BooleanFunction` (DAG form,
  used by equation generators): `extract_monomials`,
  `boolean_function_to_coef_vector`, `coef_vector_to_anf`,
  `coef_vector_to_boolean_function`.
- **`store_repr.py`** — store-level converters, the `to_anf_list`/`to_coef_matrix`
  functions used throughout `EquationSolving`. Both duck-type dispatch on
  store shape: LU path (`upper_matrix`+`solved_for`, pulls pivot rows only),
  `EqStore` path (`.equations` is an `ndarray` + `.comb_to_idx`, direct
  slice), Symbolic/Groebner path (`.equations` is a list — `SymbolicEqStore`
  reuses its existing index maps; `GroebnerEqStore` builds **fresh** index
  maps by enumerating every monomial across every stored equation, since it
  has none to begin with). Raises `TypeError` if a store matches none of
  these shapes.
- **`online_insertion.py`** — `make_online_inserter(...)` is a metaprogramming
  factory used by FAA/RAA's online phases: it inspects a store's shape
  (`isinstance(store, IndexedEqStore)`, `store.eager and store.filtering`)
  **once**, then returns specialized closures with no per-iteration type
  checks in the hot loop.

**Rule of thumb**: every conversion in `store_repr.py` is described in its
callers as strictly more expensive than working with a store's native
representation (see §5) — reach for `to_anf_list`/`to_coef_matrix` only when
a store doesn't match what you need natively.

## 3. Annihilators — a mathematical input, not a pipeline stage

An **annihilator** here is a low-degree function `g` (or more generally a
low-degree pair `(g, h)` with `h = f·g`) for a given Boolean function `f`.
This `(annihilator, multiple)` pair is exactly what RAA/FAA's offline phases
consume — see [Algebraic Attacks](../theory/Algebraic%20Attacks.md) for the terminology.

Two independent implementations, both exposing
`annihilators(input_fn, ...) -> (degrees, list[BooleanFunction])`:

- **`SparseAnnihilator.py`** — self-contained symbolic algorithm (Möller-style,
  operating directly on `BooleanANF` monomials), **no dependency on any
  EquationStore**. This is the one actually used throughout `experiments/`
  and cited in [Algebraic Attacks](../theory/Algebraic%20Attacks.md).
- **`GaussianAnnihilator.py`** — a linear-algebra approach that **does**
  depend on `EquationStores`/`EquationSolving`: it inserts degree-bounded
  candidate monomials into an `LUEqStore` (dependence check) and two
  `EqStore`s (coefficient matrices), then calls `GaussElim.reduce_matrix`
  (GF(2) RREF) on the transposed constraint matrices.

**Consumption is manual, not wired**: nothing in the library automatically
feeds annihilator output into a generator or attack. The caller invokes
`SparseAnnihilator.annihilators(f)` (or the Gaussian version) directly,
obtains `(g, h)`, and passes `g` as `annihilator` / `h` as `multiple` into
`RAA_offline`/`FAA_offline`. Every non-experiment call site does exactly this.

## 4. EquationGenerators — producing equations from register clocking

Three generators, all yielding one equation (or a list, if `output_fn` is a
list — normalized by the shared `_utils.normalize_output_fn`) per clock cycle:

| Generator | Yields | Requires | Notes |
|---|---|---|---|
| `CubeEqGenerator` | `ndarray[uint8]` coefficient vectors directly | a precomputed `var_map` (needs a `MonomialProfile` + `variable_blocks` known ahead of time) | numba-JIT "cube sum" algorithm; 2-3 orders of magnitude faster than the alternatives, but only usable when the monomial layout is known statically |
| `SubstitutionEqGenerator` | `BooleanFunction` (ANF form) | nothing precomputed — discovers monomials dynamically | composes `current_equations.compose(feedback_functions)`; only needs per-bit equations for bits actually used by `output_fn`, so it can skip irrelevant state bits |
| `SymbolicEqGenerator` | `BooleanFunction` | nothing precomputed | composes the opposite direction (`feedback_functions.compose(current_equations)`); uniquely accepts an `initialization` param to seed symbolic starting expressions (e.g. partial-key-known scenarios); needs full per-bit expressions every step, so it's the slowest when only some bits matter |

**No generator is hard-wired to a specific store.** The caller (an attack's
offline phase) decides: NAA feeds `CubeEqGenerator`/`SubstitutionEqGenerator`
output into an `LUEqStore`; FAA/RAA feed into an `EqStore` (+ a linked
`LUEqStore` for dynamic rank-tracking). Both `EqStore.insert_equation` and
`LUEqStore.insert_equation` accept either an `ndarray` or a `BooleanFunction`,
so both generator output types are store-compatible — the real constraint is
that `CubeEqGenerator` *requires* a static `var_map`, which in practice forces
a static-mode store when it's used.

## 5. EquationSolving — native store per solver, everything else via Adapters

See `EquationSolving/__init__.py` for the function + thin-wrapper-class
convention itself (module-level `solve()`, a class that stores config and
forwards to it, so attacks can swap solvers polymorphically;
`GuessSolver.guess_and_solve` is the sole bare-function exception). This
section covers which store each solver actually wants.

| Solver | Native store (fast path) | Fallback (works, but rebuilds from scratch) |
|---|---|---|
| `GaussElim` | `EqStore`-shaped (`.equations` is `ndarray` + `.comb_to_idx`) | `to_coef_matrix(store)` — "expensive but valid". `solve()` also has a third, `LUEqStore`-specific branch used only when `additional_constants` is given *and* the store has `upper_matrix`+`solved_for`: it forward-solves the constants through `L` to reconstruct the augmented system. |
| `Grob_Solver` | `GroebnerEqStore` only (returned unchanged) | `to_anf_list(store)`, re-inserted one-by-one into a **fresh** `GroebnerEqStore` and reduced from scratch |
| `LU_Solver` | `LUEqStore`-shaped (`upper_matrix`+`lower_matrix`) | `to_coef_matrix(store)`, then rebuild a brand-new `LUEqStore` by re-inserting every row — "expensive but valid" |
| `SplitGrob_Solver` | `GroebnerEqStore` only | `to_anf_list(store)`, re-`enqueue_equation`'d into a fresh `GroebnerEqStore` |
| `GuessSolver` | N/A — takes an already-extracted `base_solution`/`effect_vectors`, no store awareness at all | N/A |

**Pattern**: each solver has exactly one native store type, and reaches every
other store type through the Adapters conversion functions — always at a
described cost premium (rebuilding a store from scratch instead of reusing
one already built).

### The one hard compatibility restriction

`naive_algebraic_attack.py` explicitly rejects `GrobnerSolver`:

```python
if isinstance(solver, GrobnerSolver):
    raise ValueError(
        "NAA with GrobnerSolver is not supported: NAA's offline phase produces "
        "LU matrices with separate constants, which Gröbner solving cannot use. "
        "Use LUSolver or GaussElimSolver instead."
    )
```

**Why**: NAA's offline phase defers keystream constants entirely — it stores
only the LU-decomposed coefficient matrix and reconstructs the RHS constants
(`additional_constants`) later, at online time, from the keystream (see
[Algebraic Attacks](../theory/Algebraic%20Attacks.md)). `GrobnerSolver` has no mechanism to accept a
separately-supplied constant vector layered onto an already-decomposed
system — Gröbner reduction needs the complete equations (LHS=RHS folded
together) up front. `LUSolver`/`GaussElimSolver` both expose an
`additional_constants` parameter specifically for this.

**FAA and RAA have no equivalent restriction.** Both default to `LUSolver()`
and hard-code `LUEqStore` as their online store — not because anything else
is incompatible, but because their online-phase equations are fully combined
and consistent by construction (no deferred-constants problem exists for
them). A caller *could* pass `GrobnerSolver`/`SplitGrobnerSolver` with a
`GroebnerEqStore` as FAA/RAA's `online_store` and it should work (via
`to_anf_list` conversion where needed), just without the LU speed benefit.
This is an absence of a restriction, not a tested guarantee — nothing
exercises this combination.
