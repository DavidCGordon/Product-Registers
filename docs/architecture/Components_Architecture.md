# Cryptanalysis/Components Architecture

`src/PyPR/Cryptanalysis/Components/` has five subdirectories: `EquationStores/`,
`Adapters/`, `Annihilators/`, `EquationGenerators/` and `EquationSolving/`. This
document describes what each one provides and how they interoperate: which
representations feed into which, which store each solver works on directly and
which it reaches only through a conversion, and the one hard compatibility
restriction between an attack and a solver.

The attacks themselves are referred to by their abbreviations throughout: NAA
(naive algebraic attack), RAA (reduced algebraic attack) and FAA (fast
algebraic attack). Their offline/online phase structure is described in
[Algebraic Attacks](../theory/Algebraic%20Attacks.md), and the stores and
solvers each attack accepts are listed in
[Attack Compatibilities](Attack_Compatibilities.md). This document covers the
components, not attack-by-attack usage.

## 1. EquationStores — the equation-holding data structures

All stores share a minimal contract defined in `BaseEqStore` (`EquationStores/BaseEqStore.py`):
- `consistent` — whether contradiction-checking is enabled
- `eager` — reduce per-insertion (`True`) vs. defer for batch processing (`False`)
- `filtering` — discard redundant info as equations arrive (progressively "more solved") vs. accumulate everything as a bag
- `insert_equation` / `queue_equation` / `process_pending` / `is_determined`

`filtering` is declared on the class, not assigned per instance: it is `True`
exactly for subclasses of `FilteringEqStore`
(`EquationStores/FilteringEqStore.py`), so `store.filtering` and
`isinstance(store, FilteringEqStore)` are the same question and cannot drift
apart. That class is the second axis of the hierarchy, crossing `IndexedEqStore`
rather than nesting under it — `LUEqStore` is both, the Gröbner stores are only
filtering, `EqStore`/`SymbolicEqStore` are only indexed.

What it adds is `num_determined`: how much of the system the store has pinned
down. The two families measure it by different means — an LU store counts pivot
columns (its rank), a Gröbner store counts variables reduction has driven to a
constant — and those are **not** the same notion of "solved": rank counts
constrained dimensions, which still need back-substitution before any variable
has a known value. The quantities are therefore not comparable between stores.
What they share is the part consumers rely on: each rises only when its store
has genuinely learned something, and each reaches `num_vars` exactly when
`is_determined` becomes `True`.

That invariant is why `num_vars` has to mean "variables ever seen". The Gröbner
stores maintain it in `enqueue_equation`, where `unknown_vars` grows; solved
variables are composed out before `idxs_used()` runs, so they cannot be counted
twice, and every other mutation site moves a variable from `unknown_vars` into
`solved_vars`, which leaves the total alone.

`num_determined` is a property of the class, not something each consumer works
out by inspecting the store. The progress figure is read
by the online inserter's stall detector (`Adapters/online_insertion.py`) and by
the summary line that RAA, FAA and the cube attack log after substitution. If
each of those sites tested for store-specific
attributes itself (`hasattr(store, 'solved_vars')`, which reaches past the
interface into Gröbner-specific storage), a new filtering store would have to
be added at every one of them. With `num_determined` on the class, each store
answers for itself.

A second shared layer, `IndexedEqStore` (`EquationStores/IndexedEqStore.py`),
factors out monomial-indexing logic shared by `EqStore`, `LUEqStore`, and
`SymbolicEqStore` (but **not** `GroebnerEqStore` — see below):
- **Dynamic mode** (`comb_to_idx=None`): monomials discovered on insertion, backing arrays grow as needed.
- **Static mode** (`comb_to_idx` given): fixed monomial set known ahead of time; accepts raw `ndarray` rows directly; raises `ValueError` on an unseen monomial.
- **Linking** (`IndexedEqStore.link()`): one store propagates newly-discovered monomials to others — how RAA's offline phase keeps parallel `g`/`h` stores and their LU rank-trackers in sync (`RAA_offline` in `reduced_algebraic_attack.py`).

| Store | Base | Representation | Reduces? | Notes |
|---|---|---|---|---|
| `EqStore` | `IndexedEqStore` | dense `ndarray[uint8]` coefficient matrix, one row per equation | No — pure accumulator (`eager=True`, `filtering=False`) | Fed to a batch solver later (Gauss-elim, Gröbner via conversion, etc.) |
| `SymbolicEqStore` | `IndexedEqStore` | `list[BooleanANF]` | No — pure accumulator, like `EqStore` | Same bag semantics as `EqStore`, but keeps `comb_to_idx`/`idx_to_comb` like `EqStore`, so it interoperates with matrix stores via linking/Adapters |
| `LUEqStore` | `IndexedEqStore` + `FilteringEqStore` | incremental LU decomposition (`upper_matrix`, `lower_matrix`, `solved_for` bit-vector = rank/pivot columns, **not** variable values) | Yes — rejects linearly dependent rows (`filtering=True`, `eager=True`) | Real variable values require back-substitution via `LU_Solver`; `rank`/`is_determined` exposed directly |
| `GroebnerEqStore` | `FilteringEqStore` only (**not** `IndexedEqStore`) | `list[BooleanANF]`, no index map — monomials are `frozenset[int]`s directly | Yes — incremental Buchberger algorithm (`eager=False`, `filtering=True`); insertions are queued and only processed via `process_pending`/`_consume_queue` | Unlike the other three: a live Gröbner basis, not a coefficient-matrix accumulator. This is why Adapters needs a special case for it (§2) |

The four stores relate as follows. `EqStore` and `SymbolicEqStore` hold the same kind of object, an
unreduced bag of equations, in two representations (a coefficient matrix and a
list of ANFs). `LUEqStore` is a filtering/reducing alternative to `EqStore` for
the same indexed-matrix representation. `GroebnerEqStore` is structurally different — self-indexing,
symbolic-only — and requires monomial enumeration to interoperate with the
other three.

## 2. Adapters — converting between representations

The adapters are three modules, re-exported via `Adapters/__init__.py`:

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
- **`online_insertion.py`** — `make_online_inserter(...)` is a
  factory used by FAA/RAA's online phases: it inspects a store's shape
  (`isinstance(store, IndexedEqStore)`, `store.eager and store.filtering`)
  **once**, then returns specialized closures with no per-iteration type
  checks in the hot loop.

A conversion through `store_repr.py` is valid for any store, but a solver that
receives a converted store rebuilds its working store from scratch (§5; the
solvers' docstrings call this "expensive but valid"). Use
`to_anf_list`/`to_coef_matrix` only when a store does not already hold the
representation needed.

## 3. Annihilators — a mathematical input, not a pipeline stage

An **annihilator** of a Boolean function `f` is a nonzero function `g` with
`f·g = 0`. More generally the attacks use a low-degree pair `(g, h)` with
`h = f·g`; the annihilator is the case `h = 0`.
This `(annihilator, multiple)` pair is what RAA/FAA's offline phases
consume — see [Algebraic Attacks](../theory/Algebraic%20Attacks.md) for the terminology.

There are two independent implementations, both exposing
`annihilators(input_fn, ...) -> (degrees, list[BooleanFunction])`:

- **`SparseAnnihilator.py`** — self-contained symbolic algorithm (Möller-style,
  operating directly on `BooleanANF` monomials), **no dependency on any
  EquationStore**. This is the implementation the attack tests use
  (`test/Cryptanalysis/test_attacks.py`).
- **`GaussianAnnihilator.py`** — a linear-algebra approach that **does**
  depend on `EquationStores`/`EquationSolving`: it inserts degree-bounded
  candidate monomials into an `LUEqStore` (dependence check) and two
  `EqStore`s (coefficient matrices), then calls `GaussElim.reduce_matrix`
  (GF(2) RREF) on the transposed constraint matrices.

Annihilators are supplied by the caller: nothing in the library automatically
feeds annihilator output into a generator or attack. The caller invokes
`SparseAnnihilator.annihilators(f)` (or the Gaussian version) directly,
obtains `(g, h)`, and passes `g` as `annihilator` / `h` as `multiple` into
`RAA_offline`/`FAA_offline`.

## 4. EquationGenerators — producing equations from register clocking

There are three generators. Each yields one equation per clock cycle, or a list
of equations if `output_fn` is a list (normalized by the shared
`_utils.normalize_output_fn`):

| Generator | Yields | Requires | Notes |
|---|---|---|---|
| `CubeEqGenerator` | `ndarray[uint8]` coefficient vectors directly | a precomputed `var_map` (needs a `MonomialProfile` + `variable_blocks` known ahead of time) | numba-JIT "cube sum" algorithm; only usable when the monomial layout is known statically. [Cube Equation Generation](../theory/Cube%20Equation%20Generation.md) reports measured comparisons with symbolic composition |
| `SubstitutionEqGenerator` | `BooleanFunction` (ANF form) | nothing precomputed — discovers monomials dynamically | composes `current_equations.compose(feedback_functions)`; only needs per-bit equations for bits actually used by `output_fn`, so it can skip irrelevant state bits |
| `SymbolicEqGenerator` | `BooleanFunction` | nothing precomputed | composes the opposite direction (`feedback_functions.compose(current_equations)`); uniquely accepts an `initialization` param to seed symbolic starting expressions (e.g. partial-key-known scenarios); needs full per-bit expressions every step, so it's the slowest when only some bits matter |

No generator is hard-wired to a specific store. The caller (an attack's
offline phase) decides: NAA feeds `CubeEqGenerator`/`SubstitutionEqGenerator`
output into an `LUEqStore`; FAA/RAA feed into an `EqStore` (+ a linked
`LUEqStore` for dynamic rank-tracking). Both `EqStore.insert_equation` and
`LUEqStore.insert_equation` accept either an `ndarray` or a `BooleanFunction`,
so both generator output types are store-compatible — the real constraint is
that `CubeEqGenerator` *requires* a static `var_map`, which in practice forces
a static-mode store when it's used.

## 5. EquationSolving — native store per solver, everything else via Adapters

Each solver is a module with a module-level `solve()` that does the work and a
thin wrapper class that stores configuration and forwards to it, so attacks can
swap solvers polymorphically; `GuessSolver.guess_and_solve` is the sole
bare-function exception. The package docstring in `EquationSolving/__init__.py`
describes that convention. The table below names each solver by its module
(`Grob_Solver`, `LU_Solver`, `SplitGrob_Solver`); the wrapper classes that
attacks receive as `solver=` are `GrobnerSolver`, `LUSolver`,
`SplitGrobnerSolver` and `GaussElimSolver` (in `GaussElim`). This section
covers which store each solver works on directly.

| Solver | Native store (fast path) | Fallback (works, but rebuilds from scratch) |
|---|---|---|
| `GaussElim` | `EqStore`-shaped (`.equations` is `ndarray` + `.comb_to_idx`) | `to_coef_matrix(store)` — "expensive but valid". `solve()` also has a third, `LUEqStore`-specific branch used only when `additional_constants` is given *and* the store has `upper_matrix`+`solved_for`: it forward-solves the constants through `L` to reconstruct the augmented system. |
| `Grob_Solver` | `GroebnerEqStore` only (returned unchanged) | `to_anf_list(store)`, re-inserted one-by-one into a **fresh** `GroebnerEqStore` and reduced from scratch |
| `LU_Solver` | `LUEqStore`-shaped (`upper_matrix`+`lower_matrix`) | `to_coef_matrix(store)`, then rebuild a brand-new `LUEqStore` by re-inserting every row — "expensive but valid" |
| `SplitGrob_Solver` | `GroebnerEqStore` only | `to_anf_list(store)`, re-`enqueue_equation`'d into a fresh `GroebnerEqStore` |
| `GuessSolver` | N/A — takes an already-extracted `base_solution`/`effect_vectors`, no store awareness at all | N/A |

Every `solve()` takes `verify` (default None), forwarded to `guess_and_solve`:
a function deciding whether a candidate initial state is correct, used in place
of the default check that its keystream matches the one passed in. It exists for
attacks that reach the target only through an interface — the cube attack
passes the `test_fn` of its `access_fns` emulation, so the target's keystream
never reaches the solver. NAA, RAA and FAA leave it unset.

As the table shows, each store-aware solver has one native store type and
reaches other store types through the Adapters conversion functions, which
means rebuilding a store from scratch instead of reusing one already built.
The exceptions are `GaussElim`'s `LUEqStore` branch for `additional_constants`
and `GuessSolver`, which takes no store.

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

The reason is that NAA's offline phase defers keystream constants entirely — it stores
only the LU-decomposed coefficient matrix and reconstructs the RHS constants
(`additional_constants`) later, at online time, from the keystream (see
[Algebraic Attacks](../theory/Algebraic%20Attacks.md)). `GrobnerSolver` has no mechanism to accept a
separately-supplied constant vector layered onto an already-decomposed
system — Gröbner reduction needs the complete equations (LHS=RHS folded
together) up front. `LUSolver`/`GaussElimSolver` both expose an
`additional_constants` parameter specifically for this.

FAA and RAA have no equivalent restriction. Their online-phase equations are
fully combined and consistent by construction, so no deferred-constants
problem exists for them and no solver needs to be excluded. Both default to
`LUSolver()` with an `LUEqStore` as the online store, and the caller can
override either. A caller can pass `GrobnerSolver`/`SplitGrobnerSolver` with a
`GroebnerEqStore` as FAA/RAA's `online_store`, and it works (via
`to_anf_list` conversion where needed).
`test/Cryptanalysis/test_attacks.py` runs exactly that combination, with both
`GrobnerSolver` and `SplitGrobnerSolver`, and requires the full round trip to
recover the secret state.
