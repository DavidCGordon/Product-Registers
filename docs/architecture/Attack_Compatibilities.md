# Attack Compatibilities (Store x Solver)

See also [Components Architecture](Components_Architecture.md) for the component-level rationale (store taxonomy, solver native/fallback paths, annihilator wiring). This file is the attack-specific compatibility reference; that one is the component-specific one.

For the mathematical foundations of these attacks, see [theory/Algebraic Attacks](../theory/Algebraic%20Attacks.md).

## NAA (Naive Algebraic Attack)

- Offline phase always builds an `LUEqStore` regardless of monomial-profile mode.
- Online phase reconstructs that `LUEqStore` from the stored upper/lower matrices and layers keystream-derived `additional_constants` on top.
- Solver choice is restricted to `LUSolver` (default) or `GaussElimSolver` — both accept `additional_constants` layered onto an LU-shaped store. `GrobnerSolver` is explicitly rejected: Groebner reduction needs the complete equations up front and has no mechanism for separately-supplied constants.

## RAA (Reduced Algebraic Attack)

- Offline: `annihilator_eqs`/`multiple_eqs` are `EqStore`s — static mode (`CubeEqGenerator` + `var_map`) or dynamic mode (`SubstitutionEqGenerator`), the latter linked to a parallel pair of `LUEqStore`s purely for early rank/margin detection. Those LU stores are never solved against; they only decide when to stop generating equations.
- Online: default `online_store` is `LUEqStore`, default `solver` is `LUSolver`. Both are caller-overridable parameters.
- A `GroebnerEqStore` + `GrobnerSolver`/`SplitGrobnerSolver` combination should work here in principle (RAA's online equations are fully combined and consistent by construction — no deferred-constants problem like NAA), but this is untested.

## FAA (Fast Algebraic Attack)

- Offline: annihilator equations go into an `EqStore`; the "reduction" step is Berlekamp-Massey applied to the keystream, producing a `linear_relation` — not a second equation store the way RAA has `multiple_eqs`.
- Online: same defaults as RAA — `LUEqStore(consistent=True)` + `LUSolver`, both overridable.
- Same untested-but-plausible Groebner combination caveat as RAA applies.

## General

- All three attacks accept `solver=` (NAA) or `solver=`/`online_store=` (RAA, FAA) as override parameters. There is one enforced compatibility check (NAA rejects `GrobnerSolver`).
- Guess-and-prune (`GuessSolver.guess_and_solve`) handles remaining free variables after reduction. See [Components Architecture](Components_Architecture.md) §5.
- For equation generation strategy (cube vs composition), see [theory/Cube Equation Generation](../theory/Cube%20Equation%20Generation.md).
