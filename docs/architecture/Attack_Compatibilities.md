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
- A `GroebnerEqStore` + `GrobnerSolver`/`SplitGrobnerSolver` combination also works: RAA's online equations are fully combined and consistent by construction, so there is no deferred-constants problem like NAA's. Both pairings recover the secret state in `test/Cryptanalysis/test_attacks.py`, alongside the `LUSolver` and `GaussElimSolver` defaults.

## FAA (Fast Algebraic Attack)

- Offline: annihilator equations go into an `EqStore`; the "reduction" step is Berlekamp-Massey applied to the keystream, producing a `linear_relation` — not a second equation store the way RAA has `multiple_eqs`.
- Online: same defaults as RAA — `LUEqStore(consistent=True)` + `LUSolver`, both overridable.
- The same `GroebnerEqStore` + `GrobnerSolver`/`SplitGrobnerSolver` pairings work here as for RAA, and are covered by the same tests.
- The pair FAA is built for is not the pair RAA wants. FAA's cost is set by deg g, because the linear relation of the sequence h(s_t) cancels h whatever its degree, so a low-degree g with h = f·g nonzero is the intended input. A true annihilator (h = 0) is the degenerate case: there is no linear relation left to exploit. The tests run FAA on a nondegenerate pair (g = x0, h = x0x1 + x0x2 for f = x0x1 + x2) across every online combination, and on the true-annihilator pair separately, so that FAA's own path is exercised and not only the degenerate one.

## Cube Attack

- Offline (`cmpr_cube_attack_offline`): cube candidates come from the output's monomial profile, each with a superpoly degree bound, and are tried lowest degree first (`max_degree` caps them). Each superpoly is recovered exactly up to its bound by cube sums on a simulated register, with the known bits (`known_bits`, the public IV) at their values. A dynamic `LUEqStore` over the linearized monomials is the rank tracker: it keeps an equation only if it is independent of those kept so far, and a candidate is skipped once every monomial it could produce is a pivot. Like RAA's dynamic-path LU stores, it is never solved against. The output is a list of `(cube, t, monomials, constant)` plus a `comb to idx map` holding the constant, every state bit, and each monomial seen.
- Online (`cube_attack_online`): each cube is summed once on the target through `access_fn`, and every equation becomes a coefficient vector over that map. The known bits go in first as the equations x_i = v_i, so the system pins the whole state. Default `online_store` is `LUEqStore(consistent=True)` and default `solver` is `LUSolver`; the `GaussElimSolver` and `GroebnerEqStore` + `GrobnerSolver`/`SplitGrobnerSolver` pairings RAA and FAA accept work too. Groebner is the natural fit when the superpolys are nonlinear, since LU treats each monomial as an independent column.
- The online phase reaches the target only through the `access_fns` interface, which emulates an attacker's view of a real cipher: `access_fn` (run with chosen IV bits) and `test_fn` (does a candidate state reproduce the target's keystream). The solver verifies candidates with `test_fn` through its `verify` hook (see General below), so the target's keystream and initialization rounds stay inside the interface, and a test cannot recover a secret it was handed.
- The IV values must be the same in both phases: the offline superpolys are computed at them. `known_bits` is stored in the attack data and the online phase reads it from there.
- Covered end to end in `test/Cryptanalysis/test_cube_attacks.py`, on a hand-designed 17-bit CMPR, with linear cubes across every pairing and with nonlinear cubes (degree up to 5) across every pairing.

## General

- Every attack accepts `solver=` (NAA) or `solver=`/`online_store=` (RAA, FAA, cube) as override parameters. There is one enforced compatibility check (NAA rejects `GrobnerSolver`), and it is exercised directly.
- Every NAA, RAA and FAA combination in this file is covered end to end in `test/Cryptanalysis/test_attacks.py`: each runs a full offline/online round trip against a small chained CMPR and must return the secret state it was seeded from. Each attack is run on both offline paths — the monomial-profile path (`CubeEqGenerator`, variable count known up front, no rank tracker) and the dynamic path (`SubstitutionEqGenerator`, rank tracked by `LUEqStore`s) — since they build `attack_data` differently. A keystream for RAA or FAA must be `uint8` — its online phase XORs the keystream into a `uint8` coefficient vector in place, which numpy refuses from a wider dtype.
- Every solver's `solve()` takes `verify` (default None) and forwards it to guess-and-prune: a function deciding whether a candidate initial state is correct, in place of comparing its keystream with the one passed in. NAA, RAA and FAA leave it unset; the cube attack passes `test_fn` and no keystream.
- Guess-and-prune (`GuessSolver.guess_and_solve`) handles remaining free variables after reduction. It clocks the register compiled when the feedback function has been compiled and uncompiled otherwise; the dynamic offline path never compiles it, so assuming compiled made NAA and RAA fail there. See [Components Architecture](Components_Architecture.md) §5.
- For equation generation strategy (cube vs composition), see [theory/Cube Equation Generation](../theory/Cube%20Equation%20Generation.md).
