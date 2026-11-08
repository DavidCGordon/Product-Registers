# Attack Compatibilities (Store x Solver)

This document lists, for each attack in `PyPR.Cryptanalysis.Attacks`, which equation stores and solvers it uses by default, which others it accepts, and which combinations the test suite runs end to end. Each attack has an **offline phase** and an **online phase**. The offline phase uses the register and functions derived from the output function (and, for the cube attack, the known public bits), but no keystream. The online phase solves for one keystream. An **equation store** holds the system of equations a phase builds, and a **solver** recovers the initial state from a store.

[Components Architecture](Components_Architecture.md) is the companion reference organized by component: the kinds of store, the store each solver works on directly and the conversions it falls back to, and how annihilators are supplied. The mathematics of the attacks is in [Algebraic Attacks](../theory/Algebraic%20Attacks.md).

## NAA (Naive Algebraic Attack)

The offline phase always builds an `LUEqStore`, regardless of monomial-profile mode. The online phase reconstructs that `LUEqStore` from the stored upper/lower matrices and layers keystream-derived `additional_constants` on top.

`LUSolver` is the default solver and `GaussElimSolver` is the supported alternative: both accept `additional_constants` layered onto an LU-shaped store. `GrobnerSolver` is rejected with a `ValueError`, because Gröbner reduction needs the complete equations up front and has no mechanism for separately supplied constants. That is the only solver check NAA enforces.

## RAA (Reduced Algebraic Attack)

In the offline phase, `annihilator_eqs` and `multiple_eqs` are `EqStore`s, in static mode (`CubeEqGenerator` + `var_map`) or dynamic mode (`SubstitutionEqGenerator`). In dynamic mode they are linked to a parallel pair of `LUEqStore`s used purely for early rank/margin detection. Those LU stores are never solved against; they only decide when to stop generating equations.

In the online phase, the default `online_store` is `LUEqStore` and the default `solver` is `LUSolver`. Both are parameters the caller can override. A `GroebnerEqStore` with `GrobnerSolver` or `SplitGrobnerSolver` also works: RAA's online equations are fully combined and consistent by construction, so there is no deferred-constants problem like NAA's. Both pairings recover the secret state in `test/Cryptanalysis/test_attacks.py`, alongside the `LUSolver` and `GaussElimSolver` pairings.

## FAA (Fast Algebraic Attack)

In the offline phase, annihilator equations go into an `EqStore`. The "reduction" step is Berlekamp-Massey applied to the multiple's output sequence h(s_t), generated on a randomly seeded register, producing a `linear_relation`; no keystream is involved, and there is no second equation store corresponding to RAA's `multiple_eqs`.

In the online phase, the defaults are the same as RAA's: `LUEqStore(consistent=True)` and `LUSolver`, both overridable. The same `GroebnerEqStore` with `GrobnerSolver` or `SplitGrobnerSolver` pairings work here as for RAA, and are covered by the same tests.

The pair FAA is built for is not the pair RAA wants. The linear relation of the sequence h(s_t) cancels h from the online equations whatever its degree, so the degree of the online system is deg g, and a low-degree g with h = f·g nonzero is the intended input. The degree of h still matters elsewhere: the length L of the relation is the linear complexity of h(s_t), each online equation spans a window of L + 1 consecutive keystream bits, and finding the relation takes about 2L terms of h(s_t) offline, so h sets the data requirement and the offline cost. A true annihilator (h = 0) is the degenerate case: there is no linear relation left to exploit. The tests run FAA on a nondegenerate pair (g = x0, h = x0x1 + x0x2 for f = x0x1 + x2) across every online combination, and on the true-annihilator pair separately.

## Cube Attack

In the offline phase (`cmpr_cube_attack_offline`), cube candidates come from the output's monomial profile, each with a superpoly degree bound, and are tried lowest degree first (`max_degree` caps them). Each superpoly is recovered exactly up to its bound by cube sums on a simulated register, with the known bits (`known_bits`, the public IV) at their values. A dynamic `LUEqStore` over the linearized monomials is the rank tracker: it keeps an equation only if it is independent of those kept so far, and a candidate is skipped once every monomial it could produce is a pivot. Like RAA's dynamic-path LU stores, it is never solved against. The output is a `CubeAttackData`: `equations`, a list of `(cube, t, monomials, constant)`, plus a `comb_to_idx` map holding the constant, every state bit, and each monomial seen.

In the online phase (`cube_attack_online`), each cube is summed once on the target through `access_fn`, and every equation becomes a coefficient vector over that map. The known bits go in first as the equations x_i = v_i, so the system pins the whole state. The default `online_store` is `LUEqStore(consistent=True)` and the default `solver` is `LUSolver`; the `GaussElimSolver` and `GroebnerEqStore` + `GrobnerSolver`/`SplitGrobnerSolver` pairings RAA and FAA accept work too. The pairings differ when the superpolys are nonlinear: LU linearizes, treating each monomial as an independent unknown, so it can leave variables undetermined that the polynomial system determines.

The online phase reaches the target only through the `access_fns` interface, which emulates an attacker's view of a real cipher: `access_fn` (run with chosen IV bits) and `test_fn` (does a candidate state reproduce the target's keystream). The solver verifies candidates with `test_fn` through its `verify` hook (see General below), so the target's keystream and initialization rounds stay inside the interface, and a test cannot recover a secret it was handed.

The IV values must be the same in both phases: the offline superpolys are computed at them. `known_bits` is stored in the attack data and the online phase reads it from there.

The cube attack is covered end to end in `test/Cryptanalysis/test_cube_attacks.py`, on a hand-designed 17-bit CMPR, with linear cubes across every pairing and with nonlinear cubes (degree up to 5) across every pairing.

## General

The points in this section apply across the attacks.

Every attack accepts `solver=` (NAA) or `solver=`/`online_store=` (RAA, FAA, cube) as override parameters. There is one enforced compatibility check (NAA rejects `GrobnerSolver` with a `ValueError`), and `test/Cryptanalysis/test_attacks.py` tests that it is raised.

Every NAA, RAA and FAA combination in this file is covered end to end in `test/Cryptanalysis/test_attacks.py`: each runs a full offline/online round trip against a small chained CMPR and must return the secret state it was seeded from. Each attack is run on both offline paths — the monomial-profile path (`CubeEqGenerator`, variable count known up front, no rank tracker) and the dynamic path (`SubstitutionEqGenerator`, rank tracked by `LUEqStore`s) — since they build `attack_data` differently.

A keystream for RAA or FAA must be `uint8`: the online phase XORs the keystream into a `uint8` coefficient vector in place, which numpy refuses from a wider dtype.

Every solver's `solve()` takes `verify` (default None) and forwards it to guess-and-prune: a function deciding whether a candidate initial state is correct, in place of comparing its keystream with the one passed in. NAA, RAA and FAA leave it unset; the cube attack passes `test_fn` and no keystream.

Guess-and-prune (`GuessSolver.guess_and_solve`) handles remaining free variables after reduction. It clocks the register compiled when the feedback function has been compiled and uncompiled otherwise. Both cases occur: NAA's and RAA's dynamic offline paths do not compile the feedback function, and FAA compiles it on both offline paths. See [Components Architecture](Components_Architecture.md) §5.

The choice between the two equation generation strategies (cube sums and composition) is discussed in [Cube Equation Generation](../theory/Cube%20Equation%20Generation.md).
