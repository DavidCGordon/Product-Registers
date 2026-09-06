"""GF(2) equation solvers used by the algebraic attacks.

Convention: each solver module (`GaussElim`, `Grob_Solver`, `LU_Solver`,
`SplitGrob_Solver`) defines a module-level `solve()` function that does the
actual work, plus a thin class wrapper (`GaussElimSolver`, `GrobnerSolver`,
`LUSolver`, `SplitGrobnerSolver`) whose `__init__` stores per-instance config
(e.g. `simplify_mode`, `additional_constants`) and whose `.solve()` just
forwards to the function.

The class exists because attack code (`Cryptanalysis/Attacks/*_algebraic_attack.py`)
treats "solver" as a pluggable strategy object: it holds a `solver` instance,
calls `solver.solve(...)` polymorphically, and sometimes does
`isinstance(solver, ...)` checks to reject incompatible combinations.

`GuessSolver.guess_and_solve` is the one exception -- a bare function with no
class wrapper. It's never selected by attack code as an interchangeable
strategy; every other solver's `solve()` calls it directly by name as a
shared terminal step, so no polymorphism is needed.

When adding a new solver: use the function + thin-wrapper-class pattern if
attacks should be able to swap it in; use a bare function if it's an
internal helper always invoked by name.
"""
