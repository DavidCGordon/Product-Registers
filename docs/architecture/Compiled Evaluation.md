# Compiled Evaluation

This document describes how a `FeedbackFunction` becomes machine code, how
`FeedbackRegister` clocks with it, what each path costs, and why the compiled
loops are written the way they are. That last part includes a design
alternative, inlining the update into the loops, which was measured and
rejected. Source: `FeedbackFunction.generate_python` / `compile` in
`src/PyPR/FeedbackFunctions/FeedbackFunction.py`, and the loops at the end of
`src/PyPR/FeedbackRegister.py`.

All timings here were measured 2026-10-01 on one Windows 11 development machine
and are minima over repeated, interleaved runs. That machine's run-to-run noise
is about ±7% at these durations, and its speed drifted by a factor of two or
more between sessions, so compare figures within a table, not across tables.

## 1. From `fn_list` to compiled functions

`FeedbackFunction.compile()` wraps the statements of
`generate_python(optimize=True)` in two functions and hands each to
`numba.njit`:

- `_compiled(state)` returns the next state as a new array;
- `_compiled_inplace(state, buffer)` writes the next state into `buffer`.

The state is one `uint8` per bit. `generate_python` and its `optimize` flag are
described in [Logic Optimization](Logic%20Optimization.md) §6: the optimized
form is used when it saves operations, and otherwise each bit is written as it
stands in `fn_list`. numba compiles lazily, so the cost of compiling is paid on
the first call of each function, not in `compile()`.

`BooleanFunction.compile()` does the same for a single function, giving
`_compiled(array)`.

## 2. How a register clocks

`FeedbackRegister` keeps two buffers, `_state` and `_prev_state`, and swaps
them on every clock so the update always reads one and writes the other. The
generated code writes bits before it has finished reading, so the two must be
distinct arrays; the swap guarantees it.

| Method | What runs | Cost per clock |
|--------|-----------|----------------|
| `clock()` | A Python call into `_compiled_inplace` | Dominated by the Python call |
| `run(n)` | A numba generator that clocks and yields to Python each cycle | About 2 µs (1.7 to 2.6 measured), almost all of it Python resuming the generator |
| `period()` | A loop compiled as a whole (section 3) | Tens of nanoseconds |

`run` hands control back to Python once per clock, which costs roughly a
hundred times what the update does. Nothing about the update's code changes
that. Work that needs many clocks without Python seeing each state belongs in a
compiled loop of its own, as `period` is.

## 3. The compiled loops

`_period_compiled_unsafe_`, `_period_compiled_safe_` and `_run_compiled_` are
module-level `@numba.njit` functions that take the register's
`_compiled_inplace` as an argument. Each is declared with an explicit
signature in which the update has the numba type
`FunctionType(void(uint8[:], uint8[:]))`. numba therefore compiles each loop
once, when the module is imported, and that one compiled loop serves every
register: the update arrives as a pointer to a function, and each clock is a
call through it.

### The state comparison does not allocate

Both period loops compare two states on every clock. The direct way to write
that, `np.all(a == b)`, builds a temporary boolean array each time. Measured in
isolation that costs 77 ns for 20 bytes and 123 ns for 61 bytes, several times
the update itself. The loops call `_same_state(a, b)` instead, which returns
the same truth value from a loop that stops at the first differing bit: 6 ns
and 17 ns in the same test, where the states differed halfway along. Two
different states of a register usually differ much earlier. Using `_same_state`
in place of `np.all` costs nothing at first use and made every period loop 2.2
to 2.8× faster (the "Original" and "Comparison only" columns below).

### Why the update is not inlined into the loops

The call through a pointer is not free: an update that does nothing measured
18.9 ns per clock in such a loop, against 23.8 ns for a whole 31-bit MPR. It
can be removed without writing any loop twice. Declared with a bare
`@numba.njit` and no signature, a loop has its argument types inferred at the
first call, and a compiled function passed as an argument has a type of its
own. numba then compiles the loop again for each update function it is given,
the call inside it is to a known function, and LLVM can inline it.

That design was implemented and measured, and the library does not use it. The
table gives nanoseconds per clock in `period` for three versions of the loops:
with `np.all` and the signatures ("Original"), with `_same_state` and the
signatures (the current code), and with `_same_state` and the signatures
removed:

| Register | Loop | Original | Comparison only (current) | Comparison and inlining | Inlining alone |
|----------|------|---------:|--------------------------:|------------------------:|---------------:|
| Fibonacci 20 | unsafe | 117.0 | 41.6 | 9.5 | 4.4× |
| MPR 19 | unsafe | 118.5 | 44.6 | 10.4 | 4.3× |
| TFunction(20) | unsafe | 103.7 | 39.4 | 35.9 | 1.10× |
| Fibonacci 20 | safe | 128.4 | 58.6 | 49.0 | 1.20× |
| MPR 19 | safe | 130.5 | 48.8 | 42.9 | 1.14× |
| TFunction(20) | safe | 118.5 | 51.0 | 46.3 | 1.10× |

For the safe rows the figure is time divided by twice the period, which
understates the clocks actually run (see the end of this section), so read
those rows as ratios.

The reasons for rejecting it:

- **It does not happen for larger updates.** The gain is the call, and LLVM
  inlines a small body but not a large one. The 20-bit LFSR and MPR gained
  about 32 ns per clock in the unsafe loop; TFunction(20), whose update is
  larger, gained 3.5 ns. A register big enough for its update to dominate the
  clock is exactly the one that gets nothing.
- **It is less than the headline figure where it counts.** A bare loop with
  nothing but the update showed 6×. In the safe loop, which is the default, it
  was 1.1 to 1.2×: forcing numba to inline (`inline="always"`) measured the
  same, so the call is not what remains there.
- **It competes with compile time.** Each loop is compiled again for every
  feedback function, on first use: 0.5 to 0.9 s more for the safe loop and
  0.16 to 0.28 s for the unsafe loop, across four registers. A register is
  usually asked for its period once, so the saving has to repay that within
  one call. From the table, the safe loop saves at best about 19 ns per state on
  the cycle (the Fibonacci row; the MPR and TFunction rows give about 12 and
  9), which repays 0.5 to 0.9 s at a period of 26 to 47 million; the
  default limit for that loop is $2^{25} \approx 33.5$ million. So within its
  default limit the safe loop roughly breaks even at best. The unsafe loop on
  a small register repays at 5 to 9 million clocks and saves about 2 s at its
  $2^{26}$ limit; that is the one case that wins.
- **`run` pays and gains nothing.** It returns to Python every clock
  (section 2), so the update's call is a hundredth of its cost, yet it would
  still compile a loop per register.
- **It makes the types harder to pin down.** With no signature the loop's
  types are whatever the first caller passes. In the measured implementation
  the limit had to be converted to `np.uint64` at each call site so that a
  plain Python `int` did not change how the arithmetic in the loop was typed,
  and the annotations had to admit that type. A typing mistake would surface
  at the first call for some register instead of when the module is imported.
- **Memory.** Every compiled copy is kept for as long as the loop exists: one
  per loop per compiled feedback function. This is linear and was not the
  deciding factor.

A design between the two is possible, and the library does not implement it:
start in the shared loop and, if the
period has not been found after a few million clocks, continue from that state
in a loop specialized for the register. Only long searches would then pay for
the compilation, and only they would benefit. It would require the unsafe loop
to be resumable.

The safe loop is Brent's cycle detection (`Brent1980`), and it clocks the
register more than once per state. With no preperiod and period $\lambda$, the
first phase tries windows of 1, 2, 4, … clocks and succeeds in the first window
of size $P \ge \lambda$, after $\lambda$ clocks of it. The failed windows before
it total $P - 1$ clocks, and $\lambda \le P < 2\lambda$, so the first phase
takes between $2\lambda - 1$ and $3\lambda$ clocks. The second phase advances
one pointer $\lambda$ clocks. The safe loop therefore runs between $3\lambda$
and $4\lambda$ clocks, against $\lambda$ for the unsafe loop, and compares
states on each of them.

## 4. How close this is to the machine's limit

The compiled loops are not close to that limit, and the gap is in the design,
not the compiler. numba compiles
through LLVM, the backend clang uses, so the same program written in C is
expected to produce much the same machine code; that expectation was not
tested, and no C build was timed. The table compares four ways of
clocking a 31-bit Fibonacci LFSR, each variant checked to give the same output
stream:

| Variant | ns per clock |
|---------|-------------:|
| Update called through a pointer once per clock, as the library's loops do | 46.4 |
| The same byte-per-bit update inlined in the loop | 6.5 |
| The whole state packed into one 64-bit word, parity by shifts | 6.6 |
| Bit-sliced: each array element a 64-bit word carrying 64 registers | 6.3, or 0.098 per register |

The possible improvements are these. None of them is implemented:

- **Bit-slicing.** The generated code is `&`, `^` and `|` on array elements,
  so it works unchanged if each element is a 64-bit word holding one bit from
  each of 64 independent registers. One pass then clocks 64 registers. It
  needs `uint64` state arrays and NOT written as a bitwise complement, where
  `generate_python` now writes `1 - x`. It helps workloads that clock many
  states (cube sums, guess verification) and not `period`, which follows one.
- **Packed state.** One word per register removes the per-bit loads and
  stores. The measurement above used a shift-based parity and came out level
  with the inlined byte version; a population-count instruction would do
  better, but that was not measured.
- **Inlining the update into the loops.** Measured and rejected; see section 3.
- **The cube equation generator's loops.**
  `Cryptanalysis/Components/EquationGenerators/CubeEqGenerator.py` clocks with
  `_compiled`, which allocates a new array per clock where `_compiled_inplace`
  would not. The effect has not been measured.

Cache layout is not a factor at these sizes: a register's state is tens of
bytes and stays in the first-level cache.
