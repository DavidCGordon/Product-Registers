# Trajectory LC vs Column LC

These are two different quantities that arise in the analysis of [algebraic attacks](Algebraic%20Attacks.md) on CMPRs. They have an inverse relationship that is counterintuitive and critical for FAA rank predictions.

## Definitions

**Trajectory LC** of a state bit or monomial: the linear complexity of $g(s_t)$ when $g$ is evaluated along a trajectory from some specific initial state. This can be estimated via the bit's [RootExpression](Root%20Expressions%20and%20LC%20Estimation.md) via `re.upper()`.

**Column LC** of an initial-state monomial $T_I$: the linear complexity of the column sequence $c_I(t)$ in the equation matrix — i.e., how fast the coefficient of $T_I$ in $g(s_t)$ oscillates as $t$ increases. This is NOT the same as trajectory LC.

## The Inverse Relationship (Cube-Sum Mechanism)

The [cube-sum identity](Cube%20Equation%20Generation.md) connects the two:

$$c_I(t) = \bigoplus_{v \in \{0,1\}^{|I|}} g(s_t^v)$$

where $s_t^v$ is the trajectory starting from the state where bits in $I$ are set to $v$, all others zero.

This identity creates an **inverse relationship** between trajectory LC and column LC across the CMPR block structure:

- **Upstream/simple bits** (e.g. M7 in M7→M5→M3): small trajectory LC (they are "pure" source bits, period matches their block). BUT their column LC is **large** — the M7 bit propagates through many nonlinear chaining stages to the output, creating highly complex coefficient sequences.

- **Downstream/accumulated bits** (e.g. M3): large trajectory LC (they have accumulated all the upstream nonlinear effects in their steady-state behavior). BUT their column LC is **small** — when seeded from $e_j$ (an M3 initial bit), the cube-sum cancels much of the M5/M7 influence, leaving only the direct (low-degree) contribution of that bit.

## The All-Zeros State

The all-zeros state is **NOT** a fixed point for CMPR — it transitions to a nonzero state immediately (because chaining terms can produce nonzero output even from zero inputs, depending on the chaining structure). This means the cube-sum XORs $g(s_t^{e_j})$ with a nontrivial sequence $g(s_t^0)$, not a constant. The resulting column LC can be smaller than trajectory LC when $g(s_t^{e_j})$ and $g(s_t^0)$ share roots that cancel.

## Practical Implications

**Do not use `bit_re.upper()` or `product_of_bit_res.upper()` to predict column LC.** These overestimate downstream column LCs and underestimate upstream column LCs. When analyzing FAA rank loss or equation matrix structure for CMPR, use empirical column LC data rather than RootExpression upper bounds. The [monomial profile](Monomial%20Profile%20Theory.md) of a bit tells you which monomials *can* appear but not how their coefficient sequences (columns) behave over time.

## Open Question

What is the exact formula for column LC as a function of block sizes, annihilator degree, `max_and` chaining parameter, and monomial profile? The "lower bound" pattern is what we're trying to characterize. See `experiments/scratch/lc_table.txt` for systematic data.

Understanding this is essential for bounding the rank of the FAA equation matrix for CMPR and predicting attack complexity.
