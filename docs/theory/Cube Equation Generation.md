# Cube-Based Equation Generation

## What this is for

For algebraic attacks on CMPRs, we need to generate a large system of ANF equations — one equation per observed output bit, expressing that bit as a polynomial in the initial state variables. The naive approach is to compose ANFs symbolically through the feedback function tree. This is correct but very slow for large CMPRs, because the internal ANFs grow to exponential size.

The cube-based approach computes the same result (the coefficients of each monomial in the output ANF) via a completely different route — cube summation over register evaluations — that is far cheaper for CMPRs. The tradeoff is that it does not give you the ANF as a symbolic object; it only gives you the numerical coefficients.

## Background: The Cube Attack Identity

The starting point is this identity from the cube attack literature:

$$p(x_1, \ldots, x_n) = T_I \cdot p_{S(I)}(x_1, \ldots, x_n) + Q(x_1, \ldots, x_n)$$

where:
- $T_I = \prod_{i \in I} x_i$ is the monomial over index set $I$
- $p_{S(I)}$ is the **superpoly** of $p$ with respect to $I$ — a polynomial that does not depend on any variables in $I$
- $Q$ is a remainder containing no monomial divisible by $T_I$

The key consequence: if you sum $p$ over all $2^{|I|}$ assignments to the variables in $I$ (a "cube sum"), every term in $Q$ vanishes (each appears an even number of times), and every term in $T_I \cdot p_{S(I)}$ that has any additional variables in $I$ also vanishes. The only surviving term is $T_I$ itself — multiplied by the constant $p_{S(I)}(0, \ldots, 0)$. So:

$$\sum_{v \in \{0,1\}^{|I|}} p(v, x_{|I|+1}, \ldots, x_n) = p_{S(I)}(0, \ldots, 0)$$

This gives a direct formula: the coefficient of any monomial $T_I$ in $p$ equals the cube sum of $p$ over $I$. Applied to our register system, $p$ is the output function evaluated after $t$ clock cycles, and the variables $x_i$ are the bits of the initial state.

## The Two Key Optimizations

### 1. Precomputed Evaluations

A cube sum for index set $I$ requires summing $p$ over all $2^{|I|}$ assignments to variables in $I$, with the remaining variables fixed. Many different cube sums share the same evaluation of $p$ at the same initial state. The first optimization is to precompute all needed evaluations and reuse them.

In practice: maintain one register for each initial state we will ever need (across all cube sums). On each clock cycle, update all registers in parallel, and evaluate the output function on all of them in parallel, filling a vector of evaluations. All cube sums are computed by summing subsets of this precomputed vector.

This completely separates "running registers" from "summing coefficients."

### 2. Graded Cube Splitting

If cube sums are computed in **graded order** (degree-1 first, then degree-2, etc.), any degree-$d$ sum can be computed from already-known smaller sums, cutting the cost from $O(2^d)$ to $O(2^{d/2})$ per monomial.

Here is the argument in detail. Let:
- $V$ be the set of all binary vectors (initial states)
- $V_I \subseteq V$ be the $2^{|I|}$ initial states used for the cube over $I$
- $F(V') = \bigoplus_{v \in V'} p(v)$ be the XOR-sum of evaluations over a set of states

The key property is that $F$ is linear (over GF(2)), so:

$$\forall A, B \subseteq V: \quad F(A) \oplus F(B) = F(A \cup B) \oplus F(A \cap B)$$

Now pick any subset $I' \subseteq I$. Consider the collection of subcubes $V_{I \setminus \{i\}}$ for each $i \in I'$ (the cubes that "drop" one variable from $I'$ at a time).

Every vector in $V_I$ that does *not* have all bits in $I'$ set to 1 appears in at least one of these subcubes. The vectors that have all bits in $I'$ set form a remaining set $V_{\text{remaining}}$ of size $2^{|I| - |I'|}$. So:

$$F(V_I) = F\!\left(\bigcup_{i \in I'} V_{I \setminus \{i\}}\right) \oplus F(V_{\text{remaining}})$$

The second term costs $2^{|I| - |I'|}$ to compute directly. For the first term, apply inclusion-exclusion with the identity $V_{I \setminus S_1} \cap V_{I \setminus S_2} = V_{I \setminus (S_1 \cup S_2)}$:

$$F\!\left(\bigcup_{i \in I'} V_{I \setminus \{i\}}\right) = \bigoplus_{S \subseteq I',\, S \neq \emptyset} F(V_{I \setminus S})$$

Since we are working in graded order, every $F(V_{I \setminus S})$ on the right-hand side was already computed (it's a smaller cube sum). There are $2^{|I'|} - 1$ such subsets.

**Total cost** = $(2^{|I'|} - 1) + 2^{|I| - |I'|}$

This is minimized by choosing $|I'| \approx |I|/2$. With $d = |I|$:

$$\text{cost per monomial} = 2^{\lfloor d/2 \rfloor} + 2^{\lceil d/2 \rceil} - 1 = O(2^{d/2})$$

Compare to the naive cost of $2^d$ — this is exponentially better in the degree.

### 3. Loop Precomputation and JIT

Because the same index patterns are reused every clock cycle (the cube structure is fixed), all the indices into the evaluation vector can be precomputed once. The entire coefficient computation then reduces to a deeply nested loop of XOR operations on binary arrays — no symbolic monomial objects, no data structure overhead.

This structure is:
- Easily parallelized
- Entirely JIT-compilable (Numba in this library)
- Cache-friendly

The coefficient of the algorithm is dramatically smaller than for ANF composition.

## Asymptotic Cost Estimates

Both methods' costs are measured in **coefficient flips per ANF generated** — this unit is hardware-agnostic and allows a fair comparison of scaling without the implementation advantage of the cube method distorting the picture.

### ANF Composition Cost (Upper-Bounded Favorably)

To make the comparison as favorable as possible for composition, assume:
- Updating a monomial coefficient costs $O(1)$ (requires a hashmap or similar structure)
- It is possible to iterate over only the monomials that could possibly be nonzero (requires metadata about the ANF's monomial profile)

Neither of these is true in the current library implementation, so this analysis *understates* the true composition cost.

Under these assumptions:
- **Multiplying** two ANFs with monomial counts $|M_1|$ and $|M_2|$: cost $(|M_1| \cdot |M_2|)/2$
- **Adding** two ANFs: cost $\min(|M_1|, |M_2|)$

These are evaluated bottom-up on the binarized function tree using dynamic programming. Each node computes its own monomial profile and cost from its children's profiles. The total is the cost at the root, summed over all feedback functions and the output function.

### Cube-Based Cost

As derived above, for each output monomial of degree $d$, the cost is:

$$2^{\lfloor d/2 \rfloor} + 2^{\lceil d/2 \rceil} - 1$$


To get the total cost, iterate over all possible degree vectors (one degree per register block) using the degree rollover logic from `MonomialProfile.get_monomials()` (see [Monomial Profile Theory](Monomial%20Profile%20Theory.md) for how profiles bound the monomial set):
- Maintain a vector of current degrees, one per register.
- Increment one entry at a time; when it exceeds the block's maximum, "roll over" (reset and carry to the next block).
- For each degree vector, multiply the per-monomial cost by the number of monomials with those degrees (binomial coefficients per block, multiplied across independent blocks).

## Practical Comparison

### Asymptotic estimates

The cost estimators in `CostEstimation.py` count *coefficient flips* — a hardware-agnostic unit that strips away the implementation advantages of one method over the other so that scaling behavior is visible. The estimates are deliberately favorable to composition (O(1) coefficient lookup, iterating only over nonzero monomials) so that any advantage for cube is real, not an artifact of assuming a bad composition implementation.

### Measured speedup on CMPRs

For M7+M5+M3 (15 bits, `arman_template`, VAR(0) output):

| Quantity | Value |
|---|---|
| var_map size (states maintained by cube) | 3,046 |
| Estimated comp cost | 405,859 coefficient flips |
| Estimated cube cost | 41,792 coefficient flips |
| Estimated ratio (comp/cube) | 9.7× |
| Cube steady-state (measured) | 3.5 ms/eq |
| Symbolic at saturation (measured) | 69 ms/eq |
| Measured speedup | ≈ 20× |

The measured speedup (20×) exceeds the estimated ratio (9.7×). The difference comes from two compounding factors: (1) the JIT-compiled cube loop processes array XORs far faster than Python's per-object ANF operations, and (2) the composition cost estimate is already favorable, so the true symbolic cost is higher than counted. Both push the measured advantage beyond what the coefficient-flip model predicts.

The asymptotic ratio of ~6,000× reported from an earlier large-scale benchmark (10,000 equations in 5 s vs. 10 equations in 30 s) remains valid at those scales; the 20× figure here is for a smaller register where the symbolic ANFs haven't fully saturated.

### The crossover

A crucial nuance: cube is not *always* faster for CMPRs. In the very early equations (t ≈ 1–15 for a 15-bit register), the symbolic ANFs are still small and Python operations on them run faster than the JIT overhead from managing 3000 register states. The cube method's advantage only appears once the symbolic ANFs saturate at their maximum size.

![M7+M5+M3 crossover: per-equation time vs equation number](../../../test/Tools/figures/cmpr_speed_crossover.png)

*Generated by `generate_cmpr_crossover_plot()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py). Requires ~8 s.*

The plot shows per-equation time (log scale) for both methods as a function of clock cycle. Symbolic starts near 0.5 ms/eq, grows rapidly as the chaining ANFs fill in, and saturates around equation 15 at ~70 ms/eq. Cube is roughly constant at 3.5 ms/eq throughout. The crossover occurs near equation 15, which means for any algebraic attack that needs more than O(register size) equations, cube is the right choice.

### Steady-state comparison: all three generators

The library provides three symbolic generators. Their performance at CMPR steady state (post-ANF-saturation) differs dramatically:

- **Cube**: constant cost per equation regardless of clock cycle. The JIT-compiled state array dominates; cost is proportional to var_map size.
- **Symbolic** (`SymbolicEqGenerator`): composes small feedback functions *forward* through the large state. At saturation the state ANFs are fixed in size, so per-equation cost levels off.
- **Substitution** (`SubstitutionEqGenerator`): tracks only the bits used by the output function, but composes the *large* saturated ANF *backward* through all n feedback functions at every step — significantly more expensive than forward composition at saturation. Note: the Substitution generator's docstring claims it is "generally slightly faster than the symbolic generator", but this claim does not hold for CMPRs at steady state.

| Config | Cube | Symbolic | Substitution |
|---|---|---|---|
| M7+M5+M3 (15b) | 8 ms/eq | 51 ms/eq | 326 ms/eq |
| M7+M5+M3+M2 (17b) | 15 ms/eq | 2,335 ms/eq | >4,000 ms/eq (impractical) |

The estimated ratio for 17b is 37.9× (asymptotic); the measured cube speedup over symbolic is 153×.

![CMPR steady-state: cube vs symbolic vs substitution](../../../test/Tools/figures/cmpr_steady_state.png)

*Generated by `generate_cmpr_steady_state_comparison()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py). Requires ~60 s.*

## When Cube Helps (and When It Doesn't)

The cube method is not universally better. It wins or loses depending on how the internal ANF complexity compares to the output ANF complexity. For the generator taxonomy and how CubeEqGenerator fits alongside SymbolicEqGenerator and SubstitutionEqGenerator, see [Components Architecture](../architecture/Components_Architecture.md) section 4.

**CMPRs:** The cube method helps dramatically, and helps more as size grows. CMPRs have a complex output (the chaining functions make the internal ANFs grow to exponential size over time, as bounded by the [monomial profiles](Monomial%20Profile%20Theory.md)). Composition requires multiplying through all those internal ANFs at every clock cycle. The cube method skips all of that — it precomputes a fixed set of register evaluations once and then performs constant-cost cube sums per equation. The per-equation cost is determined entirely by the output's monomial structure, not the internal complexity.

The plot below sweeps all CMPR configurations whose component sizes appear in `MPR_MAP` (sizes 2, 3, 5, 7, 13, 17, 19, 31, 61), over total register sizes up to 100 bits, using the `fast_template` with `max_and=4`. Each point is one configuration; color encodes total register size. The $y$-axis is $\log_2(\text{comp}/\text{cube})$, so positive means cube is cheaper. The dashed line is a linear fit showing that cube's advantage grows roughly linearly with $\log_2(\text{comp cost})$, i.e., grows polynomially with register size.

![CMPR cost ratio: log2(comp/cube) vs log2(comp cost)](../../../test/Tools/figures/cmpr_cost_ratio.png)

*Generated by `generate_cmpr_ratio_plot()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py).*

**Filter generators** (simple MPR + nonlinear output function): The cube method *hurts* both asymptotically and in practice. The internal state update for an MPR is linear: each state bit at time $t$ is a linear combination (XOR) of the $n$ initial state bits, and this never changes. So the symbolic ANFs for individual state bits stay at degree 1 forever. Composing a degree-$d$ AND output through these linear functions is cheap: it involves multiplying $d$ linear polynomials, with a result that has at most $\binom{n}{d}$ monomials. The cube method, by contrast, must maintain one register state per entry in var\_map -- which for an AND of degree $d$ on $n$ variables grows as $\sum_{k \leq d}\binom{n}{k}$ (sized via the output function's [root expression](Root%20Expressions%20and%20LC%20Estimation.md)), an exponential in $d$. The cost of updating that state array per equation grows linearly with var\_map size.

The asymptotic plot below sweeps each MPR in `MPR_MAP` against AND-output degrees from 1 up to the register size. Points cluster below the zero line — cube is more expensive — and the gap widens as comp cost grows, which is the opposite of the CMPR trend.

![Filter generator cost ratio: log2(comp/cube) vs log2(comp cost)](../../../test/Tools/figures/filter_cost_ratio.png)

*Generated by `generate_filter_ratio_plot()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py).*

The measured benchmark below confirms this for M13 (13-bit MPR) at AND degrees 2–10. Symbolic stays faster at every degree. The estimated comp/cube ratio (green, right axis) crosses 1.0 near degree 10, predicting that cube should *start* to become competitive there. The measured data shows symbolic is still ~1.5× faster even at degree 10, because the JIT overhead from maintaining 8,100 register states is a large constant that the asymptotic argument does not capture. A genuine measured crossover in favor of cube would require a much larger MPR ($n \approx 31+$) at high degree, where the symbolic polynomial products finally become more expensive than the JIT state-array overhead.

![M13 filter generator: measured cube vs symbolic per-equation time](../../../test/Tools/figures/filter_speed_comparison.png)

*Generated by `generate_filter_speed_comparison()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py). Requires ~30 s.*

**T-functions (triangular feedback, binary-counter base):** The cube method's advantage is mixed and size-dependent. A T-function has $n$ single-bit blocks — one per bit — so `var_map` grows exponentially with size (TF(11) → 1025 states). Unlike a CMPR, the output ANF never permanently saturates: per-equation symbolic cost oscillates, peaking at every power-of-2 clock cycle as the carry chain reaches maximum depth. Peaks grow: the $t = 2^k$ peak is larger than $t = 2^{k-1}$. The correct characterization is the average over many equations spanning multiple oscillation cycles. For small T-functions ($n \le 9$), symbolic wins on average because the oscillation peaks are modest relative to the cube state-array overhead. For larger T-functions, the peaks grow enough that the picture is less clear.

**CrossJoin registers (single block, paired AND terms):** The entire register is one monomial-profile block, so `var_map` stays small (~14 states for a 13-bit CrossJoin). Symbolic saturates within the first 20 equations and then costs only a few ms/eq. Cube incurs JIT state-array overhead on top of its intrinsically low per-equation cost, but the var_map is so small that the overhead is not amortized. Symbolic wins throughout.

The plot below shows measured $\log_2(\text{sym}/\text{cube})$ vs register size for all three architecture families. Positive $y$ means cube is faster; negative means symbolic is faster. CMPR grows rapidly with size (symbolic saturates to exponential ANFs); CrossJoin and TFunction stay negative or near zero (small var\_map, non-saturating or cheap symbolic cost).

![Architecture sweep: measured log2(sym/cube) vs size](../../../test/Tools/figures/architecture_sweep.png)

*Generated by `generate_architecture_sweep()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py). Requires ~120 s.*

The per-equation trace below shows these three qualitatively different shapes side by side: CMPR's monotone saturation, CrossJoin's fast flat plateau, and TFunction's periodic oscillation with growing power-of-2 peaks.

![Qualitative trace: per-equation symbolic cost for three architectures](../../../test/Tools/figures/qualitative_trace.png)

*Generated by `generate_qualitative_trace_plot()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py). Requires ~25 s.*

### How chaining degree affects the CMPR advantage

The asymptotic advantage of cube over composition for a CMPR grows with the `max_and` chaining parameter. Higher `max_and` means more variables per AND gate in the chaining template, which drives internal ANFs to larger polynomials and makes composition proportionally more expensive. The cube cost, determined only by the output's monomial profile, increases much more slowly.

![CMPR asymptotic cube advantage vs max_and](../../../test/Tools/figures/cmpr_chaining_sweep.png)

*Generated by `generate_cmpr_chaining_sweep()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py).*

### How output AND degree affects the filter-generator disadvantage

For filter generators (linear MPR + AND output), the degree sweep below shows that composition remains cheaper than cube across all tested degrees and register sizes. Lines stay below $y = 0$ throughout; the gap widens as degree grows.

![Filter degree sweep: asymptotic comp/cube vs AND degree](../../../test/Tools/figures/filter_degree_sweep.png)

*Generated by `generate_filter_degree_sweep()` in [`test/Tools/cost_estimation_test.py`](../../../test/Tools/cost_estimation_test.py).*

**General heuristic:** If the internal ANF growth is the bottleneck (large feedback functions, many components, deep chaining), the cube method wins at steady state. If the internal state update is linear (no chaining) or if you only need a small number of equations (before the symbolic ANFs saturate), composition is better.

## API

```python
from PyPR.Tools.CostEstimation import estimate_cost_comp, estimate_cost_cube

mps = C.monomial_profiles()
comp = estimate_cost_comp(C, output_fn, mps)
cube = estimate_cost_cube(C, output_fn, mps)
```

- `estimate_cost_comp(feedback_fn, output_fn, monomial_profiles)` — estimates composition cost in coefficient flips, walking the binarized function tree bottom-up.
- `estimate_cost_cube(feedback_fn, output_fn, monomial_profiles)` — estimates cube-based cost using the graded-split formula, iterating over all degree vectors via rollover logic.
- Both accept any `FeedbackFunction` (not just CMPR) and any `BooleanFunction` as output.
- `monomial_profiles` is the list returned by `C.monomial_profiles()` — one `MonomialProfile` per register bit.
