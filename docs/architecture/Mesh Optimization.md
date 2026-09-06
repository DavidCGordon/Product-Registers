# Mesh Optimization

The **mesh optimization** is an alternative algorithm for computing root expressions and monomial profiles in CMPRs with well-separated block sizes. Instead of composing root expressions gate-by-gate through the chaining function's ANF (the default path), it directly enumerates all possible weight distributions across component fields — an approach analogous to the multinomial theorem.

The optimization lives in `src/PyPR/Tools/RootCounting/MeshOptimization.py` and is invoked automatically by `CMPR.root_expressions()` and `CMPR.monomial_profiles()` when the eligibility conditions are met. Pass `force_default=True` to bypass it.

## The Multinomial Analogy

Consider a CMPR with $q$ MPR blocks of distinct sizes $s_1, \ldots, s_q$. Each block's roots live in a different extension field ($\mathbb{F}_{2^{s_i}}$), and the cross-field product theorem ([Root Expressions — Proposition 3](../theory/Root%20Expressions%20and%20LC%20Estimation.md#the-five-root-expression-propositions)) guarantees that root choices from distinct fields are independent.

The chaining function at block $i$ has some degree $d_i$ as a polynomial in the preceding block's state. When the chaining is simple (see eligibility below), the root expression for the driven block is determined entirely by how the available weight distributes across the upstream fields. This parallels the multinomial expansion: given root expressions $\text{RE}_1, \ldots, \text{RE}_{i-1}$ from the upstream blocks, the chaining computes

$$\left(\sum_{j < i} \text{RE}_j\right)^{d_i}$$

and the resulting terms correspond to weight distribution vectors $(k_1, \ldots, k_{i-1})$ with $\sum_j k_j \leq d_i$, where $k_j$ factors come from block $j$. By Proposition 2 (same-field product), taking $k_j$ factors from a block of size $s_j$ gives coset weight $\min(s_j, k_j)$. So each distribution vector maps directly to a product term $\prod_j \langle s_j \cdot \min(s_j, k_j) \rangle$ without tracing through individual gates.

The multinomial coefficients are irrelevant over GF(2) — the root expression is an upper-bound object that tracks which terms *can* appear, not how many copies are generated. What matters is which distribution vectors produce non-dominated terms.

### Worked Example

For a CMPR with blocks of sizes $[7, 5, 3]$ and chaining degrees $[1, 3, 3]$ (block 0 is the source, block 2 is driven with degree-3 chaining), the mesh for block 2 enumerates weight distributions across all three fields:

| Weight vector | Root expression term |
|---|---|
| $[0,\; 0,\; 1]$ | $\langle 3 \cdot 1 \rangle$ |
| $[0,\; 3,\; 0]$ | $\langle 5 \cdot 3 \rangle$ |
| $[3,\; 2,\; 0]$ | $\langle 7 \cdot 3 \rangle \times \langle 5 \cdot 2 \rangle$ |
| $[6,\; 1,\; 0]$ | $\langle 7 \cdot 6 \rangle \times \langle 5 \cdot 1 \rangle$ |
| $[9,\; 0,\; 0]$ | $\langle 7 \cdot 7 \rangle$ |

The first row is the driven block's own $\langle 3 \cdot 1 \rangle$ root — the starting point. Each subsequent row represents one way the chaining degree can distribute weight into the upstream fields. Weight in field $\mathbb{F}_{2^7}$ saturates at 7 (the field degree), so $[9, 0, 0]$ maps to $\langle 7 \cdot 7 \rangle$, not $\langle 7 \cdot 9 \rangle$.

## Why Distinct Sizes Are Required

When two blocks share the same size $s$, their roots live in the same field $\mathbb{F}_{2^s}$. In this case:

- **Cross-field independence breaks down.** Proposition 3 requires distinct field degrees. Two blocks of the same size contribute roots to the same axis in the [coset region grid](../conventions/Coset%20Regions.md), so their products are same-field products (Proposition 2) that add coset weights, not independent axes whose counts multiply.
- **The weight distribution is no longer axis-separable.** The multinomial analogy assigns each block its own axis and distributes weight independently along each. With shared sizes, the axes collapse, and the correct result requires tracking how same-field products interact through Jordan partitions (see [Root Multiplicities — Jordan Decomposition and Products](../theory/Root%20Multiplicities%20and%20Jordan%20Decomposition.md#jordan-decomposition-and-products)).

This is the core reason T-functions cannot use the mesh optimization — every block uses the same polynomial $(1+x)$ with size 1.

## The Iterator

The mesh iterator (`_re_mesh_iterator` in `MeshOptimization.py`) enumerates weight distribution vectors using a depth-first tree traversal with integrated pruning. The iterator is numba-jitted for performance.

### Degree Extraction

Before the mesh runs, the CMPR driver extracts the chaining degree for each block. The pipeline is:

1. `MonomialProfile.from_merged()` computes the block's merged chaining function as a monomial profile.
2. `.to_BooleanFunction()` converts to an ANF representation where each ANF term becomes an `AND` of `VAR(block_id)` references, with repeated VARs representing same-field powers (the code comment "Repeated Vars => can't use normal degree" notes this).
3. The degree is `max(len(term.args))` across all ANF terms — the maximum number of factors in any single product.

This degree captures the algebraic degree of the chaining as a function of its upstream block.

### Tree Traversal

The tree is organized from the rightmost block (the driven block) to the leftmost (the source):

1. **Initialization:** Weight 1 is placed at the rightmost block: `arr[0][-1] = 1`. This initial unit represents the driven block's own $\langle s \cdot 1 \rangle$ root — the resolvent contribution that exists regardless of chaining.

2. **Trading:** At each depth, weight from the current block is "traded" to the block to its left at an exchange rate of $d_i : 1$ — one unit of weight at block $i$ produces $d_i$ units at block $i-1$. The exchange rate reflects the chaining degree: a degree-$d$ chaining expands one root into $d$ cross-field factors. The `num_traded` is at least 1 (to make progress) but jumps to `weight - size` when the block is oversaturated, since excess beyond saturation has no effect on the coset class and should be redistributed.

3. **Saturation:** A block with weight $\geq s_i$ (its field degree) is saturated — additional weight does not increase the coset class (Proposition 2 caps at the field degree). The excess is available for further trading leftward.

4. **Accumulated rates:** The `values` array stores $\texttt{values}[i] = d_0 \cdot d_1 \cdots d_i$, the total weight one unit at position $i$ is worth when fully expanded toward the source. This is used by the fillable-space pruning check.

### Pruning

Two conditions prevent the iterator from yielding dominated or duplicate terms:

**Leftmost saturation rule.** When multiple blocks are saturated (weight $\geq$ size), the excess weight can sit in any of them without changing the resulting coset classes (since weight beyond saturation is invisible). To avoid generating the same effective term multiple times, the iterator only yields vectors where all excess is concentrated in the leftmost saturated block.

**Fillable space rule.** If a partially filled block ($0 < \text{weight} < s_i$) could be increased by redistributing excess from oversaturated blocks (converting via the accumulated exchange rates in `values`), then the current vector is strictly dominated by that improved version and is skipped. This check uses `overshoot * values / values` to convert excess weight at each position into equivalent weight at the partially-filled position.

Together, these ensure the iterator yields exactly the set of maximal, non-redundant weight distribution vectors — the same set that the default algorithm would produce after expensive pairwise dominance checks.

## Output Construction

Each yielded weight vector is converted to a `JordanSet` or `TermSet`:

- **Root expressions** (`re_compute_single_mesh`): For each nonzero entry $k_j$ in the vector, the field degree $s_j$ and coset weight $\min(s_j, k_j)$ are recorded. The Jordan multiplicity set is $\{1\}$ — the mesh does not track higher Jordan lengths because the distinct-size precondition eliminates the same-field interactions that produce nontrivial Jordan structure.
- **Monomial profiles** (`mp_compute_single_mesh`): The analogous conversion to `TermSet` objects, recording block IDs and saturated weights.

After iteration, the CMPR driver handles two things the mesh does not:

1. **Constants:** If the chaining function can produce a constant term (`CONST(1)` in the ANF), the driver adds `RootExpression.logical_one()` (or `MonomialProfile.logical_one()`). Constants are detected during degree extraction and stored in the `constants_possible` array.
2. **Locked blocks:** The `locked_list` parameter zeros out sizes for unlocked blocks (`sizes *= locked_list`), so weight assigned to them immediately saturates at 0 and is available for trading — effectively making unlocked blocks invisible in the root expression.

## Eligibility Conditions

The mesh optimization is used when all of the following hold (checked in `CMPR.monomial_profiles()` and `CMPR.root_expressions()`):

| Condition | Reason |
|---|---|
| No 1-bit MPR blocks | Size-1 blocks produce roots in $\mathbb{F}_2$ ($\alpha = 1$ only), which interacts with every other field's $\alpha = 1$ root. The weight distribution model does not account for this shared-root interaction. |
| No repeated block sizes | Distinct fields are required for the cross-field independence (Proposition 3) that makes the multinomial distribution valid. Repeated sizes collapse to same-field products requiring Jordan partition handling. |
| Simple chaining | Each block's chaining function must reference only the immediately preceding block (not arbitrary upstream blocks). This ensures a single degree $d_i$ captures the full chaining structure — complex topologies would require tracking degree per upstream block rather than a single scalar. |

When any condition fails, `CMPR` falls back to the default gate-by-gate composition algorithm (`_re_default` / `_mp_default`), which handles arbitrary field overlaps, Jordan interactions, and chaining topologies at higher computational cost.

## Performance

The default composition algorithm has cost $O(nk^2 N^4)$ where $N \leq \prod s_i / \max(s_i)$ is the maximum number of terms in any root expression (see [Root Expressions — Runtime Complexity](../theory/Root%20Expressions%20and%20LC%20Estimation.md#runtime-complexity-of-root-counting)). The gate-by-gate composition builds many intermediate expression objects and performs expensive dominance-pruning at each step.

The mesh avoids intermediate objects entirely — it generates the final set of maximal terms directly via the numba-jitted iterator. The pruning is integrated into the tree traversal (checked before yielding) rather than applied as a post-processing pass over $O(N^2)$ term pairs. For CMPRs with many blocks of distinct sizes and high chaining degrees, this can be orders of magnitude faster.

## Connections

- [Root Expressions and LC Estimation](../theory/Root%20Expressions%20and%20LC%20Estimation.md) — the five propositions that the mesh relies on (especially Props 2 and 3)
- [Coset Regions](../conventions/Coset%20Regions.md) — the geometric picture of weight distributions as grid positions
- [Root Multiplicities and Jordan Decomposition](../theory/Root%20Multiplicities%20and%20Jordan%20Decomposition.md) — why same-field products (repeated sizes) require Jordan partitions that the mesh cannot provide
