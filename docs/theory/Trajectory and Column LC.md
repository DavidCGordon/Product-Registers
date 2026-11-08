# Trajectory LC vs Column LC

Two different linear complexities arise in the analysis of [algebraic attacks](Algebraic%20Attacks.md) on CMPRs (registers composed of several MPR blocks, with chaining from upstream blocks into downstream ones): the **trajectory LC** of a signal and the **column LC** of a monomial in the attack's equation matrix. Across the block structure of a CMPR the two have been observed to run in opposite directions: where one is small the other is large. The relationship matters for predicting the rank of the equation matrix in the Fast Algebraic Attack (FAA). This document defines the two quantities, separates what is derived about their relationship from what is observed, and records what is still open.

Throughout, $s_t$ is the register state at clock cycle $t$, a polynomial function of the initial state $s_0$, and $g$ is a Boolean function of the state. For a set $I$ of initial-state bit positions, $T_I = \prod_{i \in I} x_i$ is the corresponding monomial in the initial-state variables. The columns of an attack's equation matrix are indexed by the initial-state monomials that the attack keeps as unknowns, and reading down a column gives a sequence indexed by $t$.

## Definitions

The **trajectory LC** of a state bit or monomial is the linear complexity of the sequence $g(s_t)$ when $g$ is evaluated along a trajectory from some specific initial state. It can be estimated from the bit's [RootExpression](Root%20Expressions%20and%20LC%20Estimation.md) via `re.upper()`.

The **column LC** of an initial-state monomial $T_I$ is the linear complexity of the column sequence $c_I(t)$ in the equation matrix, where $c_I(t)$ is the coefficient of $T_I$ in $g(s_t)$. It measures how fast that coefficient oscillates as $t$ increases. Column LC is not the same quantity as trajectory LC.

The [cube-sum identity](Cube%20Equation%20Generation.md) expresses the column sequence through trajectories:

$$c_I(t) = \bigoplus_{v \in \{0,1\}^{|I|}} g(s_t^v)$$

where $s_t^v$ is the trajectory starting from the state where bits in $I$ are set to $v$, all others zero. For a single bit, $I = \{j\}$, the sum has two terms, $c_j(t) = g(s_t^{e_j}) \oplus g(s_t^{0})$, where $e_j$ is the state with only bit $j$ set.

## The Opposite Trends Across Blocks

The examples below refer to a CMPR with blocks of 7, 5, and 3 bits, written M7→M5→M3, with chaining flowing in the direction of the arrows: M7 is the upstream (source) block and M3 the most downstream.

**The root-expression bound on trajectory LC is small upstream and grows downstream.** Root expressions are upper bounds, so the statement is about the bound, not about the trajectory LC itself. It follows from the block-by-block propagation of root expressions ([Root Expressions — The CMPR Root Expression Algorithm](Root%20Expressions%20and%20LC%20Estimation.md#the-cmpr-root-expression-algorithm)). Write $n_1$ for the size of the source block and $n_i$ for the size of the $i$-th block. A bit of the source block has root expression $\langle n_1 \cdot 1 \rangle$, the block's own roots. A bit of each later block has the root expression of its chaining input, computed from the upstream blocks, plus its own block's $\langle n_i \cdot 1 \rangle$. The bound for a downstream bit therefore accumulates the upstream contributions.

**Column LC is large upstream and small downstream.** This is an observation from the project's experiments, not a derived result: columns of monomials in upstream bits (M7 in the example) were observed to have large LC, and columns of monomials in downstream bits (M3) small LC.

**One case can be derived.** Let $j$ be a bit of the most-downstream block $B$, with update matrix $U$. No other block receives chaining input from $B$, so the trajectories from $e_j$ and from $0$ have identical upstream state at every time, and $B$ receives the same chaining input in both. The difference $\Delta[t]$ of the two $B$-states therefore satisfies $\Delta[t+1] = U\Delta[t]$ with $\Delta[0] = e_j$, so $\Delta[t] = U^t e_j$. If $g$ is an $\mathbb{F}_2$-linear function of the bits of $B$ alone, then

$$c_j(t) = g(s_t^{e_j}) \oplus g(s_t^{0}) = g(U^t e_j),$$

which is an output sequence of the bare MPR with update $U$. Its LC is at most the size of the block.

## The All-Zeros State

The cube sum for a single bit involves the trajectory from the all-zeros state, so it matters whether that state is a fixed point. Block 0, the source, receives no chaining input and stays at zero. The first block to leave zero can do so only if its chaining function is nonzero at the zero input; a later block can also leave zero because a block upstream of it already has. If every chaining function vanishes at zero, the all-zeros state is a fixed point of the CMPR and $g(s_t^0)$ is constant. Otherwise the cube sum XORs $g(s_t^{e_j})$ with a sequence $g(s_t^0)$ that need not be constant, and the resulting column LC can be smaller than the trajectory LC when $g(s_t^{e_j})$ and $g(s_t^0)$ share roots that cancel.

## Practical Implications

`bit_re.upper()` and `product_of_bit_res.upper()` do not predict column LC: in the project's experiments they overestimate downstream column LCs and underestimate upstream column LCs. An analysis of FAA rank loss or of equation matrix structure for a CMPR should use empirical column LC data rather than RootExpression upper bounds. The [monomial profile](Monomial%20Profile%20Theory.md) of a bit tells which monomials *can* appear, but not how their coefficient sequences (columns) behave over time.

## Open Question

This documentation gives no formula for column LC as a function of block sizes, annihilator degree, the `max_and` chaining parameter, and monomial profile. Such a formula is needed to bound the rank of the FAA equation matrix for a CMPR and to predict attack complexity.
