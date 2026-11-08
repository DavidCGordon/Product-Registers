# Division Property

The division property is a way to describe what algebraic information remains in a multiset of intermediate values without enumerating every value or every monomial. It grew out of integral cryptanalysis: instead of recording only that an XOR sum is balanced, it records a family of monomial moments that are guaranteed to vanish.

This document uses the bit-subset order and matrix orientation from [Matrix Indexing](../conventions/Matrix%20Indexing.md), [Algebraic Normal Form](Algebraic%20Normal%20Form.md), and [Möbius Inversion on the Bit-Subset Lattice](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md). The original division-property reference is Todo's EUROCRYPT 2015 paper, [Structural Evaluation by Generalized Integral Property](../conventions/bibliography.bib).

## History and Motivation

Integral attacks study a collection of inputs, often called a structure, and test XOR sums of selected output coordinates. The classical integral statement is therefore a statement about moments of a multiset. Todo's 2015 division property generalized this viewpoint from a single integral distinction to a partially ordered family of algebraic degrees. A division-property label records which monomial moments are forced to vanish, allowing that information to be propagated through several rounds of a cipher.

The later literature developed several directions from this starting point:

- **Bit-based division property** labels each bit separately. A mask in $\{0,1\}^n$ says which individual coordinates participate in the tracked monomial.
- **Word-based division property** groups bits into words and tracks a weight or degree per word. This reduces the state space when a primitive is naturally word-oriented, at the cost of forgetting which particular bits inside a word occur.
- **Generalized division property** allows non-binary, usually non-negative integer labels. It can express degree or multiplicity information that the original bit mask cannot distinguish, and it contains the bit-based version as the special case of labels in $\{0,1\}$.
- **Division trails** turn propagation into a relation between input and output labels. For a linear or nonlinear component, a trail records one soundly allowed transition; a collection of trails gives an upper bound on the labels that can survive.

These names are related but are not interchangeable. Bit-based and word-based descriptions use different state spaces, while generalized descriptions change the label set and therefore the partial order. A propagation rule is meaningful only after the label type and the direction of the transition have been specified.

## The Poset Definition

Let $X$ be a multiset of vectors in $\mathbb F_2^n$. For a mask $u\in\{0,1\}^n$, define its monomial moment by

$$
M_u(X)=\bigoplus_{x\in X}\prod_{i=0}^{n-1}x_i^{u_i}.
$$

The product is 1 when every coordinate selected by $u$ is 1, and 0 otherwise. Thus $M_u(X)$ is the XOR sum of the monomial indexed by $u$ over the multiset.

Write $u\preccurlyeq v$ when $u_i\leq v_i$ for every coordinate. For bit masks this is exactly the [bit-subset relation](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md#the-bit-subset-partial-order): $u\preccurlyeq v$ means that every 1-bit of $u$ is also a 1-bit of $v$.

A multiset has the **bit-based division property** $D_\alpha$ when

$$
M_u(X)=0\qquad\text{for every }u\text{ with }\alpha\not\preccurlyeq u.
$$

Equivalently, only moments whose masks contain $\alpha$ are allowed to remain nonzero. Larger labels are stronger statements: if $\beta\preccurlyeq\alpha$, then $D_\alpha$ implies $D_\beta$. Indeed, a mask $u$ with $\beta\not\preccurlyeq u$ also has $\alpha\not\preccurlyeq u$, since $\alpha\preccurlyeq u$ would give $\beta\preccurlyeq\alpha\preccurlyeq u$; so $D_\alpha$ forces $M_u(X)=0$. The label $0^n$ imposes no vanishing condition; the label $1^n$ says that every moment except the full-product moment vanishes.

The labels are therefore ordered by $\preccurlyeq$, the coordinatewise containment of masks, and not by ordinary integer size.

## Three Vectors and Two Conversions

In this section a mask is written as the set of its 1-coordinates: $S$ stands for the mask $u$ whose 1-bits are the elements of $S$, $x^S=\prod_{i\in S}x_i$ is its monomial, and $M_S$ is the moment $M_u(X)$. The structure is the multiset $X$.

There are three vectors that are easy to confuse:

1. **Input-set vector.** Let $q_T$ indicate whether the structure contains the input mask $T$ with odd multiplicity.
2. **Moment vector.** Let $M_S$ be the XOR sum of the monomial $x^S$ over the structure.
3. **ANF coefficient vector.** Let $a_S$ be the coefficient of the monomial $x^S$ in a Boolean function.

For an input mask $T$, the monomial $x^S$ evaluates to 1 exactly when $S\subseteq T$. Therefore input multiplicities and moments are related by the upward zeta transform:

$$
M_S=\bigoplus_{T\supseteq S}q_T.
$$

In matrix form this is $M=\mathsf C_n q$, where $\mathsf C_n$ is the upward zeta transform of [Dihedral Symmetries of the Boolean Transform](Dihedral%20Symmetries%20of%20the%20Boolean%20Transform.md#the-four-basic-operations): the $2^n\times2^n$ matrix over $\mathbb F_2$, acting on column vectors indexed by masks, with $(\mathsf C_n v)_S=\bigoplus_{T\supseteq S}v_T$. The subscript is the number of variables (that document writes it $m$). Since $\mathsf C_n$ is self-inverse over $\mathbb F_2$, the same transform recovers the parity of the input-set vector from the exact moment vector.

For a Boolean function, ANF coefficients and truth-table values instead satisfy

$$
 f(T)=\bigoplus_{S\subseteq T}a_S,
$$

which is $f=\mathsf B_n a$ for the downward zeta transform $\mathsf B_n$, the matrix with $(\mathsf B_n v)_T=\bigoplus_{S\subseteq T}v_S$. This is the Möbius transform of [Algebraic Normal Form](Algebraic%20Normal%20Form.md#the-transform-matrix). Like $\mathsf C_n$, the matrix $\mathsf B_n$ is its own inverse over $\mathbb F_2$.

So there are two zeta conversions:

$$
\boxed{\text{ANF coefficients }a\xrightarrow{\mathsf B_n}\text{ truth values }f}
$$

and

$$
\boxed{\text{input-set parities }q\xrightarrow{\mathsf C_n}\text{ monomial moments }M}.
$$

The two matrices are transposes of each other. The entry of $\mathsf B_n$ in row $T$ and column $S$ is 1 exactly when $S\subseteq T$, and the entry of $\mathsf C_n$ in row $S$ and column $T$ is 1 under the same condition, so $\mathsf C_n=\mathsf B_n^{\top}$. The next section shows that $\mathsf C_n$ is also the conjugate of $\mathsf B_n$ by mask complementation.

The vectors $q$ and $M$ are equivalent descriptions of the structure $X$, since $M=\mathsf C_n q$ and $\mathsf C_n$ is invertible. The vector $a$ belongs to a Boolean function, and this section does not relate that function to $X$.

> **OPEN:** Is the Boolean function of item 3 meant to be the indicator of the structure, $f=q$, so that $a=\mathsf B_n q$? — needs author input.
> Established so far: $M=\mathsf C_n q$ and $f=\mathsf B_n a$, each derived above; nothing in this document identifies $f$ with $q$.

## Complement and the Dual Description

The antidiagonal transform $J_n$ is the permutation matrix that sends a mask $S$ to its complement $S^{\mathsf c}$, so that $(J_n v)_S=v_{S^{\mathsf c}}$. Conjugating a zeta transform by $J_n$ exchanges subsets and supersets. The entry of $J_n\mathsf B_nJ_n$ in row $S$ and column $T$ is the entry of $\mathsf B_n$ in row $S^{\mathsf c}$ and column $T^{\mathsf c}$, which is 1 exactly when $T^{\mathsf c}\subseteq S^{\mathsf c}$, that is, when $S\subseteq T$. This is the entry of $\mathsf C_n$, so

$$
\mathsf C_n=\mathsf B_n^{\top}=J_n\mathsf B_nJ_n.
$$

Consequently, a division-property calculation written in terms of present variables and one written in terms of absent variables are related by the same complement operation that appears in the [Cayley diagram](Dihedral%20Symmetries%20of%20the%20Boolean%20Transform.md#cayley-diagram) of the group generated by $\mathsf B_n$ and $J_n$. The two order-3 compositions are

$$
\mathsf A_n=\mathsf B_nJ_n,
\qquad
\mathsf D_n=J_n\mathsf B_n.
$$

The factor on the right acts first on a column vector. Thus $\mathsf A_n$ means “complement, then downward zeta,” while $\mathsf D_n$ means “downward zeta, then complement.” These descriptions are exact for coefficient, value, input-parity, or moment vectors.

A division-property calculation tracks only which coordinates of such a vector can be nonzero, and for supports the matrices give an inclusion rather than an equality. Write $\operatorname{supp}v$ for the set of coordinates where $v$ is 1 and $e_i$ for the $i$-th standard basis vector. For a 0/1 matrix $\mathsf X$ we have $\mathsf Xv=\bigoplus_{i\in\operatorname{supp}v}\mathsf Xe_i$, so

$$
\operatorname{supp}(\mathsf Xv)\subseteq\bigcup_{i\in\operatorname{supp}v}\operatorname{supp}(\mathsf Xe_i).
$$

The inclusion can be strict, because two columns $\mathsf Xe_i$ that are both 1 in the same coordinate cancel there.

> **OPEN:** What map do $\mathsf A_n$ and $\mathsf D_n$ induce on division-property labels $\alpha$? — needs author input.
> Established so far: their exact action on vectors, and the support inclusion above. No rule sending a label $D_\alpha$ to a label is stated in this document.

## Types of Division Property

### Bit-based

A bit-based label is $\alpha\in\{0,1\}^n$. It distinguishes individual coordinates and uses the Boolean-lattice order $\preccurlyeq$. It is the most direct form for bit-oriented Boolean circuits and is the form defined above.

Its main strength is precision about individual variables. Its main cost is a label space of size $2^n$, which becomes large even when the cipher operates on a small number of words.

### Word-based

Partition the state into words $W_1,\ldots,W_r$. A word-based label replaces the detailed bit mask inside each word by a quantity such as a Hamming-weight or degree bound. A label might therefore look like $(d_1,\ldots,d_r)$, where $d_j$ records how much algebraic degree may be drawn from word $W_j$.

The word-based representation is a quotient of the bit-based one: many bit masks map to the same word label. This makes propagation cheaper and often matches word-level components, but it cannot distinguish two masks with the same per-word weight.

### Generalized

A generalized label uses a larger ordered set, commonly $\mathbb Z_{\geq0}^n$ with componentwise order. The extra integer values can record repeated degree demand, multiplicity, or other information that a binary presence/absence mask cannot express. The exact meaning of a generalized label depends on the definition adopted by the particular paper or implementation; the order and vanishing predicate must be stated with it.

The bit-based property embeds into this setting by restricting every coordinate to 0 or 1. Generalization is therefore not merely a larger bit mask: it changes which labels are comparable and which moments are required to vanish.

## Propagation Rules

A propagation rule is a sound implication from an input property to an output property. At the exact polynomial level, the rules come from expanding the component's ANF. At the division-property level, one keeps every output label that could arise and usually discards coefficient cancellations.

### Constants

The constant 0 contributes no monomial. The constant 1 contributes only the empty monomial. XOR with a constant can therefore change the empty-mask coefficient, but it does not introduce a nonempty mask.

### Copying and permutations

A wire copy transports the corresponding coordinate. A permutation $\pi$ sends a mask $S$ to $\pi(S)$. In matrix language, a permutation matrix simply permutes the mask coordinates; in a division trail, this gives a deterministic relabelling of the mask.

### XOR and linear maps

For an XOR gate, expand the selected output monomial into products of input terms. A linear map is handled in the same way: for an output mask $u$, compute the ANF of the product of the selected output coordinates, then record every input mask that appears. The resulting input-to-output relation is the exact division-trail relation for that linear component.

This direction matters. When analysing a forward computation $y=L(x)$, one often asks which input masks can contribute to an output mask, so the relevant expansion is the pullback of output monomials through $L$. When analysing a reverse or dual description, the transpose/zeta orientation appears instead. The matrix convention prevents the two directions from being silently identified.

### AND and multiplication

At the monomial level, multiplying $x^S$ and $x^T$ gives $x^{S\cup T}$ because Boolean variables are idempotent. Thus an AND combines masks by set union:

$$
 x^Sx^T=x^{S\cup T}.
$$

On the bit masks $s$ and $t$ of $S$ and $T$, the union is the bitwise OR $s\mathbin{\mathrm{OR}}t$.

For generalized degree labels, the corresponding operation is usually componentwise addition followed by whatever saturation or truncation the definition prescribes. The bit-based OR rule and the generalized addition rule should not be conflated.

### Nonlinear components and S-boxes

For a nonlinear component, derive the transition relation from its coordinate ANFs or from a precomputed local table. For every output mask, expand the selected output monomial in the input variables and retain the input masks that occur. A division-property implementation may then prune transitions using symmetry, word bounds, or known impossible cases, provided the pruning is sound.

This local rule is the part that depends most strongly on the chosen division-property type. A bit-based rule may distinguish every input mask; a word-based rule merges many such transitions; a generalized rule may retain degree values above one.

## What the Transform Does Not Replace

The Boolean-lattice transforms explain how to change coordinates between input-set parities, monomial moments, ANF coefficients, and truth-table values. They do not replace component propagation rules. In particular:

- A zeta transform does not determine the nonlinear behaviour of an S-box.
- A complemented mask is not automatically the same cryptanalytic statement as the original mask; it is a statement in a different presence/absence convention.
- Support propagation is an upper bound whenever different exact coefficients can cancel.
- Word-based and generalized properties require their own label semantics: each has its own label set and its own order, as described in [Types of Division Property](#types-of-division-property).

The practical rule is to state three things whenever a division property is used: the label space, the partial order, and the direction in which a transition is being propagated. Once those are fixed, the upward/downward zeta transforms and the $J_n$ complement give a consistent way to move between the corresponding descriptions.

## References

- Yosuke Todo, “Structural Evaluation by Generalized Integral Property,” EUROCRYPT 2015. See the [bibliography entry](../conventions/bibliography.bib) for publication details.
- [Algebraic Normal Form](Algebraic%20Normal%20Form.md) for the exact coefficient-to-value transform.
- [Möbius Inversion on the Bit-Subset Lattice](Moebius%20Inversion%20on%20the%20Bit-Subset%20Lattice.md) for the poset proof of the zeta/Möbius identities.
- [Dihedral Symmetries of the Boolean Transform](Dihedral%20Symmetries%20of%20the%20Boolean%20Transform.md) for the six compositions generated by zeta transformation and complementation.
