# Coset Regions

Cross-field products of coset classes have a natural geometric interpretation: each component field defines an axis, coset weight is the coordinate, and a product of coset classes sweeps out a rectangular **region** in this multi-dimensional space. This picture makes the five root-expression propositions visually intuitive and provides a shared vocabulary for discussing root expression structure.

For the formal definitions of cyclotomic cosets, coset weight, and coset classes $\langle e \cdot w \rangle$, see [Root Expressions — Coset Classes](../theory/Root%20Expressions%20and%20LC%20Estimation.md#coset-classes). For the five propositions that govern root expression arithmetic, see [Root Expressions — Propositions](../theory/Root%20Expressions%20and%20LC%20Estimation.md#the-five-root-expression-propositions). For the exponential notation $\alpha^{p/q}$ used in the examples, see [Roots and the Algebraic Closure](../theory/Roots%20and%20the%20Algebraic%20Closure.md#the-notation-alphapq).

## The Geometric Picture

Consider a CMPR with two MPR blocks of sizes 3 and 5. Each block's roots live in $\mathbb{F}_{2^3}$ and $\mathbb{F}_{2^5}$ respectively. Since $\gcd(3,5) = 1$, these fields are independent — a root from one determines nothing about roots from the other. We draw them as two axes of a grid, with one position per coset weight class.

### The Trivial Coset and Single-Field Roots

Each axis carries $e$ positions: the identity element $\{1\}$, followed by weight classes $1, 2, \ldots, e-1$. The identity gets its own position because the weight-$e$ exponent — the all-ones bit pattern $\underbrace{11\ldots 1}_e$ — wraps to $0$ modulo $2^e - 1$, giving $\alpha^0 = 1$. In the grid we pull this element to the front and separate it from the nontrivial weights with a dashed line:

<div align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 345 175" width="580" style="max-width:100%">
  <rect width="345" height="175" fill="white" rx="5"/>
  <text x="205" y="16" text-anchor="middle" font-family="Georgia,serif" font-size="11.5" fill="#64748b" font-style="italic">F₂⁵ coset weight →</text>
  <text x="14" y="112" text-anchor="middle" font-family="Georgia,serif" font-size="11.5" fill="#64748b" font-style="italic" transform="rotate(-90,14,112)">F₂³ weight</text>
  <text x="83" y="40" text-anchor="middle" font-family="Georgia,serif" font-size="11" fill="#78350f">{1}</text>
  <text x="147" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">1</text>
  <text x="199" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">2</text>
  <text x="251" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">3</text>
  <text x="303" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">4</text>
  <text x="50" y="69" text-anchor="end" font-family="Georgia,serif" font-size="11" fill="#78350f">{1}</text>
  <text x="50" y="117" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">1</text>
  <text x="50" y="153" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">2</text>
  <line x1="115" y1="42" x2="115" y2="168" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <line x1="52" y1="89" x2="330" y2="89" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <rect x="58" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="174" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="226" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="278" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="58" y="96" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="174" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="226" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="278" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="58" y="132" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="174" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="226" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="278" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <text x="83" y="69" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">1</text>
  <text x="147" y="69" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">5</text>
  <text x="199" y="69" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">10</text>
  <text x="251" y="69" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">10</text>
  <text x="303" y="69" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">5</text>
  <text x="83" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">3</text>
  <text x="147" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">15</text>
  <text x="199" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">30</text>
  <text x="251" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">30</text>
  <text x="303" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">15</text>
  <text x="83" y="153" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#78350f">3</text>
  <text x="147" y="153" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">15</text>
  <text x="199" y="153" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">30</text>
  <text x="251" y="153" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">30</text>
  <text x="303" y="153" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">15</text>
</svg>
</div>

Each cell at position $(p_1, p_2)$ contains the product of the per-axis element counts. Position $\{1\}$ contributes $\binom{e}{e} = 1$ element (the identity), while position $w$ contributes $\binom{e}{w}$ elements.

The warm-shaded cells along the $\{1\}$ row and column are **single-field roots**: a cell at $(\{1\}, w)$ has a trivial $\mathbb{F}_{2^3}$ component — the product is just a weight-$w$ root from $\mathbb{F}_{2^5}$ alone. Symmetrically, $\{1\}$-column cells are pure $\mathbb{F}_{2^3}$ roots. The corner $(\{1\}, \{1\})$ is the single element $1 \cdot 1 = 1$. Everything past the dashed line in both directions is a genuine cross-field product.

### Coset Classes as Intervals

A single coset class $\langle e \cdot w \rangle$ is an interval on its axis — it selects all weights from $1$ up to $w$, leaving the $\{1\}$ position empty. The exception is the **full coset** $\langle e \cdot e \rangle$: since the weight-$e$ exponent wraps to the identity, the full coset crosses the gap to include $\{1\}$, covering every element of $\mathbb{F}_{2^e}^*$. On the $\mathbb{F}_{2^3}$ axis:

<div align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 310 142" width="520" style="max-width:100%">
  <rect width="310" height="142" fill="white" rx="5"/>
  <text x="79" y="18" text-anchor="middle" font-family="Georgia,serif" font-size="10" fill="#78350f">{1}</text>
  <text x="151" y="18" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">w = 1</text>
  <text x="211" y="18" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">w = 2</text>
  <line x1="115" y1="22" x2="115" y2="136" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <text x="44" y="47" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·1⟩</text>
  <rect x="50" y="28" width="58" height="30" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="28" width="58" height="30" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="28" width="58" height="30" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <text x="151" y="47" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="258" y="47" font-family="sans-serif" font-size="10" fill="#64748b">3 roots</text>
  <text x="44" y="85" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·2⟩</text>
  <rect x="50" y="66" width="58" height="30" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="66" width="58" height="30" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="66" width="58" height="30" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <text x="151" y="85" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="211" y="85" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="258" y="85" font-family="sans-serif" font-size="10" fill="#64748b">6 roots</text>
  <text x="44" y="123" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·3⟩</text>
  <rect x="50" y="104" width="58" height="30" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="122" y="104" width="58" height="30" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="104" width="58" height="30" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <text x="79" y="123" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">1</text>
  <text x="151" y="123" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="211" y="123" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="258" y="123" font-family="sans-serif" font-size="10" fill="#64748b">7 roots</text>
</svg>
</div>

These counts match Proposition 1: $|\langle 3 \cdot 1 \rangle| = \binom{3}{1} = 3$, $|\langle 3 \cdot 2 \rangle| = \binom{3}{1} + \binom{3}{2} = 6$, $|\langle 3 \cdot 3 \rangle| = 2^3 - 1 = 7$. Notice that $\langle 3 \cdot 3 \rangle$ is the only interval that fills the $\{1\}$ cell — it spans the gap to cover the full multiplicative group $\mathbb{F}_{2^3}^*$.

### Cross-Field Products as Rectangles

A product of coset classes from different fields sweeps out a rectangle. The product $\langle 3 \cdot 1 \rangle \times \langle 5 \cdot 2 \rangle$ covers all cells with $\mathbb{F}_{2^3}$ weight $\leq 1$ and $\mathbb{F}_{2^5}$ weight $\leq 2$:

<div align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 345 175" width="580" style="max-width:100%">
  <rect width="345" height="175" fill="white" rx="5"/>
  <text x="205" y="16" text-anchor="middle" font-family="Georgia,serif" font-size="11.5" fill="#64748b" font-style="italic">F₂⁵ coset weight →</text>
  <text x="14" y="112" text-anchor="middle" font-family="Georgia,serif" font-size="11.5" fill="#64748b" font-style="italic" transform="rotate(-90,14,112)">F₂³ weight</text>
  <text x="83" y="40" text-anchor="middle" font-family="Georgia,serif" font-size="11" fill="#78350f">{1}</text>
  <text x="147" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">1</text>
  <text x="199" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">2</text>
  <text x="251" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">3</text>
  <text x="303" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">4</text>
  <text x="50" y="69" text-anchor="end" font-family="Georgia,serif" font-size="11" fill="#78350f">{1}</text>
  <text x="50" y="117" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">1</text>
  <text x="50" y="153" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">2</text>
  <line x1="115" y1="42" x2="115" y2="168" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <line x1="52" y1="89" x2="330" y2="89" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <rect x="58" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="174" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="226" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="278" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="58" y="96" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="96" width="50" height="34" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1.2"/>
  <rect x="174" y="96" width="50" height="34" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1.2"/>
  <rect x="226" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="278" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="58" y="132" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="174" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="226" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="278" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <text x="147" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">15</text>
  <text x="199" y="117" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">30</text>
</svg>
</div>

$$\text{Count} = |\langle 3 \cdot 1 \rangle| \times |\langle 5 \cdot 2 \rangle| = 3 \times 15 = 45 \text{ roots}$$

This is Proposition 3: the count of a cross-field product is the product of the individual counts. The rectangle lies entirely past both dashed lines — neither coset class includes the identity, so no single-field roots appear.

### Same-Field Products Extend Intervals

The same-field product $\langle 3 \cdot 1 \rangle \cdot \langle 3 \cdot 1 \rangle = \langle 3 \cdot 2 \rangle$ extends the interval along the $\mathbb{F}_{2^3}$ axis (Proposition 2). Geometrically, multiplying two roots from the same field adds their exponents, which can increase the coset weight:

<div align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 300 100" width="500" style="max-width:100%">
  <rect width="300" height="100" fill="white" rx="5"/>
  <text x="79" y="16" text-anchor="middle" font-family="Georgia,serif" font-size="10" fill="#78350f">{1}</text>
  <text x="151" y="16" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">w = 1</text>
  <text x="211" y="16" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">w = 2</text>
  <line x1="115" y1="20" x2="115" y2="96" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <text x="44" y="39" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·1⟩</text>
  <rect x="50" y="22" width="58" height="28" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="22" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="22" width="58" height="28" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <text x="151" y="40" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="150" y="62" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#64748b">↓  · ⟨3·1⟩</text>
  <text x="44" y="83" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·2⟩</text>
  <rect x="50" y="68" width="58" height="28" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="68" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="68" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <text x="151" y="86" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="211" y="86" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
</svg>
</div>

The weight saturates at $e$: $\langle 3 \cdot 2 \rangle \cdot \langle 3 \cdot 2 \rangle = \langle 3 \cdot 3 \rangle$, not $\langle 3 \cdot 4 \rangle$. Visually, saturation means the interval has reached the full coset — it crosses the gap to include $\{1\}$:

<div align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 300 100" width="500" style="max-width:100%">
  <rect width="300" height="100" fill="white" rx="5"/>
  <text x="79" y="16" text-anchor="middle" font-family="Georgia,serif" font-size="10" fill="#78350f">{1}</text>
  <text x="151" y="16" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">w = 1</text>
  <text x="211" y="16" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#334155">w = 2</text>
  <line x1="115" y1="20" x2="115" y2="96" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <text x="44" y="39" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·2⟩</text>
  <rect x="50" y="22" width="58" height="28" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="22" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="22" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <text x="151" y="40" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="211" y="40" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="150" y="62" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#64748b">↓  · ⟨3·2⟩  (saturates)</text>
  <text x="44" y="83" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">⟨3·3⟩</text>
  <rect x="50" y="68" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="122" y="68" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <rect x="182" y="68" width="58" height="28" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1"/>
  <text x="79" y="86" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">1</text>
  <text x="151" y="86" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
  <text x="211" y="86" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#1e40af">3</text>
</svg>
</div>

Once a coset class fills $\{1\}$, that field can contribute the identity to cross-field products — meaning a root from the other field can appear on its own. This is why the full coset is the saturation point: $\langle e \cdot e \rangle = \mathbb{F}_{2^e}^*$.

## Regions in Root Expressions

A root expression is a sum of product terms: $\sum_j \prod_i \langle e_{ji} \cdot w_{ji} \rangle$. Each product term is a rectangle in the grid; the full expression is their union. The propositions give the rules for computing with these unions.

### Example: Two Overlapping Rectangles

Suppose a chaining function produces the root expression $\langle 3 \cdot 1 \rangle \langle 5 \cdot 2 \rangle + \langle 3 \cdot 2 \rangle \langle 5 \cdot 1 \rangle$. Term A ($\langle 3 \cdot 1 \rangle \langle 5 \cdot 2 \rangle$) covers weight $\leq 1$ on the $\mathbb{F}_{2^3}$ axis and weight $\leq 2$ on $\mathbb{F}_{2^5}$; Term B ($\langle 3 \cdot 2 \rangle \langle 5 \cdot 1 \rangle$) covers weight $\leq 2$ and weight $\leq 1$ respectively. Their union:

<div align="center">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 345 210" width="580" style="max-width:100%">
  <rect width="345" height="210" fill="white" rx="5"/>
  <text x="205" y="16" text-anchor="middle" font-family="Georgia,serif" font-size="11.5" fill="#64748b" font-style="italic">F₂⁵ coset weight →</text>
  <text x="14" y="112" text-anchor="middle" font-family="Georgia,serif" font-size="11.5" fill="#64748b" font-style="italic" transform="rotate(-90,14,112)">F₂³ weight</text>
  <text x="83" y="40" text-anchor="middle" font-family="Georgia,serif" font-size="11" fill="#78350f">{1}</text>
  <text x="147" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">1</text>
  <text x="199" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">2</text>
  <text x="251" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">3</text>
  <text x="303" y="40" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#334155">4</text>
  <text x="50" y="69" text-anchor="end" font-family="Georgia,serif" font-size="11" fill="#78350f">{1}</text>
  <text x="50" y="117" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">1</text>
  <text x="50" y="153" text-anchor="end" font-family="sans-serif" font-size="11" fill="#334155">2</text>
  <line x1="115" y1="42" x2="115" y2="168" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <line x1="52" y1="89" x2="330" y2="89" stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>
  <rect x="58" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="174" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="226" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="278" y="48" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="58" y="96" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="96" width="50" height="34" rx="3" fill="#c4b5fd" stroke="#7c3aed" stroke-width="1.2"/>
  <rect x="174" y="96" width="50" height="34" rx="3" fill="#93c5fd" stroke="#3b82f6" stroke-width="1.2"/>
  <rect x="226" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="278" y="96" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="58" y="132" width="50" height="34" rx="3" fill="#fef9ef" stroke="#d4a574" stroke-width="0.8"/>
  <rect x="122" y="132" width="50" height="34" rx="3" fill="#fdba74" stroke="#f97316" stroke-width="1.2"/>
  <rect x="174" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="226" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <rect x="278" y="132" width="50" height="34" rx="3" fill="#eef2f7" stroke="#94a3b8" stroke-width="0.8"/>
  <text x="147" y="115" text-anchor="middle" font-family="sans-serif" font-size="9" fill="#5b21b6">A ∩ B</text>
  <text x="199" y="115" text-anchor="middle" font-family="sans-serif" font-size="9" fill="#1e40af">A</text>
  <text x="147" y="151" text-anchor="middle" font-family="sans-serif" font-size="9" fill="#9a3412">B</text>
  <rect x="55" y="180" width="12" height="12" rx="2" fill="#93c5fd" stroke="#3b82f6" stroke-width="0.8"/>
  <text x="72" y="190" font-family="sans-serif" font-size="10" fill="#334155">A: ⟨3·1⟩⟨5·2⟩</text>
  <rect x="155" y="180" width="12" height="12" rx="2" fill="#fdba74" stroke="#f97316" stroke-width="0.8"/>
  <text x="172" y="190" font-family="sans-serif" font-size="10" fill="#334155">B: ⟨3·2⟩⟨5·1⟩</text>
  <rect x="260" y="180" width="12" height="12" rx="2" fill="#c4b5fd" stroke="#7c3aed" stroke-width="0.8"/>
  <text x="277" y="190" font-family="sans-serif" font-size="10" fill="#334155">A ∩ B</text>
</svg>
</div>

The overlap (Proposition 5) is $\langle 3 \cdot \min(1,2) \rangle \langle 5 \cdot \min(2,1) \rangle = \langle 3 \cdot 1 \rangle \langle 5 \cdot 1 \rangle$, the single cell at $(1, 1)$ — exactly the component-wise minimum of the weight bounds:

$$|A \cup B| = |A| + |B| - |A \cap B| = 45 + 30 - 15 = 60$$

## Concrete Root Example

To connect the geometry back to actual roots: consider the roots of two MPR blocks of sizes 3 and 7 (fields $\mathbb{F}_{2^3}$ and $\mathbb{F}_{2^7}$, orders 7 and 127).

In the exponential notation, $\alpha^{1/7}$ is a weight-1 root from $\mathbb{F}_{2^3}$ (exponent 1 has binary representation $001$). The root $\alpha^{5/7}$ is weight-2 (exponent 5 = $101$ in binary). Similarly, $\alpha^{1/127}$ is weight-1 from $\mathbb{F}_{2^7}$.

Their cross-field products:

| $\mathbb{F}_{2^3}$ root | position | $\mathbb{F}_{2^7}$ root | position | Product | Grid position |
|---|---|---|---|---|---|
| $1$ | $\{1\}$ | $\alpha^{1/127}$ | 1 | $\alpha^{1/127}$ | $(\{1\}, 1)$ |
| $\alpha^{1/7}$ | 1 | $\alpha^{1/127}$ | 1 | period $\text{lcm}(7, 127)$ | $(1, 1)$ |
| $\alpha^{5/7}$ | 2 | $\alpha^{1/127}$ | 1 | period $\text{lcm}(7, 127)$ | $(2, 1)$ |
| $\alpha^{1/7}$ | 1 | $\alpha^{3/127}$ | 2 | period $\text{lcm}(7, 127)$ | $(1, 2)$ |

The first row sits on the $\{1\}$ row of the grid — the $\mathbb{F}_{2^3}$ component is the identity, so the product is just a root from $\mathbb{F}_{2^7}$ alone. In a root expression, this root would only be covered if the $\mathbb{F}_{2^3}$ coset class is the full coset $\langle 3 \cdot 3 \rangle$ (the only interval that crosses the gap to include $\{1\}$).

The remaining products land in the interior of the grid — genuine cross-field roots. The coset class product $\langle 3 \cdot 1 \rangle \times \langle 7 \cdot 1 \rangle$ covers the $(1,1)$ cell — all $3 \times 7 = 21$ cross-products of weight-1 roots from each field. The product $\langle 3 \cdot 2 \rangle \times \langle 7 \cdot 1 \rangle$ extends the rectangle down to include the $(2,1)$ cell, adding $3 \times 7 = 21$ more roots.

The "region" is the rectangle that a coset class product sweeps out in this grid. It is a compact way to describe which cross-field root combinations are possible in a signal, without enumerating individual roots.

## Higher Dimensions

For a CMPR with $k$ blocks of distinct sizes, the grid has $k$ axes — each with its own $\{1\}$ position and dashed gap. A product term $\prod_{i=1}^{k} \langle s_i \cdot w_i \rangle$ is a $k$-dimensional box (hyperrectangle). The propositions generalize naturally:

- **Proposition 3** says the volume (root count) of a box is the product of its side lengths.
- **Proposition 5** says the intersection of two boxes is the box with component-wise minimum side lengths.
- **Proposition 4** says the union count uses inclusion-exclusion over the intersections.

The CMPR root expression algorithm builds up the root expression table block by block. Each block adds one axis to the grid — the new block's $\langle s_i \cdot 1 \rangle$ — and the chaining function determines how the new axis combines with the existing ones.

A product term with $\langle s_i \cdot s_i \rangle$ on one axis extends across the gap on that axis, meaning the corresponding field can contribute $1$ to the product. In a $k$-dimensional grid, faces of the hyperrectangle that touch a $\{1\}$ boundary represent roots from fewer than $k$ subfields.

## Terminology Summary

| Term | Meaning |
|---|---|
| **Axis** | One component field $\mathbb{F}_{2^{s_i}}$ of the CMPR |
| **$\{1\}$ position** | The identity element on each axis — the weight-$e$ exponent wrapping to $\alpha^0 = 1$ |
| **Gap** | The dashed separator between $\{1\}$ and the nontrivial weight positions |
| **Coordinate** | Coset weight along that axis ($1$ to $e-1$, plus the $\{1\}$ position) |
| **Cell** | A specific weight pair/tuple — the set of cross-product roots at those weights |
| **Single-field cell** | A cell on a $\{1\}$ row/column — its root lives in one subfield only |
| **Region** | A rectangle (or hyperrectangle) in the grid — the set of cells covered by a coset class product |
| **Root expression** | A union of regions — the full set of possible roots in a signal |
| **Full coset** | $\langle e \cdot e \rangle = \mathbb{F}_{2^e}^*$ — the only coset class that crosses the gap |
