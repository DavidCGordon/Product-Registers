# Writing Standards

Standards for mathematical prose in this repository — theory docs, architecture docs, docstrings, and inline comments that explain non-obvious mathematics.

The target reader is a strong undergraduate in mathematics: solid linear algebra and abstract algebra, no prior exposure to feedback registers or algebraic cryptanalysis. Writing should be rigorous and complete at that level.

## The Core Test

Every sentence that asserts a mathematical relationship **incurs a debt**. The debt is discharged by stating the relationship precisely enough that the reader could verify it, or reconstruct it, without consulting an outside source.

A sentence that asserts a relationship and does not discharge the debt is worse than no sentence at all: it leaves the reader believing something has been established when nothing has. The reader cannot tell whether the gap is theirs or the author's.

This is the test to apply throughout:

> **If I deleted this sentence, would the reader lose information, or only lose the impression that something was explained?**

If the answer is "only the impression," the sentence must either be expanded until it carries content or removed.

### The debt is about content, not vocabulary

None of what follows is a list of forbidden words. Words like *natural*, *canonical*, *dual*, and *analogous* are legitimate technical vocabulary, and in many places they are the correct word. The problem is never the word — it is using the word **in place of** the content it is supposed to summarize.

"The natural map" is fine when the map is written down nearby. It is a debt when the reader is left to guess which map is meant.

## Part 1 — Patterns That Work

These are the constructions to reach for. They are what discharging the debt looks like in practice.

### Definition before use, in place

Give a self-contained definition at the point of first use, even when another document covers it in more depth. A one-sentence definition plus a link is better than a link alone, and much better than neither.

> **Good** — [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md):
>
> "A **Boolean function** on $m$ variables is a function $f: \mathbb{F}_2^m \to \mathbb{F}_2$ — it takes $m$ binary inputs and returns a binary output."

The definition is complete, sits where it is needed, and does not require the reader to leave the page.

### Multiple equivalent definitions, with the equivalence shown

When a concept has several standard characterizations, giving more than one is valuable — different definitions make different properties obvious. This is encouraged, with one requirement: **state that they are equivalent and show why**, or point to where the equivalence is proved.

The ANF doc does this well. It defines the ANF three ways:

1. As a multilinear polynomial: $f = \bigoplus_S a_S \prod_{b \in S} x_b$
2. In the binomial basis: $f(n) = \bigoplus_i a_i \binom{n}{i} \bmod 2$
3. As a matrix identity: $\mathbf{f} = T\mathbf{a}$ with $T_{n,i} = \binom{n}{i} \bmod 2$

The equivalences are not asserted — they are derived. Definition 1 becomes definition 2 through Lucas' theorem, with the bit-subset condition $i \preccurlyeq n$ shown to be the same condition in both. Definition 2 becomes definition 3 by observing that evaluating at every input is a matrix-vector product.

The payoff is concrete: the polynomial form makes the algebraic degree obvious, the binomial form connects to Lucas and decimation identities, and the matrix form makes the self-inverse property a one-line computation. Each definition earns its place by making something easy that the others make hard.

An equivalence stated but not shown is a debt like any other. "These are equivalent" needs either a proof, a computation, or a pointer to the specific section that proves it.

### Analogies that name the correspondence

An analogy is legitimate when it comes with the map. State what corresponds to what, and what property is preserved.

> **Good** — [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md):
>
> "This is the Boolean analog of the FFT butterfly: each Kronecker factor contributes one round of pairwise operations, and the rounds compose to apply the full matrix."

The correspondence is spelled out in the same sentence — Kronecker factor ↔ round of pairwise operations. The reader can check it.

Compare a version that would incur a debt: *"This is the Boolean analog of the FFT."* True, but it teaches nothing. Which structure corresponds to which? Is it the divide-and-conquer, the butterfly, the $O(n \log n)$ bound, or the group-theoretic character sum? The reader has no way to know.

### Constructions written down, not just named

When a construction is claimed to be natural, canonical, obvious, or unique, **exhibit it**. The claim of canonicity is often the interesting part, and it is exactly what a reader cannot reconstruct alone.

> **Good** — [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md) opens by calling the ANF "the canonical way to represent a Boolean function as a polynomial over $\mathbb{F}_2$" — then proves the uniqueness claim: $2^m$ monomials give $2^{2^m}$ candidate polynomials, matching the count of Boolean functions, so uniqueness follows from linear independence, which the triangular evaluation matrix establishes.

The word "canonical" is used correctly here because the document goes on to say in what sense — unique representation — and proves it. If a claim of canonicity cannot be cashed out this way, it should be softened to a description of the specific choice being made and why it was made.

### Mechanisms, not appeals to authority

Favor explicit computation over "by [named theorem]." Naming the theorem is useful as a signpost, but the step must also be shown.

> **Good** — [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md) invokes Lucas' theorem *and* shows the consequence:
>
> "Each factor $\binom{n_b}{i_b}$ is 0 only when $i_b = 1$ and $n_b = 0$. So $\binom{n}{i} \equiv 1 \pmod{2}$ exactly when $i \preccurlyeq n$."

The theorem is named, then the specific instance is worked out. The reader who has never seen Lucas' theorem can still follow the argument.

### Contrasts that locate the boundary

Stating what a result does *not* say is often as valuable as stating what it does. Contrasts prevent overreach.

> **Good** — [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md):
>
> "Over $\mathbb{Z}$, the Möbius inverse of Pascal's triangle has alternating signs; over $\mathbb{F}_2$, those signs collapse and the matrix is its own inverse."

This tells the reader precisely which hypothesis the self-inverse property depends on. A reader who later works over $\mathbb{Z}$ knows immediately that the property fails.

The "What the Transform Does Not Replace" section of [Division Property](../theory/Division%20Property.md) is a good example at section scale — it enumerates what the zeta transforms do not give you, preventing the reader from over-applying them.

### Worked examples with real numbers

A small concrete instance, computed all the way through, catches misunderstandings that prose cannot. The ANF doc's $m = 2$ butterfly example — starting from $\mathbf{f} = (0,1,0,0)$ and running both passes to reach $\mathbf{a} = (0,1,0,1)$ — lets the reader verify their understanding of the algorithm against a known answer.

Prefer examples small enough to check by hand.

## Part 2 — Claims That Incur a Debt

These are heuristics for **locating** passages that need scrutiny. They are neither exhaustive nor decisive. A passage matching one of these patterns is not automatically wrong; it is a passage where the reviewer should ask whether the debt has been discharged. A passage matching none of them can still be vacuous.

The response to a flagged passage is one of three things:

1. **Discharge it** — add the missing correspondence, definition, or computation.
2. **Weaken it** — replace the strong claim with the weaker one that is actually supported.
3. **Delete it** — if the sentence carries no content once the unsupported part is removed.

Rewriting the sentence to sound more hedged is not one of the options. "Arguably the natural analogue" is the same debt with a disclaimer attached.

### Unjustified analogy

**Pattern:** "$X$ is the $Y$ analogue of $Z$," "$X$ is the $Y$-theoretic version of $Z$," "$X$ plays the role of $Z$" — without stating the correspondence.

**Why it misleads:** Analogy in mathematics ranges from "there is a functor" to "these two formulas have a similar shape." Naming an analogy without specifying which end of that range it occupies invites the reader to assume the strong reading.

**To discharge:** Name the map. Say what corresponds to what and what is preserved. If the honest answer is "the formulas look similar," say that instead — a noted resemblance is a legitimate observation when labeled as one.

> **Needs work** — [Division Property](../theory/Division%20Property.md):
>
> "This is the division-property analogue of the binomial/Möbius duality. The upward transform is not a new ad hoc operation: it is the same Boolean-lattice zeta transform viewed from the dual direction."
>
> "Viewed from the dual direction" is carrying the entire claim without defining what the dual direction is. The content exists — the next section gives $\mathsf{C}_m = J_m \mathsf{B}_m J_m$, which is exactly the precise statement — but it is asserted before it is shown, so at the point of reading the sentence is an appeal, not an argument. The fix is to state the conjugation identity here, or to defer the claim to the section that proves it.

### Name-drop as explanation

**Pattern:** "by Möbius inversion," "this is just the ANF," "by duality," "standard," "well-known" — invoking a name in place of the step.

**Why it misleads:** The name is a pointer to a result the reader may not have. Even a reader who knows the result may not see which instance is being applied.

**To discharge:** Name the result *and* show the step it licenses. Two extra sentences usually suffice.

### Unstated morphism

**Pattern:** "$X$ is the same as $Y$," "$X$ *is* $Y$," "these coincide" — without saying under what identification.

**Why it misleads:** Two objects can be equal, isomorphic, equal after a choice of basis, or merely in bijection. These are different claims with different consequences.

**To discharge:** State the identification. "Equal as matrices," "isomorphic as $\mathbb{F}_2$-vector spaces via the map $\ldots$," "in bijection but not canonically."

This matters acutely in this codebase, where primal/dual polynomial conventions and matrix orientation are exactly the places where a sloppy "is the same as" produces a real bug. See [Polynomial Conventions](Polynomial%20Conventions.md) and [Matrix Indexing](Matrix%20Indexing.md).

### Profundity adjectives

**Pattern:** *deep*, *elegant*, *beautiful*, *powerful*, *profound*, *remarkable*, *fundamental*, *rich*, *striking* — used to characterize a result.

**Why it misleads:** These words report the author's feeling about a result rather than its content. They tell the reader that something impressive happened without saying what.

**To discharge:** Replace with the specific property that prompted the adjective. "Deep" usually means "the proof is not obvious" or "it connects two areas that look unrelated" — say which. "Powerful" usually means "it applies in more cases than you would expect" — say which cases.

Note that *fundamental* has legitimate technical uses (fundamental theorem, fundamental domain, fundamental group). The heuristic targets the adjectival, evaluative use.

### Connection without mechanism

**Pattern:** "this connects to $X$," "this relates to $X$," "this is reminiscent of $X$," "there is a deep connection between $X$ and $Y$."

**Why it misleads:** "Connects" spans everything from "is a special case of" to "both involve triangles." Without the mechanism, the reader cannot use the connection for anything.

**To discharge:** State what structure is shared, what map realizes the connection, and what the connection buys you. If the answer to the last question is "nothing yet, but it is suggestive," label it as an open observation rather than a result.

### Hedge covering a gap

**Pattern:** *clearly*, *obviously*, *it turns out that*, *one can show*, *it is easy to see*, *of course*, *simply*.

**Why it misleads:** These sometimes mark a genuinely routine step. More often they mark the exact place where the author did not want to write out the argument. The reader cannot distinguish the two cases, and a reader who does not find it obvious is left doubting themselves rather than the text.

**To discharge:** Either write the step out, or say why it is routine ("this follows by the same induction as above"). If the step is genuinely long, state that and point to where it is done.

### Unspecified "the usual" / "the standard"

**Pattern:** "the usual loss of precision," "the standard argument," "with the obvious modifications."

**Why it misleads:** It presumes shared context the reader does not have, and it is unfalsifiable — the reader cannot check a claim whose content is "whatever you were already thinking."

**To discharge:** Say which loss, which argument, which modifications.

> **Needs work** — [Division Property](../theory/Division%20Property.md):
>
> "At the division-property level they become corresponding transformations of the label set, with the usual loss of precision caused by forgetting coefficients."
>
> Two debts, and they are of different kinds — worth separating, because they call for different fixes.
>
> "The usual loss of precision" is a **sequencing** problem, not an unsupported claim. The content appears later in the same document: §Propagation Rules says propagation "keeps every output label that could arise and usually discards coefficient cancellations," and §What the Transform Does Not Replace says "support propagation is an upper bound whenever different exact coefficients can cancel." The reader who finishes the document has been told. The reader at this sentence has not. A forward pointer or a single clause in place discharges it — no new content is needed.
>
> "Corresponding transformations" is the harder debt: which transformations $\mathsf{A}_m$ and $\mathsf{D}_m$ induce on the label set is not stated here or anywhere else in the document. That one needs content.
>
> The discipline this example teaches: check whether a debt is discharged *elsewhere* before deciding it is unsupported. The two diagnoses lead to different work, and reporting the first as the second wastes the author's time.

### Rhetorical "naturally" / "canonically"

**Pattern:** "$X$ is naturally a $Y$," "the canonical choice is $\ldots$," "this is the natural way to $\ldots$."

**Why it misleads:** *Natural* and *canonical* have precise technical meanings (independent of arbitrary choices; natural in the category-theoretic sense). Used loosely they mean "this is the choice I made," which is fine to say directly but should not be dressed as a theorem.

**To discharge:** This is the pattern where **the word is often correct and the fix is not to remove it**. If there genuinely is one obvious construction, keep the word and *write the construction down*. If the construction depends on a choice, name the choice and say why it was made.

> **Needs work** — [Division Property](../theory/Division%20Property.md):
>
> "This definition explains why the division property is naturally a poset theory."
>
> This is a small debt, and the example is worth studying for how small. Most of the supporting content is already adjacent: the same sentence's continuation names the order ("not ordered by ordinary integer size but by coordinatewise containment of masks"), and the preceding paragraph supplies monotonicity ("Larger labels are stronger statements because fewer moments are exempt from vanishing"). What remains undischarged is the phrase "a poset theory," which names a category of thing rather than a property of this thing — and "naturally," which asserts inevitability without an argument for it. The fix is to say what the order is *on* (the labels) and drop the appeal, not to add a paragraph.

### Forward reference substituting for definition

**Pattern:** Using a term with a link to another document as the only explanation.

**Why it misleads:** It breaks the reader's flow and, in aggregate, makes documents unreadable in isolation. A reader three links deep has lost the thread of the original argument.

**To discharge:** One-sentence inline definition, then the link for depth. The link is for the reader who wants the full treatment, not for the reader who needs to understand this sentence.

An internal forward reference to a later section of the same document is different and is usually fine, provided it is specific: "we prove this below via the evaluation matrix" tells the reader exactly where to look and that the debt is acknowledged.

### Overloaded significance

**Pattern:** "this reveals that," "this shows the true nature of," "at a deeper level," "what is really going on is."

**Why it misleads:** It frames a restatement as a discovery. Frequently the sentence that follows is a reformulation, which is useful, but calling it a revelation inflates it.

**To discharge:** Say what kind of statement it is. "Restating this in terms of $\ldots$ makes $\ldots$ immediate" is honest and equally informative.

## Part 3 — Applying These Standards

### For authors

Write the mathematics first and the framing second. Most debts are incurred in summary sentences written to introduce or conclude a section — the technical body is usually fine. When revising, read only the topic sentences and the closing sentences of each section, and apply the core test to each one.

### The agents

Three agents share this work, split so that **the agent that writes is never the agent that judges**.

- **[technical-writer](../../.claude/agents/technical-writer.md)** — the only generative agent. Authors documents, applies guard reports, and self-reviews its own drafts against these standards before handing them back.
- **[docs-theory-guard](../../.claude/agents/docs-theory-guard.md)** — read-only. Reviews prose: flags passages by the Part 2 heuristics, adjudicates whether the claim is true, locates where missing content already exists, and prescribes the fix. Never rewrites.
- **[code-theory-guard](../../.claude/agents/code-theory-guard.md)** — read-only. Reviews code against the theory docs.
- **[docs-steward](../../.claude/agents/docs-steward.md)** — structure, and **arbiter** when a doc and its implementation contradict each other.

```
technical-writer  ⇄  docs-theory-guard  →  author
   composes            diagnoses             decides
   self-reviews        adjudicates truth
                       sources content
```

The loop runs until the guard has no HIGH flags left or the remainder needs you.

### When a doc contradicts the code

Neither guard decides this. Each has a systematic pull toward the artifact it is responsible for — the code guard toward "the doc is stale," the docs guard toward "the code has a bug" — and neither can see its own bias.

`docs-steward` arbitrates, and is impartial precisely because it does not reason about the mathematics: it has no position to defend. It states the contradiction as two incompatible propositions, gathers historical evidence (`git log`/`blame`, tests, commit intent), then dispatches both guards independently on a neutrally framed question — neither one told which artifact is suspect, neither shown the other's answer. The result classifies as:

- both guards say the code is right → documentation defect → `technical-writer`
- both say the doc is right → implementation bug → you
- they disagree → escalated with both arguments intact; a real finding, not a failed arbitration
- both coherent, describing different things → an **intent** question → always you

It never breaks a tie itself. Which side is wrong is frequently undecidable from the artifacts alone, and no agent guesses what was intended.

### The guard's diagnosis

Each flag carries four judgments, which together tell the writer exactly what to do:

1. **Failure mode** — which Part 2 pattern.
2. **Verdict** — TRUE (prose fails to establish a true claim), FALSE (the claim is wrong — the most dangerous case, since precision makes it credible), OVERSTATED, or UNDETERMINED.
3. **Sourcing** — SOURCED (found, cited), DERIVED (proved, every step shown), STANDARD (named result, applied), or UNSOURCED (not available anywhere it looked).
4. **Prescription** — EXPAND, WEAKEN, CUT, CORRECT, or **ASK AUTHOR**.

`ASK AUTHOR` is the class that reaches you. It means the claim may well be true but nothing in the repository or the cited literature establishes it.

### The constraint that makes this trustworthy

**No agent may invent mathematics to close a gap.** The guard resolves only by citation or by derivation with every step exhibited — an unexhibited derivation is an assertion, and being the theory agent does not license assertion. The writer composes only from what the guard sourced, and is forbidden from filling in an `ASK AUTHOR` because a plausible completion occurred to it.

This is why the roles are split. A single agent that both writes and judges will rationalize its own gaps; that is the origin of the failure mode these standards exist to prevent.

### What reaches you

Open questions, marked in the text where they occur:

```
> **OPEN:** <the precise question> — needs author input.
> <what has been established so far, and where>
```

A report weighted toward open questions is a good report if that is the true distribution. Four responses are all legitimate:

- The claim is true and you can supply the justification — the writer states it.
- A weaker version is what you meant — the writer weakens it.
- The claim does not survive a precise statement — cut it.
- The guard's own reasoning is wrong — reject it and say why.

Nothing the agents produce is authoritative. Read proposed text against these standards before accepting it; a fix can discharge one debt while incurring another.

### Scope

These standards apply to mathematical prose: theory docs, architecture docs, docstring theory paragraphs, and inline comments that explain an invariant. They do not apply to parameter descriptions, API listings, changelogs, or ordinary code comments, where terseness is the correct register.

## Related

- [Notation and Terminology](Notation%20and%20Terminology.md) — library-wide vocabulary and coefficient conventions.
- [Polynomial Conventions](Polynomial%20Conventions.md) — primal/dual distinctions, where imprecise identification causes real bugs.
- [Matrix Indexing](Matrix%20Indexing.md) — orientation conventions, likewise.
- [Algebraic Normal Form](../theory/Algebraic%20Normal%20Form.md) — the reference example for the patterns in Part 1.
