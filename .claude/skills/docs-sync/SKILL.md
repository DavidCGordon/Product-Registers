---
name: docs-sync
description: Reconcile documentation to a batch of code changes the user made by hand. Maps the blast radius of a diff, checks whether the changes break invariants documented elsewhere, updates the affected docs and docstrings through the writer/guard pipeline, and arbitrates any drift the diff does not explain.
---

# Documentation Sync

Use this when **the user has made code changes manually** and the documentation needs to catch up.

This is not a general audit. It is scoped to one diff and reconciles against it. For a codebase-wide sweep with no particular change in view, use `docs-audit` instead.

## The Governing Assumption

**The user's changes were deliberate.** They wrote this code on purpose, so for anything inside the diff, the code is authoritative and the docs follow. Do not open a `docs-steward` arbitration over a mismatch the diff explains — arbitration exists for when nobody knows which side is right, and here the author does.

That assumption holds only *inside* the diff. It does not license rewriting any doc that happens to disagree with any code.

## Three Kinds of Mismatch

Everything a sweep turns up falls into one of these. Conflating them is the main way this workflow goes wrong.

**1. Explained by the diff** — the doc describes the old behavior, the user changed the behavior.
→ Update the doc. No arbitration. This is the bulk of the work.

**2. Not explained by the diff** — a contradiction in code the user did not touch.
→ **Pre-existing drift.** Do not silently rewrite it under cover of this change. Route it to `docs-steward` for arbitration, or report it and leave it alone. The user did not ask about this code and may not know it is broken.

**3. Inside the diff, but breaks an invariant documented elsewhere** — the change is deliberate, but a *different* doc records a guarantee that the change invalidates.
→ **Surface it before writing anything.** This is the highest-value finding in the whole workflow: the user changed X, some doc says X guarantees Y, and they may not have realized Y was load-bearing. Ask before documenting the new behavior as intended.

Category 3 is why "the user made the change, so update the docs" is not the whole workflow. A deliberate change can have consequences its author did not intend, and the docs are frequently where those consequences are written down.

## Workflow

### 1. Establish the diff

Ask the user what to sync against if it is not obvious: uncommitted working-tree changes, staged changes, a commit range, or a branch. Then get it precisely:

```bash
git status
git diff                    # unstaged
git diff --staged           # staged
git diff <base>..HEAD       # a range of commits
```

Read the actual changes, not just the file list. You need to know *what behavior changed*, not merely what files moved.

### 2. Map the blast radius — `docs-steward`

Spawn `docs-steward` to find everything that references the changed code:

- Theory docs and architecture docs naming the changed symbols
- Docstrings on and around the changed functions
- `docs/Theory_Map.md` entries for the changed modules
- Code examples in docs and skills that call the changed API
- Signatures quoted anywhere in `docs/`, `.claude/`, or `CLAUDE.md`

It fixes the purely mechanical drift itself — moved paths, changed signatures, stale imports. It reports what needs mathematical judgment rather than guessing at it.

### 3. Check for broken invariants — `code-theory-guard`

Spawn `code-theory-guard` on the diff. The question is not "is this code good" but:

> Does this change preserve the invariants the theory docs claim about these modules? Where it does not, name the invariant, the doc that states it, and the failure case.

This is what produces category 3. Give it the blast-radius list from step 2 so it knows which docs make claims about this code.

### 4. Triage

Sort every finding into the three categories above. Then, **before writing any prose**, report category 3 to the user:

> `Fibonacci.period()` now returns the least period rather than a multiple. `docs/theory/X.md:44` states that callers may rely on the returned value being a multiple of the register length — that guarantee no longer holds. Intended?

Do not proceed to document new behavior as intended until the user confirms it. Category 2 findings go in the same report, marked as pre-existing and not part of this change.

### 5. Update the docs — `technical-writer`

Spawn `technical-writer` with the category 1 findings and the confirmed category 3 answers. Give it:

- The diff, so it knows what actually changed
- The blast-radius list, so it knows every surface to update
- The invariant findings, so it does not restate a guarantee that no longer holds

It updates docstrings, theory docs, and architecture docs, and self-reviews before returning. It marks anything it cannot source as `**OPEN:**` rather than guessing.

### 6. Review the prose — `docs-theory-guard`

Spawn `docs-theory-guard` on what the writer produced. Loop steps 5–6 until no HIGH flags remain. `ASK AUTHOR` flags go to the user.

Skip this only for purely mechanical updates — a renamed parameter in a `:param:` block needs no theory review. Anything that restates mathematics does.

### 7. Final structural pass — `docs-steward`

Links resolve, `Theory_Map.md` covers any new or moved modules, examples still run against the current API, nothing was orphaned.

## Scaling Down

The full seven steps are for a substantial diff touching mathematical code. Cut it short when the change is smaller:

- **Renames, signature changes, moved files, no behavior change** → steps 1, 2, 7. No theory review needed; nothing mathematical changed.
- **Behavior changed in one module, docs are thin** → steps 1–5, skipping the guard loop if the writer only updated docstrings mechanically.
- **A single function's semantics changed** → the whole pipeline, because that is exactly where category 3 hides.

Judge by whether the *mathematics* changed, not by how many lines moved.

## Report

- **Synced:** docs and docstrings updated, and which change each reflects.
- **Invariants affected:** category 3 findings and what the user decided about each.
- **Pre-existing drift:** category 2 findings — reported, not fixed, unless the user asked.
- **Open:** `ASK AUTHOR` flags and `**OPEN:**` markers left in the text.
- **Verified:** what was checked and found already in sync.

## Rules

- **Never change code.** This workflow moves docs toward code, never the reverse. If a guard concludes the code has a bug, report it — that is a finding for the user, not a fix to apply.
- **Never document a bug as intended behavior.** If the change looks accidental, ask. The whole point of category 3 is catching this.
- **Do not fix pre-existing drift silently.** Report it and let the user decide whether it is in scope.
- **Do not invent mathematics.** The pipeline's standing constraint: every claim sourced or derived, and an honest `**OPEN:**` marker where neither is possible.
