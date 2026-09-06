---
name: docs-audit
description: Scan the PyPR codebase and Git history for documentation drift, mathematical corrections, stale docstrings, broken pointers, and missing theory-map entries; update the relevant docs while keeping code, tests, and documentation in parity.
---

# Documentation Audit

Use this skill when documentation may have drifted from the implementation, when a mathematical assumption has been corrected, or when a codebase-wide documentation review is requested.

**To reconcile docs against a specific set of changes**, use `docs-sync` instead — it is scoped to one diff, treats the user's deliberate changes as authoritative, and dispatches the writer/guard agents. This skill is the broad sweep with no particular change in view.

## Scope

Audit the source tree, tests, examples, `docs/`, `.claude/`, `CLAUDE.md`, and repository history. Treat verified behavior and current code as the primary evidence; use history to recover intent, identify regressions, and explain convention changes, not to preserve claims contradicted by current evidence.

## Workflow

1. Read `CLAUDE.md` and `docs/Theory_Map.md`.
2. Inventory source modules, public classes and functions, docstrings, tests, examples, theory documents, architecture documents, convention documents, skills, agents, and cross-references.
3. For each mathematically sensitive module, read every document mapped by `docs/Theory_Map.md` before evaluating its documentation.
4. Inspect Git history when available:
   - `git log --oneline -- <path>` for relevant changes and prior intent.
   - `git blame <path>` for the origin of a claim or convention.
   - `git show <commit>:<path>` to compare earlier implementations or documentation.
   - `git log -S <term> -- <path>` or `git log -G <regex> -- <path>` to find when a formula, API, or terminology changed.
5. Compare code, tests, docstrings, examples, and docs for:
   - changed signatures, imports, paths, or public behavior;
   - mathematical conventions, coefficient ordering, polynomial orientation, shift direction, field bases, recurrence identities, and invariants;
   - claims that experiments or tests disprove;
   - missing or stale links and `docs/Theory_Map.md` entries;
   - documentation files that mix unrelated topics or have become oversized.
6. Use the cheapest focused test, probe, type check, or executable example to distinguish competing interpretations. Record whether each conclusion is verified, historical, intended, or unresolved.
7. Update all affected surfaces in one focused change: implementation docstrings, theory or architecture docs, conventions, examples, tests, agent/skill guidance, and `docs/Theory_Map.md` pointers as applicable.
8. Prefer concise, topic-specific documents in `docs/theory/`, `docs/architecture/`, and `docs/conventions/`. Split catch-all documents when necessary and add cross-links from the index and affected source docstrings.
9. Validate changed examples, links, headings, and focused tests. Re-run the audit for references to renamed or corrected material.

## Editing Rules

- Correct documentation when evidence disproves it; do not silently alter code to match stale prose.
- When history and current behavior disagree, state the disagreement and use tests or experiments to determine the current contract.
- Preserve user changes and avoid unrelated rewrites.
- Keep mathematical explanations explicit about conventions and domains. Include a short derivation or counterexample when correcting a non-obvious claim.
- Keep public code examples self-contained and aligned with current APIs.
- Do not create memory files for repository knowledge; record durable findings in versioned docs.
- Do not rewrite Git history or use destructive Git commands.

## Report

Finish with:

- **Updated:** files changed and the synchronized claims or pointers they contain.
- **Verified:** focused checks, tests, links, and history inspected.
- **Open questions:** unresolved mathematical or historical ambiguity, if any.
