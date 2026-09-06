# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important: Approach to Editing This Codebase

Much of this code is highly theoretical and implements non-trivial mathematics from academic papers — algebraic attacks, finite field arithmetic, annihilator theory, resolvent polynomials, root counting via Jordan partitions, etc. The correctness of many algorithms depends on subtle invariants that are not obvious from reading the code alone.

**Before editing a file, consult `docs/Theory_Map.md` for required reading.** Load and read those theory docs before making changes. Use the `/theory-context` skill to automate this lookup.

### Safe Edits (no domain expertise required)

These changes are mechanical and can be made confidently:
- Docstring formatting (converting to split Sphinx style)
- Import reordering and cleanup
- Type annotations on function signatures
- Renaming local variables for clarity
- `__init__.py` re-exports
- Serialization boilerplate (`to_JSON` / `from_JSON` / `generate_ids`)
- Test infrastructure (fixtures, parametrization, file I/O scaffolding)
- Adding `verbose` / `print_depth` plumbing through existing call chains

### Domain-Sensitive Edits (ask first or use the code-theory-guard)

These changes require understanding the underlying mathematics. If uncertain, spawn the `code-theory-guard` agent for an independent assessment, or ask the user:
- Algorithm logic in any module under `Cryptanalysis/`, `Tools/`, or `FeedbackFunctions/`
- Equation generation strategy or monomial indexing
- Polynomial arithmetic (GF(2) operations, coefficient ordering, primal/dual conventions)
- Annihilator computation or low-degree pair selection
- Period calculation or root-counting logic
- Store/solver compatibility wiring (see `docs/architecture/Components_Architecture.md`)
- Any change that alters what a function computes, not just how it's structured

The user is the domain expert — lean on them for mathematical questions.

## Overview

PyPR is a Python research library for simulating and cryptanalyzing feedback registers (LFSRs, NFSRs, CMPRs, FCSRs, MPRs). It provides register simulation, boolean function algebra, period analysis, LFSR synthesis, and multiple algebraic attack implementations.

## Setup & Development

```bash
pip install -e .          # editable install from project root
```

Requires Python 3.12+. Key dependencies: `numba` (JIT), `numpy`, `galois` (finite fields), `python-sat` (SAT solving).

## Running Tests

There is no formal test framework. Tests are ad-hoc scripts in `test/` mirroring the `src/` structure. Run them directly:

```bash
python test/script.py
python test/Cryptanalysis/...
```

## Architecture

### Core Classes

**`FeedbackRegister`** (`src/PyPR/FeedbackRegister.py`) — the main simulation class. Takes a register width and a `FeedbackFunction`, then provides:
- `.seed(state)` — set initial state
- `.clock()` / `.run(limit)` — advance one or many cycles
- `.period()` — detect period (four variants: compiled/uncompiled × safe/unsafe, using Brent's algorithm)
- `.to_JSON()` / `.from_JSON()` / `.to_file()` / `.from_file()` — serialization

**`FeedbackFunction`** (`src/PyPR/FeedbackFunctions/FeedbackFunction.py`) — abstract base class for state update logic. Concrete types:
- `MPR` / `CMPR` — Mersenne Product Register and Composite variant (GF(2^n) operations)
- `Fibonacci` / `Galois` — LFSR variants
- `FCSR` — Feedback Carry Shift Register
- `CrossJoin` / `TFunction` — custom compositions

All feedback functions hold a `fn_list: list[BooleanFunction]` and support `.compile()` for numba JIT optimization.

**`BooleanFunction`** (`src/PyPR/BooleanLogic/BooleanFunction.py`) — DAG-based boolean function representation. Supports:
- `.eval(state)` — evaluate on a bit-state
- `.translate_ANF()` — convert to Algebraic Normal Form (`BooleanANF`)
- Code generation: C, VHDL, Python, LaTeX
- Composition via gate nodes (`XOR`, `AND`, `OR`, `NOT` in `Gates.py`)
- SAT solving integration (`SAT.py`)

### Module Map

```
src/PyPR/
├── FeedbackRegister.py          # main simulation class
├── FeedbackFunctions/           # register update logic (MPR, CMPR, FCSR, etc.)
├── BooleanLogic/                # boolean function DAG, ANF, gates, SAT
│   └── ChainingGeneration/      # templates for composing boolean functions
├── Cryptanalysis/
│   ├── Attacks/                 # cube, algebraic (fast/naive/reduced) attacks
│   └── Components/
│       ├── Annihilators/        # annihilator computation
│       ├── EquationGenerators/  # generate boolean equation systems
│       ├── EquationSolving/     # Gauss elimination, LU, Groebner solvers, guess-and-prune
│       └── EquationStores/      # manage equation systems (LU, Symbolic)
├── Tools/
│   ├── MersenneTools.py         # period analysis utilities
│   ├── ResolventSolving.py      # resolvent polynomial computations
│   ├── RegisterSynthesis/       # Berlekamp-Massey, LFSR/NLFSR/FCSR synthesis
│   ├── AlgClosure/              # discrete log, number theory
│   └── RootCounting/            # root enumeration (Jordan partition, etc.)
└── JSON_Serialization.py        # custom serialization framework
```

Documentation lives in `docs/` at the project root, organized into three categories:
- **`docs/theory/`** — mathematical foundations (polynomial duality, root expressions, monomial profiles, algebraic attacks, etc.)
- **`docs/architecture/`** — implementation details and component interactions (Components Architecture, Attack Compatibilities)
- **`docs/conventions/`** — style, terminology, patterns (Notation and Terminology)

See `docs/Theory_Map.md` for the full index mapping source modules to required reading.

### Performance Pattern

Critical paths use numba `@njit`. All feedback functions and register clocking have both compiled and uncompiled modes. Call `.compile()` on a `FeedbackFunction` before creating a `FeedbackRegister` to enable the compiled path. Pointer swapping via XOR avoids array copies in the hot loop.

### Serialization Pattern

Custom JSON framework throughout: classes implement `generate_ids()`, `_generate_JSON_entry()`, and `_parse_JSON_entry()` to support complex object graphs. Top-level API: `to_JSON()` / `from_JSON()`, `to_file()` / `from_file()`.

### EquationSolving Pattern

Solvers in `Cryptanalysis/Components/EquationSolving/` follow a function + thin-wrapper-class convention (module-level `solve()` does the work, a class stores config and forwards to it) so attacks can swap solvers polymorphically. See the package docstring in `Cryptanalysis/Components/EquationSolving/__init__.py` for the full rationale before adding a new solver.

### Components Architecture (Stores / Adapters / Annihilators / Generators / Solvers)

`Cryptanalysis/Components/` has five subdirectories that interoperate in specific, non-obvious ways — which store types each solver natively accepts vs. only reaches via conversion, what an annihilator is and how it's (manually) wired into RAA/FAA, and the one hard incompatibility (NAA rejects `GrobnerSolver`). See `docs/architecture/Components_Architecture.md` before making changes that cross these subdirectories, and `docs/architecture/Attack_Compatibilities.md` for the attack-by-attack (NAA/RAA/FAA) store+solver compatibility reference.

### Writing Style for Mathematical Prose

`docs/conventions/Writing Standards.md` is the specification for theory docs, architecture docs, and the theory paragraphs of docstrings. The governing rule: **every sentence asserting a mathematical relationship incurs a debt**, discharged by stating the relationship precisely enough that the reader could verify or reconstruct it without an outside source.

The test to apply to any framing sentence: *if I deleted this, would the reader lose information, or only lose the impression that something was explained?*

This is not a vocabulary restriction. Words like *natural*, *canonical*, and *analogous* are correct in many places — the requirement is that the construction, correspondence, or map they refer to is written down nearby. Read the standards doc before writing theory content; the `docs-theory-guard` agent enforces it, and `technical-writer` writes to it.

### Docstring Style

Sphinx format (per `.vscode/settings.json`, which is local-only and not checked in — the convention below is the authoritative statement of it). Always use the **split style**:

```
:param param_name: description
:type param_name: TypeAnnotation
:return: description
:rtype: TypeAnnotation
:raises ExceptionType: condition
```

NOT the inline style (`:param Type name: description`). Structure: short summary line, blank line, optional theory paragraph(s), blank line, then `:param`/`:type`/`:return`/`:rtype`/`:raises` blocks. Type annotations belong on the function signature itself, not only in the docstring.

Frame docstrings in mathematical/theoretical language — refer to the underlying algebra, finite fields, or register theory rather than describing implementation steps. For `RootExpression` and `MonomialProfile`, lead with the linear-complexity / Binet-analogy framing (what users care about), not the Jordan-block implementation (how it's computed).

### Test Style

No single-use helpers. Inline helper functions rather than extracting them unless they are used in 5+ places. Tests should be readable without scrolling to a helper definition. Code that depends on non-obvious theory (polynomial conventions, field arithmetic, Jordan partitions) must carry inline comments explaining the invariant being tested.

## Agents and Skills

### Scope: CLAUDE.md vs. Agent Files

`CLAUDE.md` holds global workflow patterns, project-wide conventions, and the module map — things every session needs regardless of which agent is running. Agent-specific behavior (writing standards, output formats, what an agent does and doesn't do) belongs in the agent's own file under `.claude/agents/`. If a guideline applies only when a specific agent is active, put it there, not here.

### Custom Agents (`.claude/agents/`)

Two **guards** review and never write. One **writer** writes and never signs off on itself. Two support agents handle structure and experiments.

- **`code-theory-guard`** — Read-only. Validates that a change to mathematical code preserves correctness and stays aligned with the theory docs. Reports SAFE / CONCERN / INCORRECT with the invariant named and the failure case described. Reports code/doc divergence but does not adjudicate it — that goes to `docs-steward`.
- **`docs-theory-guard`** — Read-only. Validates mathematical prose. Flags passages that assert a relationship without establishing it, **adjudicates whether the claim is actually true**, searches for where the missing content already exists, and prescribes the fix. Every flag carries a failure mode, a verdict (TRUE / FALSE / OVERSTATED / UNDETERMINED), a sourcing label (SOURCED / DERIVED / STANDARD / UNSOURCED), and a prescription (EXPAND / WEAKEN / CUT / CORRECT / ASK AUTHOR). Enforces `docs/conventions/Writing Standards.md`. Never rewrites.
- **`technical-writer`** — The only generative agent. Authors theory docs, architecture docs, and docstring theory paragraphs; applies guard reports; and self-reviews its own drafts against the standards before handing them back. Never signs off on its own work — always returns to `docs-theory-guard`.
- **`docs-steward`** — Documentation structure **and drift arbiter**. Links resolve, documented signatures match the code, `Theory_Map.md` is current, cross-references hold. Tracks coverage and liveness: orphaned docs, uncovered modules, dangling references. When code and docs contradict each other, it detects the contradiction, gathers historical evidence (`git log`/`blame`, tests, commit intent), and dispatches **both guards independently on a neutrally framed question** to determine which side is wrong. Fixes mechanical drift directly. Does **not** judge writing quality or mathematics — including never breaking a tie between the guards.
- **`experimenter`** — Designs and runs probe scripts when a claim needs empirical verification. The guards do not run code; they name this agent instead.

### Recommended Workflow for Domain-Sensitive Changes

1. Use `/theory-context` to load the relevant theory docs before editing.
2. Make the change.
3. If uncertain about correctness, spawn `code-theory-guard` for an independent assessment.
4. If it reports that the code and a theory doc contradict each other, hand that to `docs-steward` to arbitrate — not to `technical-writer` to "fix the doc." Which side is wrong is not yet known.
5. After the change is finalized, spawn `docs-steward` to verify documentation structure stays in sync.

### Recommended Workflow for Syncing Docs to Manual Code Changes

Use the `docs-sync` skill. It is scoped to one diff and treats the user's changes as deliberate, so the code is authoritative *within the diff* and no arbitration is opened for mismatches the diff explains.

The distinction that matters: a sweep turns up three kinds of mismatch, and they get different treatment. Ones **explained by the diff** are ordinary doc updates. Ones **not explained by the diff** are pre-existing drift in code the user did not touch — report or arbitrate, never rewrite silently. Ones **inside the diff that break an invariant documented elsewhere** get surfaced to the user before anything is written, because a deliberate change can have consequences its author did not intend, and the docs are often where those consequences were recorded.

### Recommended Workflow for Writing Documentation

1. Read `docs/conventions/Writing Standards.md`.
2. Spawn `technical-writer` to author or revise the document. It self-reviews before returning.
3. Spawn `docs-theory-guard` to review the result.
4. Loop 2–3 until no HIGH flags remain. The writer applies EXPAND / WEAKEN / CUT / CORRECT itself; only `ASK AUTHOR` flags need you.
5. Answer the open questions. They are marked `**OPEN:**` in the text where they occur.
6. Spawn `docs-steward` to wire the doc into `Theory_Map.md` and verify links resolve.

### Why the Roles Are Split This Way

Each agent answers exactly one question:

- `docs-steward` — *is this connected and current, and where artifacts disagree, what process settles it?* Mechanical, verifiable.
- `docs-theory-guard` — *is this claim stated precisely enough to be checked, and is it true?*
- `code-theory-guard` — *does this code implement what the theory docs claim?*
- `technical-writer` — *how should this be said?*

**The agent that writes is never the agent that judges.** A single agent doing both will rationalize its own gaps — that is the origin of the failure mode the writing standards exist to prevent. Bundling structure with writing quality had the same effect earlier: the standards got treated as a checklist item rather than the core task.

**Neither guard adjudicates a code/doc contradiction.** Each has a systematic pull toward its own artifact — the code guard toward "the doc is stale," the docs guard toward "the code has a bug" — and neither can see its own bias. `docs-steward` arbitrates precisely *because* it does not reason about the mathematics: it has no position to defend. It dispatches both guards independently, on a question framed neutrally and without showing either the other's answer, then classifies the result as a documentation defect, an implementation bug, a genuine disagreement, or an intent question. It never breaks a tie itself.

**Which side is wrong is often undecidable from the artifacts.** Code and docs can each be internally coherent while describing different intentions. That is an intent question, and it always goes to the user — no agent guesses what was meant.

**No agent may invent mathematics to close a gap.** The guard resolves only by citation or by derivation with every step exhibited; an unexhibited derivation is an assertion, and being the theory agent does not license assertion. The writer composes only from what the guard sourced, and may not fill in an `ASK AUTHOR` because a plausible completion occurred to it. When nothing establishes a claim, the honest output is an open question — and a report weighted toward open questions is a good report if that is the true distribution.

**Agents hand off rather than stretch.** Each agent file names the specialist for work outside its scope, and reports that handoff instead of attempting the work. Agents cannot spawn each other; the main session reads the recommendation and dispatches.

### Skills (`.claude/commands/`)

- **`/pypr-api`** — Comprehensive API reference for all PyPR classes and functions.
- **`/pypr-probe`** — Guide for writing and running quick experiment scripts.
- **`/theory-context`** — Load relevant theory docs before editing a module.
- **`/attack-setup`** — Guided workflow for constructing attack experiments with correct component wiring.
- **`docs-audit`** (`.claude/skills/docs-audit/SKILL.md`) — Scan code, docs, docstrings, tests, and Git history for drift or corrected mathematical assumptions, then synchronize the relevant documentation and pointers.

## Knowledge Base

Project knowledge lives in `docs/` — not in memory files. When you discover a convention gotcha, a non-obvious invariant, or a mathematical subtlety while working on this codebase:

1. **Add it to the appropriate doc** in `docs/theory/`, `docs/architecture/`, or `docs/conventions/`.
2. **Add inline cross-links** (wiki-style) to related docs where the new content is referenced.
3. **Update `docs/Theory_Map.md`** if the new content adds a required-reading entry for a source module.

Do not create memory files. All persistent project knowledge belongs in the versioned docs.

### Documentation Parity and Organization

- Update the appropriate versioned documentation in the same change whenever an experiment, code inspection, test, or review reveals new information. This is especially important for corrections to mathematical assumptions, polynomial conventions, coefficient ordering, shift direction, field representations, recurrence identities, or algorithmic invariants.
- Keep implementation, docstrings, theory docs, architecture docs, examples, and tests in parity. When a public signature, behavior, convention, or invariant changes, inspect and update every affected surface.
- Do not preserve a documented claim after evidence shows it is false or incomplete. Replace it with the verified relationship and record important limitations or counterexamples.
- Keep mathematical foundations in `docs/theory/`, implementation and component relationships in `docs/architecture/`, and terminology or representation conventions in `docs/conventions/`. Split documents when they become catch-alls or mix independent topics; avoid oversized files.
- Keep `docs/Theory_Map.md` current when source modules gain prerequisites, topics move, or new required-reading documents are added. Add pointers from relevant module and class docstrings to the documents they require.
- Use focused probes or tests to establish mathematical claims, and distinguish verified behavior from intended behavior or untested conjecture.
