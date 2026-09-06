---
description: Load the relevant theory docs and constraints before editing a module. Invoke with the file path being edited.
---

# Theory Context

Before editing mathematically sensitive code, load the theory docs that govern the module you are about to change.

## Workflow

1. Read `docs/Theory_Map.md` to find which theory docs apply to the file being edited.
2. Read each listed theory doc.
3. Summarize the key constraints and invariants the edit must respect.
4. Proceed with the edit, keeping those constraints in mind.

## When to Invoke

Before any edit to a file under:
- `Cryptanalysis/` — always
- `Tools/` — always
- `FeedbackFunctions/` — always for algorithm changes; skip for mechanical edits (formatting, imports)

Not needed for:
- `BooleanLogic/` — self-contained, no external theory dependencies
- `FeedbackRegister.py` — simulation harness, not mathematically sensitive
- `JSON_Serialization.py` — mechanical serialization code
- Pure formatting, import reordering, or docstring style changes anywhere

## How to Use the Theory Map

The Theory Map uses source paths relative to `src/PyPR/`. Match the file you're editing against the most specific entry. If editing `Cryptanalysis/Components/EquationSolving/LU_Solver.py`, match against:
1. `Cryptanalysis/Components/EquationSolving/*` (most specific)
2. `Cryptanalysis/Components/*` (fallback)

Read **all** listed docs for the matched entry. If multiple entries match, read the union of all listed docs.

## After Loading

Once you have read the theory docs, summarize in your response:
- The key mathematical invariants that apply to this module.
- Any compatibility constraints (e.g., which stores this solver accepts natively vs. via conversion).
- Any convention choices that affect correctness (polynomial orientation, coefficient ordering, etc.).

Then proceed with the edit.
