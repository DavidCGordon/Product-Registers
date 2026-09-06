# Theory Prerequisites Map

Maps each source module to the docs an agent should read before editing.
Skills like `/theory-context` parse this file. Agents can also consult it directly.

Paths are relative to `docs/`.

## Documentation Categories

- **`theory/`** — Pure mathematical foundations: roots and the algebraic closure, polynomial types, root expressions, root multiplicities and Jordan decomposition, decimation/expansion and LC, monomial profiles, division property, resolvent analysis, cube equation generation, algebraic attack theory, trajectory vs column LC
- **`architecture/`** — Library structure and implementation details: component interactions, store/solver/generator compatibility, attack implementation specifics, mesh optimization
- **`conventions/`** — Style, terminology, patterns: notation, matrix indexing, bit ordering, coefficient representation, vocabulary, polynomial conventions, coset region diagrams, writing standards

## Writing Documentation

```
Any theory/ or architecture/ doc     → conventions/Writing Standards.md
Docstring theory paragraphs          → conventions/Writing Standards.md
```

Required reading before writing or substantially revising mathematical prose. Written to by the `technical-writer` agent; enforced by `docs-theory-guard`.

## Cryptanalysis

```
Cryptanalysis/Attacks/*                → theory/Algebraic Attacks.md, architecture/Attack_Compatibilities.md
Cryptanalysis/Components/*             → architecture/Components_Architecture.md
Cryptanalysis/Components/Annihilators/*        → architecture/Components_Architecture.md §3
Cryptanalysis/Components/Adapters/*            → architecture/Components_Architecture.md §2
Cryptanalysis/Components/EquationGenerators/*  → architecture/Components_Architecture.md §4, theory/Cube Equation Generation.md
Cryptanalysis/Components/EquationSolving/*     → architecture/Components_Architecture.md §5
Cryptanalysis/Components/EquationStores/*      → architecture/Components_Architecture.md §1
```

## Tools

```
Tools/RootCounting/*            → theory/Roots and the Algebraic Closure.md, theory/Root Expressions and LC Estimation.md, theory/Root Multiplicities and Jordan Decomposition.md, theory/Monomial Profile Theory.md
Tools/RootCounting/MeshOptimization.py → architecture/Mesh Optimization.md
Tools/ResolventSolving.py       → theory/Resolvent Analysis.md
Tools/RegisterSynthesis/*       → conventions/Polynomial Conventions.md (Berlekamp-Massey section)
Tools/AlgClosure/*              → theory/Roots and the Algebraic Closure.md, theory/Algebraic Closure.md, theory/Decimation Expansion and Linear Complexity.md
Tools/CostEstimation.py         → theory/Cube Equation Generation.md
```

## FeedbackFunctions

```
FeedbackFunctions/Fibonacci.py  → conventions/Polynomial Conventions.md
FeedbackFunctions/Galois.py     → conventions/Polynomial Conventions.md
FeedbackFunctions/MPR.py        → conventions/Polynomial Conventions.md (MPR section)
FeedbackFunctions/CMPR.py       → conventions/Notation and Terminology.md, theory/Root Expressions and LC Estimation.md, theory/Root Multiplicities and Jordan Decomposition.md, architecture/Mesh Optimization.md
FeedbackFunctions/TFunction.py  → conventions/Notation and Terminology.md (T-Function Bit Order section), theory/Root Multiplicities and Jordan Decomposition.md
FeedbackFunctions/CrossJoin.py  → conventions/Notation and Terminology.md, theory/Root Multiplicities and Jordan Decomposition.md
```

## BooleanLogic

```
BooleanLogic/*                  → (no required theory docs — self-contained)
```

## Division Property

```
Division-property analysis       → theory/Division Property.md, theory/Algebraic Normal Form.md, theory/Moebius Inversion on the Bit-Subset Lattice.md
```

## Other

```
Interleaved and clock-controlled constructions → theory/Interleaving and Clock-Controlled Registers.md, theory/Decimation Expansion and Linear Complexity.md, theory/Root Multiplicities and Jordan Decomposition.md
```

```
FeedbackRegister.py             → (no required theory docs — consult /pypr-api for API details)
JSON_Serialization.py           → (no required theory docs — mechanical serialization code)
```
