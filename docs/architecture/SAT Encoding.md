# SAT Encoding

`PyPR.BooleanLogic.SAT` turns a `BooleanFunction` DAG into CNF with the Tseytin
transform, and every SAT-based method (`sat`, `enum_models`,
`functionally_equivalent`) goes through that encoding. This document states the
invariant the encoder is built on, derives why it makes the whole CNF correct,
and describes what it allows — fusion nodes in particular.

## The invariant: every node owns its labels and the clauses implementing it

For every node $N$ the encoder reaches:

1. **Labels.** `node_labels[N]` is a non-empty list of solver variables, and its
   last entry is $N$'s **result wire**. Any earlier entries are $N$'s interior;
   their layout is $N$'s own business.
2. **Clauses.** $N$ contributes exactly the clauses returned by
   `N._tseytin_clauses(label_map)`. They mention only $N$'s own labels, the
   result wires of its arguments (`label_map[arg][-1]`), and the constant
   literal $\pm 1$.
3. **Correctness, per class.** For every assignment of the arguments' result
   wires, $N$'s clauses are satisfiable, and every assignment satisfying them
   puts on $N$'s result wire the value `N.eval` gives on those arguments. The
   interior wires are whatever the clauses allow.

The CNF of a function is the union of its nodes' clauses together with the unit
clause $(1)$, which fixes variable 1 to true. Nothing else emits clauses.

**Why this makes the whole CNF correct.** Fix an assignment of the input
variables. The walks visit nodes in post-order, so a node is reached only after
its arguments. Suppose every argument's result wire already carries its
subfunction's value in every model. By requirement 3, $N$'s clauses can be
satisfied, and only with $N$'s result wire carrying $N$'s value; by requirement
2 they constrain no other node's wires, so satisfying them cannot disturb what
was established beneath. By induction from the leaves, every assignment of the
inputs extends to a model, and in every model every node's result wire carries
that node's value. `sat` and `enum_models` add the root's result wire as a unit
clause, so their models are exactly the input assignments on which the function
is true.

### Leaves

- `VAR(i)` is labelled `[variable_labels[i]]`. Every `VAR` node with the same
  index shares that one wire, because they denote the same input. It emits no
  clauses.
- `CONST(1)` is labelled `[1]` and `CONST(0)` `[-1]`: the reserved variable, or
  its negation, which the unit clause $(1)$ pins. They emit no clauses.
- New wires start at 2. Wire 0 is never used: 0 ends a clause in DIMACS.

### Gates

The default `_tseytin_labels` claims $\max(1, k-1)$ wires for a node with $k$
arguments: the intermediates of a chain of two-input gates, the result last.
Each gate's `_tseytin_clauses` relates consecutive links of that chain.
`test/BooleanLogic/test_gate_encodings.py` checks requirement 3 directly for
every gate and arity from 1 to 5: the models of a gate's clauses, and of its
negation's, must be exactly the inputs `eval` accepts and rejects.

## The two walks

Both walks traverse `args` only, in post-order, and call only the two hooks.
Neither knows how any class encodes itself.

- `tseytin_labels` labels a node once all of its arguments are labelled. A node
  already in `node_labels` is skipped along with everything beneath it. This is
  what *threading* rests on: passing a previous `(clauses, node_labels,
  variable_labels)` back into `tseytin` encodes a second circuit on the same
  wires wherever it shares nodes with the first. The next free wire is one past
  the highest label so far, and never below 2.
- `tseytin_clauses` visits each node once and collects its clauses, keeping one
  copy of each clause in the order first emitted.

## What the invariant allows

- **Reading and constraining any wire.** A node's value in a model is at
  `node_labels[N][-1]`. A caller can add clauses on any node's wires, for
  example to pin a subexpression, and pass the result to a solver directly.
- **Sharing.** A node used in several places, or in several threaded circuits,
  is labelled and encoded once.
- **Local change.** A class can change its encoding without touching any other
  class, and its correctness is checked in isolation. The encoding of one-input
  XNOR, NAND and NOR, and of NAND and NOR over three or more inputs, was once
  wrong. It was fixed class by class, and no other node's encoding moved.
- **Fusion**, below: a node that encodes a whole subtree its own way.

## Fusion nodes

`TseytinFuse(template)(*args)` builds a `_TseytinFuse`, a node denoting
`template.compose({i: args[i]})` for a template over `VAR(0) .. VAR(n-1)`. It
changes the encoding of that subtree, never the function it computes, and it is
never introduced automatically.

It meets the invariant like any other node:

- **Labels.** Its list holds every wire its encoding needs, with the result
  last. The wires interior to the fused domain have no node of their own, so
  they are published only here. A caller can still read or constrain them, but
  cannot key them by node.
- **Clauses.** Its clauses are computed once, at construction, over placeholder
  wires (the *skeleton*). At emission they are renumbered onto its own labels
  and its arguments' result wires. Constant literals pass through unchanged.

For a template that is a single $n$-input AND or OR over `VAR(0) .. VAR(n-1)`
in order, the skeleton is the direct encoding: one wire and $n + 1$ clauses,
where a chain of two-input gates spends $n - 1$ wires and $3(n-1)$ clauses.
Any other template uses its own ordinary encoding, which is always correct.

Requirement 3 is not automatic for a hand-written encoding.
`node.verify()` checks it by a SAT equivalence call between the node and its
expansion (`node.expand()`). It is kept out of emission, since its cost is
unbounded and it would encode the node it is checking.

The template is copied when the node is constructed and again when the node is
copied, and the walks never enter it. A fusion node therefore shares no objects
with the surrounding DAG. Like any other node it is built from arguments that
already exist, so it adds no way for a DAG to become cyclic short of deliberate
mutation.
