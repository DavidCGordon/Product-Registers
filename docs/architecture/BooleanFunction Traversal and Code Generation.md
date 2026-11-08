# BooleanFunction Traversal and Code Generation

A `BooleanFunction` is a DAG of gate and leaf nodes. This document describes
how whole-function methods on `BooleanFunction` walk that DAG, how the code
generators (`generate_python`, `generate_c`, `generate_VHDL`, `generate_latex`)
share one driver, and what a new node class has to provide. Source:
`src/PyPR/BooleanLogic/BooleanFunction.py`, `Gates.py`, `FunctionInputs.py`,
`Latex.py`, and `FeedbackFunctions/FeedbackFunction.py` for the register-level
writers.

## 1. The walk: `postorder(stop=None)`

Every method that visits a whole function (evaluation, ANF translation,
copying, composition, `merge_redundant`, Tseytin labelling, serialization ids,
code generation, cost estimation) is a loop over `postorder()` that fills a dict
keyed by node:

```python
values = {}
for node in fn.postorder():
    values[node] = node._eval(values, array)
```

The contract the loops rely on:

- **Arguments first.** A node is yielded only after every one of its arguments,
  so `values[arg]` is always present when the node arrives.
- **Each node once.** A node reached along several paths is yielded the first
  time only. That is what makes shared structure cost one evaluation, and what
  makes a copy keep its sharing.
- **Identity, not equality.** Nodes are distinguished by object identity: two
  separate `VAR(0)` objects are two nodes.
- **Order.** Nodes come in the order a left-to-right depth-first search
  finishes them. Code generation and serialization ids depend on this order
  being deterministic, which it is for a given DAG.
- **Iterative.** No recursion, so depth is not limited by Python's recursion
  limit.
- **`stop(node)`** treats a node as already handled: it is neither yielded nor
  descended into, and the caller supplies its value. The library uses it to
  resume from an earlier call (nodes that already have a serialization id, a
  SAT label, or a cost), and to cut the DAG at subfunction and override
  boundaries in code generation.

What differs between methods (the leaf case, the gate case, where results go,
where to stop) stays in the loop body and the `stop` predicate. Nothing is
expressed as a hook on the walk itself.

**Performance.** `postorder` runs under every `eval`, so it allocates nothing
per node and makes no Python call per node when `stop` is `None`; passing a
dict's own `__contains__` as `stop` also avoids a Python-level call. Measured
against a separate traversal loop written inside each method, it is no slower.
A change to `postorder` is a change to the cost of every method listed above.

## 2. Per-node hooks

A gate or leaf class supplies one method per whole-function operation; the
whole-function method calls it with the dict of results for its arguments.

| Hook | Called by | Default |
|------|-----------|---------|
| `_eval(values, array)` | `eval` | none (required) |
| `_eval_ANF(values, array)` | `eval_ANF`, `translate_ANF` | none (required) |
| `_generate_python` / `_generate_c` / `_generate_VHDL(strings, array_name)` | the matching `generate_*` | none (required for that language) |
| `_generate_latex(terms, style)` | `generate_latex` | looks the class name up in `style.operators` |
| `_copy(child_copies)` | `copy`, `compose`, `merge_redundant` | rebuilds `type(self)(*copied_args, arg_limit=...)` |
| `_merge_redundant(cache, subfunctions, in_place)` | `merge_redundant` | drops a one-argument gate; XOR, AND and OR override it to flatten, and XNOR, NAND and NOR to turn a one-argument gate into a NOT |
| `_binarize(cache)` | `binarize` | none (raises `NotImplementedError`) |
| `_tseytin_clauses(label_map)` | `tseytin` | none (raises `NotImplementedError`); see [SAT Encoding](SAT%20Encoding.md) |

**`merge_redundant` and deferred arguments.** XOR, AND and OR share one hook,
`Gates._merge_associative`, which flattens nested gates of the same kind into
one. It does not build its result's arguments when called. It records them in
`cache.pending` (result → the node it was made for, and its arguments as an
insertion-ordered dict), and `merge_redundant` writes every pending result's
`args` after the walk. A same-gate argument that is not a subfunction has only
one parent, so the parent takes that argument's pending terms over instead of
copying them, folding the smaller set into the larger. This makes a nested
chain of n gates cost O(n), where building each intermediate gate in full
costs O(n²). Three consequences follow:

- A hook that reads a child result's `args` must read `cache.pending[result]`
  for a pending child.
- The order of a merged gate's arguments is unspecified.
- A result is taken over only when it was made for that argument itself. A
  one-argument gate passes its child's result through, and that result may
  belong to a node with other parents.

A node class carrying data beyond its arguments (a gate's `arg_limit`, a fused
SAT node's template, a `CONST`'s value) overrides `_copy` to carry that data
over, and lists any field that may hold another serializable object in
`_JSON_references` (a `CONST` whose value is a `BooleanANF` is the example in
the library). Serialization gives an object held in such a field its id before
the node's, and stores the field as that id under the key `__refs__`, so an
object shared between nodes stays shared when loaded back. A field holding a
plain value is stored as it is.

## 3. Code generation: `_generate_assignments`

All four `generate_*` methods are thin wrappers over
`_generate_assignments(hook, context, output_name, subfunction_name, overrides,
inline_subfunctions)`. It returns `(name, rendered value)` pairs, dependencies
first; the wrapper only lays each pair out as a line in its language.

- **Assignments.** Each subfunction (a node with more than one parent, from
  `subfunctions()`) gets its own assignment, named by `subfunction_name(i)` with
  `i` 1-based, in topological order, followed by the function itself under
  `output_name`. With `inline_subfunctions=True` there is a single assignment.
- **Regions.** An assignment's expression is rendered over the nodes reachable
  from its root without passing through another assigned subfunction or an
  override. The walk stops at those boundaries.
- **What a hook sees.** Inside a region, an assigned subfunction stands for its
  name and an override node for its override string, so the hook receives a
  plain string for both and never sees the node behind it.
- **Overrides.** `overrides` maps a node to the text written in its place. An
  override node gets no assignment of its own, even if it is a subfunction or
  the root.

`context` is the hook's second argument: the state array's name for the code
languages, the `LatexStyle` for LaTeX.

### Register level

`FeedbackFunction.write_VHDL` and `FeedbackFunction.generate_latex` call the
per-function generator once per bit, highest bit first, and pass the **same
`overrides` dict** to every call. Each call adds the subfunctions it names, so a
subfunction shared between bits is defined once, in the first bit that uses it,
and later bits refer to it by that name. `write_VHDL` then declares one signal
per assignment.

Only nodes that are subfunctions *within* some bit's function are split out. A
node shared between two bits but used once inside each is not a subfunction of
either, so it is written inline in both.

### Compiled code

`FeedbackFunction.generate_python` is a third register-level writer, and
`compile` wraps its statements in a function. By default it calls
`BooleanFunction.generate_python` per bit with a shared `overrides` dict,
exactly as above. With `optimize=True`, which is what `compile` passes, and an
update the optimizer can improve, it emits the optimized DAG with
`_straight_line_program`, which gives a local to every gate with more than one
user across *all* bits. See [Logic Optimization](Logic%20Optimization.md).

## 4. LaTeX notation: `LatexStyle`

`generate_latex` renders bottom-up into `LatexTerm(text, precedence, symbol)`.
The node decides only *which* operator it is; everything about how it looks
comes from the `LatexStyle` passed in:

- **Operators.** `style.operators` maps a gate class name to
  `LatexOperator(join, precedence, negated=False, associative=True)`. A new gate
  class needs an entry here, not a `_generate_latex` override, unless it is not
  an operator joining its arguments.
- **Parentheses.** An argument is parenthesized when its precedence is lower
  than its parent's. An argument joined by the same associative symbol is
  flattened into its parent instead. XOR and OR share a precedence (`SUM`), so
  where they meet the argument is always parenthesized. No convention for how
  they bind relative to each other is assumed.
- **Negation.** XNOR, NAND, NOR and NOT are written with the `negation`
  template. `negation_delimits` says whether that template already groups its
  argument (`\overline{$expr}`, which then binds like an atom) or not (a prefix
  `\neg $expr`, which parenthesizes a compound argument and binds at `UNARY`).
- **Templates.** The templates (`variable`, `constant`, `negation`, `line`, and
  the `subfunction_name`/`output_name` arguments) are `string.Template` strings,
  so LaTeX braces need no escaping. A callable taking the same fields as
  keywords works anywhere a template does.

`LatexStyle()` is the algebraic notation used for ANF; `LatexStyle.logical()`
is the propositional one.

## 5. Copying

`copy()`, `copy.copy` and `copy.deepcopy` all produce the same deep copy: the
DAG is rebuilt bottom-up through `_copy`, so the copy has exactly the sharing
of the original and shares no node with it. `__deepcopy__` threads `memo`, so a
function deep-copied as part of a larger object (a `FeedbackFunction`'s
`fn_list`, say) keeps any sharing with the other objects copied in the same
call.
