"""Generated code and equations must compute the function they were generated from.

`test_gate_encodings.py` checks each gate on its own. These tests check whole
functions, where the generators' shared machinery matters: splitting a DAG into
one assignment per subfunction in dependency order, referring to an assigned
subfunction by name, honouring overrides, and (for LaTeX) placing parentheses
by precedence instead of around every gate. Each backend's output is executed
or evaluated and compared against `eval` on every input.

No C compiler or VHDL simulator is assumed. C is checked by rewriting it to
Python -- the generators only ever emit `!` as `(!(...))`, so it maps onto
`(not (...))` without precedence changes -- and additionally compiled and run
when a C compiler is on the PATH (it is on CI). VHDL is checked by an evaluator
that implements the VHDL rules for this subset and rejects what VHDL rejects.
"""
import itertools
import random
import re
import shutil
import subprocess

import pytest

from PyPR.BooleanLogic import AND, CONST, NOT, OR, VAR, XOR
from PyPR.BooleanLogic.Gates import NAND, NOR, XNOR
from PyPR.BooleanLogic.Latex import ATOM, PRODUCT, SUM, LatexOperator, LatexStyle
from PyPR.BooleanLogic.SAT import TseytinFuse

from PyPR.FeedbackFunctions import FeedbackFunction, Fibonacci

NVARS = 4
x = [VAR(i) for i in range(NVARS)]
shared = AND(x[0], x[1])
majority = TseytinFuse(OR(AND(VAR(0), VAR(1)), AND(VAR(1), VAR(2)), AND(VAR(0), VAR(2))))

FUNCTIONS = {
    "xor chain": XOR(x[0], x[1], x[2]),
    "anf": XOR(AND(x[0], x[1]), AND(x[1], x[2], x[3]), x[3], CONST(1)),
    "negated gates, 3 args": XOR(XNOR(x[0], x[1], x[2]), NAND(x[1], x[2], x[3]), NOR(x[0], x[3], x[2])),
    "negated gates, 4 args": AND(XNOR(*x), OR(NAND(*x), NOR(*x))),
    "single-argument gates": XOR(XNOR(x[0]), NAND(x[1]), NOR(x[2]), AND(x[3]), OR(x[0]), NOT(x[1])),
    "shared subfunction": OR(XOR(shared, x[2]), AND(shared, x[3]), shared),
    "argument repeated": XOR(shared, shared, x[0]),
    "mixed nesting": AND(OR(XOR(x[0], x[1]), x[2]), XOR(OR(x[1], x[3]), x[0]), NOT(AND(x[2], x[3]))),
    "nested negation": NOT(NOT(XOR(NOT(x[0]), NAND(x[1], NOT(x[2]))))),
    "fused": XOR(majority(x[0], x[1], x[2]), majority(x[1], x[2], shared)),
    "constant root": CONST(1),
    "variable root": VAR(2),
}

# Seeded random DAGs on top of the hand-written cases. Drawing arguments from a
# growing pool gives shared nodes, repeated arguments and every gate at several
# arities; the Random instance is local, so the functions do not depend on
# test order.
_rng = random.Random(2024)
for _case in range(40):
    _pool = [*x, VAR(_rng.randrange(NVARS)), CONST(0), CONST(1)]
    for _ in range(_rng.randrange(2, 20)):
        _gate = _rng.choice([XOR, AND, OR, XNOR, NAND, NOR, NOT])
        _arity = 1 if _gate is NOT else _rng.choice([1, 2, 2, 3, 4])
        _pool.append(_gate(*(_rng.choice(_pool) for _ in range(_arity))))
    FUNCTIONS[f"random {_case}"] = XOR(*_pool[-3:])

ALL_INPUTS = list(itertools.product([0, 1], repeat=NVARS))


# ── Python ───────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name", FUNCTIONS)
def test_generated_python_agrees_with_eval(name):
    """Executed top to bottom, the lines only work if every subfunction is
    assigned before it is used, so this also checks the dependency order."""
    fn = FUNCTIONS[name]
    lines = fn.generate_python()
    for bits in ALL_INPUTS:
        namespace = {"array": list(bits)}
        exec("\n".join(lines), namespace)
        assert namespace["output"] == int(fn.eval(list(bits))), f"{name}: differs at {bits}"

def test_python_splits_out_each_subfunction_once():
    fn = FUNCTIONS["shared subfunction"]
    assert fn.generate_python() == [
        "fn_1 = (array[0] & array[1])",
        "output = ((fn_1 ^ array[2]) | (fn_1 & array[3]) | fn_1)",
    ]

def test_python_names_are_configurable():
    lines = FUNCTIONS["shared subfunction"].generate_python(
        output_name="out", subfunction_prefix="tmp", array_name="state"
    )
    assert lines[0] == "tmp_1 = (state[0] & state[1])"
    assert lines[-1].startswith("out = ")

def test_python_override_replaces_a_node_everywhere():
    """An override stands in for its node wherever it occurs, and the node gets
    no assignment of its own -- this is how FeedbackFunction.compile reuses a
    subfunction computed for an earlier bit."""
    fn = FUNCTIONS["shared subfunction"]
    lines = fn.generate_python(overrides={shared: "known"})
    assert lines == ["output = ((known ^ array[2]) | (known & array[3]) | known)"]

    for bits in ALL_INPUTS:
        namespace = {"array": list(bits), "known": int(shared.eval(list(bits)))}
        exec("\n".join(lines), namespace)
        assert namespace["output"] == int(fn.eval(list(bits)))

def test_python_override_of_the_root_emits_nothing():
    fn = FUNCTIONS["anf"]
    assert fn.generate_python(overrides={fn: "already_known"}) == []

def test_overrides_are_not_modified():
    overrides = {shared: "known"}
    FUNCTIONS["shared subfunction"].generate_python(overrides=overrides)
    assert overrides == {shared: "known"}


# ── C ────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name", FUNCTIONS)
def test_generated_c_agrees_with_eval(name):
    fn = FUNCTIONS[name]
    lines = fn.generate_c()
    assert all(line.endswith(";") for line in lines)

    # C's `!` on a 0/1 int is logical negation; every `!` is emitted as `(!(`,
    # so the operand is already parenthesized and `not` binds the same way
    python = "\n".join(line[:-1].replace("(!(", "(not (") for line in lines)
    for bits in ALL_INPUTS:
        namespace = {"array": list(bits)}
        exec(python, namespace)
        assert int(namespace["output"]) == int(fn.eval(list(bits))), f"{name}: differs at {bits}"

@pytest.mark.skipif(shutil.which("cc") is None, reason="no C compiler on the PATH")
def test_generated_c_compiles_and_agrees_with_eval(tmp_path):
    """Compile every function into one program which prints each function's
    truth table, so a real compiler checks syntax and semantics together."""
    bodies = []
    for fn in FUNCTIONS.values():
        lines = fn.generate_c()
        names = [line.split(" = ")[0] for line in lines]
        bodies.append(
            "        { int " + ", ".join(names) + ";\n"
            + "".join(f"          {line}\n" for line in lines)
            + '          printf("%d", output != 0); }\n'
        )
    program = (
        "#include <stdio.h>\n"
        "int main(void) {\n"
        f"    int array[{NVARS}];\n"
        f"    for (int mask = 0; mask < (1 << {NVARS}); mask++) {{\n"
        f"        for (int i = 0; i < {NVARS}; i++) array[i] = (mask >> i) & 1;\n"
        + "".join(bodies)
        + '        printf("\\n");\n'
        "    }\n"
        "    return 0;\n"
        "}\n"
    )
    source, binary = tmp_path / "generated.c", tmp_path / "generated"
    source.write_text(program)
    subprocess.run(["cc", "-std=c99", "-o", str(binary), str(source)], check=True)
    rows = subprocess.run([str(binary)], check=True, capture_output=True, text=True).stdout.split()

    for mask, row in enumerate(rows):
        # the program sets array[i] from bit i of mask
        bits = [(mask >> i) & 1 for i in range(NVARS)]
        expected = "".join(str(int(fn.eval(bits))) for fn in FUNCTIONS.values())
        assert row == expected, f"compiled C differs from eval at {bits}"


# ── VHDL ─────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name", FUNCTIONS)
def test_generated_vhdl_is_legal_and_agrees_with_eval(name):
    """Evaluate the VHDL by VHDL's own rules for logical expressions: an
    unparenthesized sequence may use only one operator, nand and nor may not be
    chained at all, sequences associate left to right, and `not` applies to a
    single primary. Anything VHDL would reject fails the test."""
    fn = FUNCTIONS[name]
    lines = fn.generate_VHDL()
    assert all(line.endswith(";") for line in lines)

    def evaluate(expr, bits, signals):
        tokens = re.findall(r"'[01]'|\d+|[A-Za-z_]\w*|[()]", expr)
        pos = 0

        def take():
            nonlocal pos
            pos += 1
            return tokens[pos - 1]

        def primary():
            token = take()
            if token == "(":
                value = sequence()
                assert take() == ")", f"unbalanced parentheses in {expr}"
                return value
            if token == "NOT":
                return 1 - primary()
            if token.startswith("'"):
                return int(token[1])
            if token == "array":
                assert take() == "("
                index = int(take())
                assert take() == ")"
                return bits[index]
            return signals[token]

        def sequence():
            value = primary()
            operator = None
            while pos < len(tokens) and tokens[pos] in {"AND", "OR", "XOR", "XNOR", "NAND", "NOR"}:
                token = take()
                assert operator in (None, token), f"mixed {operator}/{token} without parentheses: {expr}"
                assert operator is None or token not in ("NAND", "NOR"), f"chained {token}: {expr}"
                operator = token
                rhs = primary()
                value = {
                    "AND": value & rhs, "OR": value | rhs, "XOR": value ^ rhs,
                    "XNOR": 1 - (value ^ rhs), "NAND": 1 - (value & rhs), "NOR": 1 - (value | rhs),
                }[token]
            return value

        value = sequence()
        assert pos == len(tokens), f"trailing tokens in {expr}"
        return value

    for bits in ALL_INPUTS:
        signals = {}
        for line in lines:
            target, expr = line[:-1].split(" <= ")
            signals[target] = evaluate(expr, bits, signals)
        assert signals["output"] == int(fn.eval(list(bits))), f"{name}: differs at {bits}"

@pytest.mark.parametrize(("gate", "expected"), [
    (XNOR, "output <= (NOT(array(0) XOR array(1) XOR array(2)));"),
    (NAND, "output <= (NOT(array(0) AND array(1) AND array(2)));"),
    (NOR,  "output <= (NOT(array(0) OR array(1) OR array(2)));"),
], ids=["XNOR", "NAND", "NOR"])
def test_vhdl_negated_gates_past_two_arguments(gate, expected):
    """VHDL reads `a xnor b xnor c` as `(a xnor b) xnor c`, which is a xor b xor
    c -- not the negation of it -- and rejects chained nand/nor outright. Both
    used to be emitted; past two arguments the negated form is written instead."""
    assert gate(x[0], x[1], x[2]).generate_VHDL() == [expected]

def test_vhdl_negated_gates_keep_the_native_operator_for_two_arguments():
    assert XNOR(x[0], x[1]).generate_VHDL() == ["output <= (array(0) XNOR array(1));"]

def test_write_vhdl_declares_every_signal_it_assigns(tmp_path):
    """A subfunction split out of a bit becomes a signal of its own. Those
    signals used to be assigned without being declared, so any register with a
    shared subexpression produced VHDL that would not compile."""
    # bit 0 uses `shared` twice, so it is split out as signal fn_0_1
    F = FeedbackFunction([XOR(shared, x[2], shared), AND(shared, x[3]), x[1], x[0]])
    path = tmp_path / "fpr.vhd"
    F.write_VHDL(str(path))
    vhdl = path.read_text()

    declared = set()
    for line in vhdl.splitlines():
        if line.strip().startswith("signal "):
            names = line.strip()[len("signal "):].split(":")[0]
            declared |= {name.strip() for name in names.split(",")}
    assigned = {
        line.strip().split(" <= ")[0] for line in vhdl.splitlines()
        if " <= " in line and not line.strip().startswith(("curr_state", "output", "next_state("))
    }
    assert assigned == {"fn_0_1"}
    assert assigned <= declared
    assert "output <= curr_state;" in vhdl


# ── LaTeX ────────────────────────────────────────────────────────────────────

# A style whose "LaTeX" is Python with the same precedences: & binds tighter
# than ^ and |, and the style never lets ^ and | meet unparenthesized, so
# Python's own ordering of those two never matters. Executing the lines checks
# the precedence engine -- every parenthesis it omits must be one the
# precedences make unnecessary.
PYTHON_STYLE = LatexStyle(
    operators={
        "XOR": LatexOperator(" ^ ", SUM),
        "AND": LatexOperator(" & ", PRODUCT),
        "OR": LatexOperator(" | ", SUM),
        "XNOR": LatexOperator(" ^ ", SUM, negated=True),
        "NAND": LatexOperator(" & ", PRODUCT, negated=True),
        "NOR": LatexOperator(" | ", SUM, negated=True),
        "NOT": LatexOperator("", ATOM, negated=True),
    },
    variable="v[$index]",
    negation="(1 - ($expr))",
    parentheses=("(", ")"),
    line="$name = $expr",
)

@pytest.mark.parametrize("inline", [False, True], ids=["subfunctions", "inline"])
@pytest.mark.parametrize("name", FUNCTIONS)
def test_latex_precedence_preserves_meaning(name, inline):
    fn = FUNCTIONS[name]
    lines = fn.generate_latex(style=PYTHON_STYLE, subfunction_name="g$index", inline_subfunctions=inline)
    if inline:
        assert len(lines) == 1
    for bits in ALL_INPUTS:
        namespace = {"v": list(bits)}
        exec("\n".join(lines), namespace)
        assert int(namespace["f"]) == int(fn.eval(list(bits))), f"{name}: differs at {bits}"

def test_latex_default_anf_notation():
    """Products by juxtaposition bind tighter than the sum, so an ANF needs no
    parentheses at all."""
    assert FUNCTIONS["anf"].generate_latex() == [
        r"f &= x_{0} x_{1} \oplus x_{1} x_{2} x_{3} \oplus x_{3} \oplus 1"
    ]

def test_latex_flattens_nested_associative_gates():
    fn = XOR(x[0], XOR(x[1], XOR(x[2], x[3])), AND(x[0], AND(x[1], x[2])))
    assert fn.generate_latex() == [r"f &= x_{0} \oplus x_{1} \oplus x_{2} \oplus x_{3} \oplus x_{0} x_{1} x_{2}"]

def test_latex_parenthesizes_where_precedence_requires():
    # a sum inside a product needs grouping; XOR and OR share a level, so
    # wherever they meet the inner one is grouped rather than left ambiguous
    fn = AND(XOR(x[0], x[1]), OR(x[2], XOR(x[3], x[0])))
    assert fn.generate_latex() == [
        r"f &= \left(x_{0} \oplus x_{1}\right) \left(x_{2} \vee \left(x_{3} \oplus x_{0}\right)\right)"
    ]

def test_latex_negation_forms():
    fn = XOR(NOT(x[0]), XNOR(x[1], x[2]), NAND(x[1], XOR(x[2], x[3])))
    # the overline groups its argument, so nothing inside it is parenthesized
    assert fn.generate_latex() == [
        (r"f &= \overline{x_{0}} \oplus \overline{x_{1} \oplus x_{2}} \oplus "
         r"\overline{x_{1} \left(x_{2} \oplus x_{3}\right)}")
    ]
    # a prefix negation binds tightly, so a compound argument is grouped
    assert fn.generate_latex(style=LatexStyle.logical()) == [
        (r"f &= \neg x_{0} \oplus \neg \left(x_{1} \oplus x_{2}\right) \oplus "
         r"\neg \left(x_{1} \wedge \left(x_{2} \oplus x_{3}\right)\right)")
    ]

def test_latex_single_argument_gates_are_their_argument():
    assert XOR(AND(x[0]), OR(x[1])).generate_latex() == [r"f &= x_{0} \oplus x_{1}"]
    assert XNOR(x[0]).generate_latex() == [r"f &= \overline{x_{0}}"]

def test_latex_splits_subfunctions_and_names_them():
    fn = FUNCTIONS["shared subfunction"]
    assert fn.generate_latex() == [
        r"g_{1} &= x_{0} x_{1}",
        r"f &= \left(g_{1} \oplus x_{2}\right) \vee g_{1} x_{3} \vee g_{1}",
    ]
    assert fn.generate_latex(output_name="y", subfunction_name=r"s^{($index)}")[0] == r"s^{(1)} &= x_{0} x_{1}"
    assert fn.generate_latex(inline_subfunctions=True) == [
        r"f &= \left(x_{0} x_{1} \oplus x_{2}\right) \vee x_{0} x_{1} x_{3} \vee x_{0} x_{1}"
    ]

def test_latex_override_is_an_atom():
    fn = AND(XOR(x[0], x[1]), x[2])
    assert fn.generate_latex(overrides={fn.args[0]: r"\sigma"}) == [r"f &= \sigma x_{2}"]

def test_latex_style_fields_are_independent():
    fn = XOR(AND(x[0], NOT(x[1])), CONST(1))
    style = LatexStyle(
        variable="s_{$index}^{(t)}",
        constant=lambda value: r"\mathbf{1}" if value else r"\mathbf{0}",
        parentheses=("(", ")"),
        line="$name = $expr",
    )
    style.operators["XOR"] = LatexOperator(" + ", SUM)
    style.operators["AND"] = LatexOperator(r" \cdot ", PRODUCT)
    assert fn.generate_latex(style=style) == [r"f = s_{0}^{(t)} \cdot \overline{s_{1}^{(t)}} + \mathbf{1}"]

    # the default style is unaffected by edits to another instance's operators
    assert XOR(x[0], x[1]).generate_latex() == [r"f &= x_{0} \oplus x_{1}"]

def test_latex_callable_templates():
    names = "abcd"
    style = LatexStyle(variable=lambda index: names[index])
    lines = FUNCTIONS["shared subfunction"].generate_latex(
        style=style, subfunction_name=lambda index: rf"t_{index}"
    )
    assert lines[0] == "t_1 &= a b"

def test_latex_fused_node_renders_its_template():
    fn = XOR(majority(x[0], x[1], x[2]), x[3])
    assert fn.generate_latex() == [
        r"f &= \left(x_{0} x_{1} \vee x_{1} x_{2} \vee x_{0} x_{2}\right) \oplus x_{3}"
    ]

def test_latex_unknown_gate_names_the_fix():
    class MUX(XOR):
        pass
    with pytest.raises(KeyError, match="no operator for gate 'MUX'"):
        MUX(x[0], x[1]).generate_latex()

    style = LatexStyle()
    style.operators["MUX"] = LatexOperator(r" \diamond ", SUM)
    assert MUX(x[0], x[1]).generate_latex(style=style) == [r"f &= x_{0} \diamond x_{1}"]

def test_feedback_function_latex():
    F = Fibonacci(3, [1, 1, 0, 1])
    body = F.generate_latex(environment=None)
    assert F.generate_latex() == "\\begin{align*}\n" + body + "\n\\end{align*}"
    # every bit gets one equation, highest bit first, over state at time t
    lines = F.generate_latex(environment=None).split(" \\\\\n")
    assert len(lines) == 3
    assert lines[1] == "x_{1}[t+1] &= x_{2}[t]"
    assert lines[2] == "x_{0}[t+1] &= x_{1}[t]"

def test_feedback_function_latex_shares_subfunctions_across_bits():
    """A subfunction is defined by the first (highest) bit that splits it out;
    later bits refer to it by that name instead of defining it again."""
    F = FeedbackFunction([AND(shared, x[3]), XOR(shared, x[2], shared), x[0]])
    body = F.generate_latex(environment=None, subfunction_name="g_{$bit.$index}")
    assert body.split(" \\\\\n") == [
        "x_{2}[t+1] &= x_{0}[t]",
        "g_{1.1} &= x_{0}[t] x_{1}[t]",
        r"x_{1}[t+1] &= g_{1.1} \oplus x_{2}[t] \oplus g_{1.1}",
        "x_{0}[t+1] &= g_{1.1} x_{3}[t]",
    ]
