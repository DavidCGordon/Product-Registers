"""Notation for rendering BooleanFunction DAGs as LaTeX.

`BooleanFunction.generate_latex` walks a function bottom-up and asks each node to
render itself; everything about *how* the result looks lives in a `LatexStyle`.
A style says which symbol joins each gate's arguments, how tightly each gate
binds, how negation, variables and constants are written, and how an equation
line is laid out. Every field can be replaced independently, so changing the
emitted equations never requires touching the gates themselves::

    style = LatexStyle(variable="s_{$index}")
    style.operators["XOR"] = LatexOperator(" + ", precedence=1)
    fn.generate_latex(style=style)

Template fields are `string.Template` strings, so a placeholder is written
`$name` and LaTeX braces need no escaping (write `$$` for a literal dollar
sign). Anywhere a template is accepted, a callable taking the same fields as
keyword arguments may be given instead, for notation a template cannot express.
"""
from collections.abc import Callable
from dataclasses import dataclass, field
from string import Template
from typing import Any, NamedTuple

# Binding strengths. A higher value binds more tightly; an argument is
# parenthesized when it binds more loosely than the gate that contains it.
SUM = 1       # XOR, OR and their negations
PRODUCT = 2   # AND and its negation
UNARY = 3     # a prefix negation such as \neg
ATOM = 100    # variables, constants, names, and anything already delimited


class LatexTerm(NamedTuple):
    """A rendered sub-expression and what its parent needs to know to place it.

    :ivar text: The LaTeX for the sub-expression.
    :ivar precedence: How tightly the outermost operator of `text` binds.
    :ivar symbol: The join symbol of the outermost operator when it is
        associative, so an argument joined by the same symbol can be flattened
        into its parent without parentheses. None otherwise.
    """
    text: str
    precedence: int = ATOM
    symbol: str | None = None


@dataclass
class LatexOperator:
    """How one gate class is written.

    :ivar join: The LaTeX placed between consecutive arguments, spaces included.
        An empty or whitespace join writes the arguments side by side, which is
        how a product reads in ANF (`x_{0} x_{1}`).
    :ivar precedence: How tightly the gate binds, compared against its
        parent's precedence to decide on parentheses (see `SUM`, `PRODUCT`).
    :ivar negated: Whether the gate is the negation of the joined arguments,
        as XNOR, NAND and NOR are. The negation is written with the style's
        `negation` template.
    :ivar associative: Whether nested uses of the same `join` may be flattened,
        e.g. `a \\oplus (b \\oplus c)` written as `a \\oplus b \\oplus c`.
    """
    join: str
    precedence: int
    negated: bool = False
    associative: bool = True


def _algebraic_operators() -> dict[str, LatexOperator]:
    return {
        "XOR": LatexOperator(r" \oplus ", SUM),
        "AND": LatexOperator(" ", PRODUCT),
        "OR": LatexOperator(r" \vee ", SUM),
        "XNOR": LatexOperator(r" \oplus ", SUM, negated=True),
        "NAND": LatexOperator(" ", PRODUCT, negated=True),
        "NOR": LatexOperator(r" \vee ", SUM, negated=True),
        "NOT": LatexOperator("", ATOM, negated=True),
    }


def fill_template(template: str | Callable[..., str], **fields: Any) -> str:
    """Fill a template: a `string.Template` string, or a callable taking the fields as keywords.

    :param template: The template to fill.
    :type template: str | Callable[..., str]
    :param fields: The placeholder values.
    :type fields: Any
    :return: The filled template.
    :rtype: str
    """
    if callable(template):
        return template(**fields)
    return Template(template).substitute(**{k: str(v) for k, v in fields.items()})


@dataclass
class LatexStyle:
    """The notation `generate_latex` writes a function in.

    The default is the algebraic notation used for ANF: exclusive-or as
    `\\oplus`, conjunction by juxtaposition, disjunction as `\\vee`, and negation
    as an overline. `LatexStyle.logical()` gives the propositional-logic
    alternative. To change one aspect, pass that field (or use
    `dataclasses.replace` on an existing style) and leave the rest.

    XOR and OR are given the same precedence on purpose: there is no settled
    convention for how they bind relative to each other, so wherever the two
    meet the argument is parenthesized rather than left to the reader.

    :ivar operators: The `LatexOperator` for each gate, keyed by class name. A
        gate class with no entry here must override `_generate_latex` itself.
    :ivar variable: Template for a variable, with field `$index`.
    :ivar constant: Template for a constant, with field `$value`.
    :ivar negation: Template for a negated expression, with field `$expr`.
    :ivar negation_delimits: Whether the negation template already groups its
        argument, as `\\overline{$expr}` does. If so, the argument is never
        parenthesized and the result binds like an atom. If not, as for a prefix
        `\\neg $expr`, a compound argument is parenthesized and the result binds
        at `UNARY` precedence.
    :ivar parentheses: The opening and closing delimiters used when an argument
        needs grouping.
    :ivar line: Template for one equation, with fields `$name` and `$expr`. The
        default aligns on the equals sign for an `align` environment.
    """
    operators: dict[str, LatexOperator] = field(default_factory=_algebraic_operators)
    variable: str | Callable[..., str] = "x_{$index}"
    constant: str | Callable[..., str] = "$value"
    negation: str | Callable[..., str] = r"\overline{$expr}"
    negation_delimits: bool = True
    parentheses: tuple[str, str] = (r"\left(", r"\right)")
    line: str | Callable[..., str] = "$name &= $expr"

    @classmethod
    def logical(cls, **changes: Any) -> "LatexStyle":
        """The propositional-logic notation: `\\oplus`, `\\wedge`, `\\vee` and a prefix `\\neg`.

        :param changes: Further fields to set on the returned style.
        :type changes: Any
        :return: A style using logical connectives.
        :rtype: LatexStyle
        """
        operators = {
            "XOR": LatexOperator(r" \oplus ", SUM),
            "AND": LatexOperator(r" \wedge ", PRODUCT),
            "OR": LatexOperator(r" \vee ", SUM),
            "XNOR": LatexOperator(r" \oplus ", SUM, negated=True),
            "NAND": LatexOperator(r" \wedge ", PRODUCT, negated=True),
            "NOR": LatexOperator(r" \vee ", SUM, negated=True),
            "NOT": LatexOperator("", ATOM, negated=True),
        }
        fields = {"operators": operators, "negation": r"\neg $expr", "negation_delimits": False}
        fields.update(changes)
        return cls(**fields)

    # rendering, called by the node hooks:
    def render_variable(self, index: int) -> LatexTerm:
        """Render the variable with the given index as an atom.

        :param index: The variable's index.
        :type index: int
        :return: The rendered variable.
        :rtype: LatexTerm
        """
        return LatexTerm(fill_template(self.variable, index=index))

    def render_constant(self, value: Any) -> LatexTerm:
        """Render a constant as an atom.

        :param value: The constant's value.
        :type value: Any
        :return: The rendered constant.
        :rtype: LatexTerm
        """
        return LatexTerm(fill_template(self.constant, value=value))

    def render_gate(self,
        name: str,
        args: list["LatexTerm | str"]
    ) -> LatexTerm:
        """Render a gate over its already-rendered arguments.

        An argument given as a plain string -- a subfunction's name or an
        override -- is treated as an atom. A gate with one argument is that
        argument (negated, for a negated gate), so `XOR(x)` is written `x`.

        :param name: The gate's class name, looked up in `operators`.
        :type name: str
        :param args: The rendered arguments, in order.
        :type args: list[LatexTerm | str]
        :raises KeyError: If `operators` has no entry for `name`.
        :return: The rendered gate.
        :rtype: LatexTerm
        """
        if name not in self.operators:
            raise KeyError(
                f"LatexStyle has no operator for gate '{name}': add one to "
                f"style.operators, or override _generate_latex on the class"
            )
        op = self.operators[name]
        terms = [LatexTerm(a) if isinstance(a, str) else a for a in args]

        if len(terms) == 1:
            inner = terms[0]
        else:
            inner = LatexTerm(
                op.join.join(self._operand(t, op) for t in terms),
                op.precedence,
                op.join if op.associative else None,
            )

        return self.negate(inner) if op.negated else inner

    def negate(self, term: LatexTerm) -> LatexTerm:
        """Render the negation of an expression.

        :param term: The expression to negate.
        :type term: LatexTerm
        :return: The negated expression.
        :rtype: LatexTerm
        """
        if self.negation_delimits:
            return LatexTerm(fill_template(self.negation, expr=term.text))
        text = term.text if term.precedence >= UNARY else self._group(term.text)
        return LatexTerm(fill_template(self.negation, expr=text), UNARY)

    def format_line(self, name: str, expr: str) -> str:
        """Lay out one equation.

        :param name: The left-hand side.
        :type name: str
        :param expr: The right-hand side.
        :type expr: str
        :return: The equation, laid out by the `line` template.
        :rtype: str
        """
        return fill_template(self.line, name=name, expr=expr)

    def _operand(self, term: LatexTerm, parent: LatexOperator) -> str:
        # flatten a nested use of the same associative operator; otherwise an
        # argument binding no more tightly than its parent gets grouped
        if term.precedence > parent.precedence:
            return term.text
        if parent.associative and term.symbol == parent.join:
            return term.text
        return self._group(term.text)

    def _group(self, text: str) -> str:
        return f"{self.parentheses[0]}{text}{self.parentheses[1]}"


def partial_name(template: str | Callable[..., str], **fields: Any) -> str | Callable[..., str]:
    """Fill some of a naming template's fields, leaving the rest for later.

    `FeedbackFunction.generate_latex` uses this to fix `$bit` in a subfunction
    template before each bit's function fills in `$index`.

    :param template: A `string.Template` string or a callable.
    :type template: str | Callable[..., str]
    :param fields: The placeholder values to fix now.
    :type fields: Any
    :return: The partially filled template, of the same kind as the input.
    :rtype: str | Callable[..., str]
    """
    if callable(template):
        fn = template
        return lambda **rest: fn(**fields, **rest)
    return Template(template).safe_substitute(**{k: str(v) for k, v in fields.items()})

