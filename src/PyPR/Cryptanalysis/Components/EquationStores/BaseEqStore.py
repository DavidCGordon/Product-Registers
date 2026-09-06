"""Minimal shared interface for all equation stores.

Every equation store (EqStore, LUEqStore, GroebnerEqStore, SymbolicEqStore)
inherits from this class, ensuring a common ``insert_equation`` contract
and declaring universal store properties.
"""
from __future__ import annotations

from typing import Any

from PyPR.BooleanLogic import BooleanFunction


class BaseEqStore:
    """Minimal shared interface for all equation stores.

    Subclasses define what equation types they accept and how
    they store them internally. The contract includes:

    - ``num_vars``: number of tracked variables/monomials
    - ``num_eqs``: number of stored equations
    - ``insert_equation``: accept an equation and store it

    Universal store properties (set in subclass ``__init__``):

    - ``consistent``: whether the store enforces consistency
      checking (contradictions raise ValueError)
    - ``eager``: whether the store reduces on every insertion
      (True) or defers reduction for batch processing (False).
      Eager stores support per-insertion ``is_determined`` checks.
    - ``filtering``: whether the store discards redundant
      information as equations are inserted, becoming increasingly
      solved over time. Non-filtering stores are pure accumulators.

    :param consistent: Whether to enable consistency checking.
    :type consistent: bool
    """

    num_vars: int
    num_eqs: int
    consistent: bool
    eager: bool
    filtering: bool

    def __init__(self, *, consistent: bool = False):
        self.consistent = consistent
        self.eager = True
        self.filtering = False

    @property
    def is_determined(self) -> bool:
        """Whether the stored equations fully determine all variables.

        Only meaningful for filtering stores. For eager stores, this
        is live after each insertion. For deferred stores, only valid
        after :meth:`process_pending`.

        :rtype: bool
        """
        return False

    def queue_equation(
        self,
        equation: Any,
        identifier: Any = None,
        translate_ANF: bool = True,
    ):
        """Accept an equation without triggering batch processing.

        For eager stores, equivalent to :meth:`insert_equation` (the
        "queue" is always immediately drained). For deferred stores,
        stages the equation for later :meth:`process_pending`.

        :param equation: The equation to insert.
        :param identifier: Optional metadata (clock time, cube, etc.).
        :param translate_ANF: Whether to convert to ANF internally.
        """
        return self.insert_equation(equation, identifier, translate_ANF)

    def process_pending(self, *, verbose: bool = False, batch_size: int | None = None) -> None:
        """Finalize any deferred reduction. No-op for eager/passive stores.

        :param verbose: Whether to print progress.
        :type verbose: bool
        :param batch_size: Process in batches of this size,
            returning control between batches (for progress reporting).
            ``None`` processes everything at once.
        :type batch_size: int | None
        """
        pass

    def insert_equation(
        self,
        equation: Any,
        identifier: Any = None,
        translate_ANF: bool = True,
    ):
        """Accept an equation. Subclasses define what types they accept.

        :param equation: The equation to insert.
        :param identifier: Optional metadata (clock time, cube, etc.).
        :param translate_ANF: Whether to convert to ANF internally.
        """
        raise NotImplementedError
