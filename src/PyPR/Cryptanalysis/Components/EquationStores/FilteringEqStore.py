"""Base class for stores that converge toward a solution as equations arrive.

Separates the stores that make *progress* from those that merely accumulate,
and gives that progress a single name every consumer can read regardless of
how the store represents it internally.

See `docs/architecture/Components_Architecture.md` §1 for how this axis crosses
`IndexedEqStore`, and for why the progress figure is a member of the contract
rather than a duck-typed check at each call site.
"""
from __future__ import annotations

from PyPR.Cryptanalysis.Components.EquationStores.BaseEqStore import BaseEqStore


class FilteringEqStore(BaseEqStore):
    """Base for stores that discard redundant information as they reduce.

    A non-filtering store is a pure accumulator: every equation inserted is
    still there afterwards, and the store is no closer to a solution than the
    equations it was handed. A filtering store reduces instead — it rejects
    what its existing contents already imply, so the portion of the system it
    has pinned down grows monotonically toward the whole.

    That pinned-down portion is what :attr:`num_determined` reports, and the
    two families measure it by different means: an LU store counts pivot
    columns (its rank), a Groebner store counts variables its reduction has
    driven to a constant. The quantities are not the same notion of "solved" —
    rank counts constrained dimensions, which need back-substitution before any
    variable has a known value — and they are not comparable between stores.
    What they share is the only thing consumers rely on: both increase only
    when the store has genuinely learned something, and both reach
    :attr:`~BaseEqStore.num_vars` exactly when the system is determined.

    Subclasses must implement :attr:`num_determined` and
    :attr:`~BaseEqStore.is_determined`, and must maintain
    :attr:`~BaseEqStore.num_vars` as the count of variables the store has ever
    seen, so that the completion invariant below holds.

    :cvar filtering: Always ``True`` for this family. Declared on the class
        rather than assigned per instance so that it cannot disagree with the
        type; ``store.filtering`` and ``isinstance(store, FilteringEqStore)``
        are the same question.
    :vartype filtering: bool
    """

    filtering = True

    @property
    def num_determined(self) -> int:
        """How much of the system this store has pinned down.

        Monotone non-decreasing across insertions into the same store, and
        equal to :attr:`~BaseEqStore.num_vars` exactly when
        :attr:`~BaseEqStore.is_determined` is ``True``. Meaningful as a
        progress measure for one store over time; **not** comparable across
        stores of different families, which count different things.

        :return: The number of variables the store has determined, by
            whatever measure the concrete store reduces with.
        :rtype: int
        :raises NotImplementedError: If not implemented for the store class.
        """
        raise NotImplementedError
