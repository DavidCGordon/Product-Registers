"""Base class for equation stores that map monomials to column indices.

Handles dynamic/static mode, comb_to_idx/idx_to_comb maps, linking,
and monomial expansion — all currently duplicated across EqStore,
LUEqStore, and SymbolicEqStore.

Subclasses implement ``_expand_storage()`` to resize their specific
backing arrays when the monomial set grows.
"""
from __future__ import annotations

from typing import Any

from PyPR.Cryptanalysis.Components.EquationStores.BaseEqStore import BaseEqStore


class IndexedEqStore(BaseEqStore):
    """Base for stores that map monomials to column indices.

    **Dynamic mode** (``comb_to_idx=None``): monomials are discovered
    on insertion. Backing storage starts overallocated and doubles
    as needed. Cannot accept raw ndarray input.

    **Static mode** (``comb_to_idx`` provided): fixed monomial set.
    Can accept ndarray directly (indices are known). Cannot add
    unknown monomials.

    :param comb_to_idx: Monomial-to-index mapping.
        ``None`` for dynamic mode, a dict for static mode.
    :type comb_to_idx: dict[tuple[int, ...], int] | None
    :param consistent: Whether to enable consistency checking.
        Passed to :class:`BaseEqStore`; enforcement is store-specific.
    :type consistent: bool
    """

    linked_stores: set[IndexedEqStore]
    comb_to_idx: dict[tuple[int, ...], int]
    idx_to_comb: dict[int, tuple[int, ...]]
    dynamic: bool

    def __init__(
        self,
        comb_to_idx: dict[tuple[int, ...], int] | None = None,
        consistent: bool = False,
    ):
        super().__init__(consistent=consistent)
        self.linked_stores: set[IndexedEqStore] = set([self])
        self.equation_ids: dict[int, Any] = {}
        self.num_eqs = 0

        if comb_to_idx is None:
            self.dynamic = True
            self.comb_to_idx = {}
            self.idx_to_comb = {}
            self.num_vars = 0
        else:
            self.dynamic = False
            self.comb_to_idx = {k: v for k, v in comb_to_idx.items()}
            self.idx_to_comb = {v: k for k, v in comb_to_idx.items()}
            self.num_vars = len(self.comb_to_idx)

    def link(self, other_store: IndexedEqStore) -> None:
        """Create a link to another store for monomial propagation.

        When this store discovers new monomials, they are propagated
        to all linked stores. This is directional: linking sends all
        currently known monomials to ``other_store`` immediately.

        Static stores should be source nodes in the linkage graph
        (adding an unknown monomial to a static store raises ValueError).

        :param other_store: The store to receive this store's monomial updates.
        :type other_store: IndexedEqStore
        """
        self.linked_stores.add(other_store)
        for monomial in self.comb_to_idx:
            other_store._update_known_monomials(monomial)

    def _update_known_monomials(self, comb: tuple[int, ...]) -> None:
        """Register a new monomial, expanding backing storage if needed.

        If the monomial is already known, this is a no-op. Otherwise:
        1. Calls ``_expand_storage()`` if capacity is exceeded
        2. Inserts into index maps
        3. Propagates to all linked stores (only when the monomial is
           genuinely new, preventing infinite recursion)

        :param comb: The monomial as a sorted tuple of variable indices.
        :type comb: tuple[int, ...]
        :raises ValueError: If called on a static store with an unknown monomial.
        """
        if comb in self.comb_to_idx:
            return

        if not self.dynamic:
            # More descriptive error messages for static and consistent stores
            # with no constant column (so no equation can contain a constant term):
            if self.consistent and comb == tuple():
                raise ValueError(
                    "Inserted equation has a constant term, but this store was marked "
                    "consistent without a constant column. All equations must have zero "
                    "constant term for consistency to hold."
                )
            # Normal error for any unknown monomial in a static store:
            raise ValueError(
                f"Monomial {comb} not in index maps for static {type(self)}"
            )

        # expands only when needed - not every time
        self._expand_storage()

        self.idx_to_comb[self.num_vars] = comb
        self.comb_to_idx[comb] = self.num_vars
        self.num_vars += 1

        for store in self.linked_stores:
            if store is not self:
                store._update_known_monomials(comb)

    def _expand_storage(self) -> None:
        """Resize backing arrays/matrices when capacity is exceeded.

        Called by ``_update_known_monomials`` before inserting a new
        monomial. Subclasses override this to resize their specific
        backing data structures.
        """
        raise NotImplementedError
