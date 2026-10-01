"""Online-phase insertion adapter for algebraic attacks.

The online phase streams coefficient vectors into a store one at a time.
Different stores require different input formats (ndarray vs BooleanFunction),
different submission strategies (eager insert vs deferred queue), and different
control flow (early stopping vs drain-then-process).

Rather than checking store properties on every iteration of the hot loop,
:func:`make_online_inserter` inspects the store once and returns specialized
closures that bake in the correct behavior. This is a metaprogramming pattern:
the factory "generates" the right function at setup time so the inner loop
runs without isinstance/hasattr/property checks.

The two dispatch axes are:

1. **Equation format** (``isinstance(store, IndexedEqStore)``):
   Indexed stores accept ndarray coefficient vectors directly.
   Non-indexed stores (GroebnerEqStore) need conversion to BooleanFunction.

2. **Early stopping** (``store.eager and store.filtering``):
   Eager+filtering stores (LUEqStore) know when the system is fully
   determined after each insert, enabling early loop termination.
   All other stores either defer processing (so determinacy isn't known
   yet) or don't filter (so they can never be "determined").

These axes are orthogonal, yielding four closure variants. Each variant
contains only the logic relevant to that combination — no dead branches.
"""
from PyPR.Reporting import get_logger

from PyPR.Cryptanalysis.Components.Adapters.equation_repr import (
    coef_vector_to_boolean_function,
)
from PyPR.Cryptanalysis.Components.EquationStores.IndexedEqStore import IndexedEqStore

log = get_logger(__name__)


def make_online_inserter(store, idx_to_comb, *, total_eqs, num_vars,
                         stall_limit=100):
    """Inspect store properties once, return specialized hot-loop closures.

    The returned ``insert_fn`` and ``finalize_fn`` capture the store and
    all loop-invariant state (index maps, the progress meter) in their
    closure. The caller uses them as::

        insert_eq, finalize = make_online_inserter(store, idx_to_comb, ...)
        for eq_idx in range(n):
            coef_vector = ...
            if insert_eq(coef_vector, eq_idx):
                break
        finalize()

    Progress is reported through a meter opened here, so it belongs to the
    step the caller is in when it builds the inserter and closes with it.

    :param store: The online equation store.
    :param idx_to_comb: Column-index-to-monomial-tuple mapping.
    :param total_eqs: Total equations expected (for progress display).
    :type total_eqs: int
    :param num_vars: Number of variables (for progress display).
    :type num_vars: int
    :param stall_limit: Stop early when the filtering metric (rank, basis
        size, etc.) does not change for this many consecutive steps.
        ``None`` disables stall detection.
    :type stall_limit: int | None
    :return: ``(insert_fn, finalize_fn)`` — insert_fn(coef_vector, eq_idx) -> bool
        (True = early stop), finalize_fn() -> None.
    :rtype: tuple[Callable, Callable]
    """
    # --- Dispatch axis 1: equation format ---
    is_indexed = isinstance(store, IndexedEqStore)

    # --- Dispatch axis 2: early stopping ---
    can_early_stop = store.eager and store.filtering

    # --- Stall metric ---
    # Only a filtering store makes progress to measure; each one reports it as
    # num_determined, by whatever measure it reduces with. Bound unconditionally
    # so the closures below always have something to call -- a non-filtering
    # store never reaches the stall check, and reads a constant if it did.
    _get_solved = (lambda: store.num_determined) if store.filtering else (lambda: 0)

    # A store that queues equations for later only enqueues here; the work
    # happens in finalize_fn.
    substituted = log.progress(
        "Equations substituted" if store.eager else "Equations enqueued",
        total=total_eqs, unit="eq",
    )

    def _stop(reason, *args):
        substituted.close()
        log.info(reason, *args)
        return True

    # ---------------------------------------------------------------
    # Build the insert closure.
    #
    # Four variants from the two boolean axes; all type checks and
    # property lookups are resolved here at factory time.
    #
    # Stall detection (for eager+filtering): if the solved-variable
    # count doesn't increase for stall_limit consecutive insertions,
    # stop early — remaining equations are likely all redundant.
    # ---------------------------------------------------------------

    if is_indexed and can_early_stop:
        # Typical case: LUEqStore.  Insert ndarray, check rank.
        _stall_count = 0
        _last_solved = _get_solved()

        def _insert_indexed_stopping(coef_vector, eq_idx):
            nonlocal _stall_count, _last_solved
            store.queue_equation(coef_vector, identifier=eq_idx)

            current = _get_solved()
            if current != _last_solved:
                _stall_count = 0
                _last_solved = current
            else:
                _stall_count += 1

            substituted.update()
            if substituted.shown:
                substituted.set_status(f"Rank: {store.rank}/{num_vars}")
            if store.is_determined:
                return _stop("System determined after %d/%d equations", eq_idx + 1, total_eqs)
            if stall_limit is not None and _stall_count >= stall_limit:
                return _stop("Stalled for %d insertions -- stopping (solved %d/%d)",
                             stall_limit, _last_solved, num_vars)
            return False

        insert_fn = _insert_indexed_stopping

    elif is_indexed:
        # Indexed but passive (EqStore / SymbolicEqStore as online store).
        # Insert ndarray, never stop early.
        def _insert_indexed_passive(coef_vector, eq_idx):
            store.queue_equation(coef_vector, identifier=eq_idx)
            substituted.update()
            return False

        insert_fn = _insert_indexed_passive

    elif can_early_stop:
        # Non-indexed, eager+filtering.  Convert to BooleanFunction.
        _stall_count = 0
        _last_solved = _get_solved()

        def _insert_converting_stopping(coef_vector, eq_idx):
            nonlocal _stall_count, _last_solved
            bf = coef_vector_to_boolean_function(coef_vector, idx_to_comb)
            store.queue_equation(bf)

            current = _get_solved()
            if current != _last_solved:
                _stall_count = 0
                _last_solved = current
            else:
                _stall_count += 1

            substituted.update()
            if substituted.shown:
                substituted.set_status(f"Solved: {_last_solved}/{num_vars}")
            if store.is_determined:
                return _stop("System determined after %d/%d equations", eq_idx + 1, total_eqs)
            if stall_limit is not None and _stall_count >= stall_limit:
                return _stop("Stalled for %d insertions -- stopping (solved %d/%d)",
                             stall_limit, _last_solved, num_vars)
            return False

        insert_fn = _insert_converting_stopping

    else:
        # Typical case: GroebnerEqStore.  Convert to BooleanFunction,
        # queue for batch processing, never stop early.
        def _insert_converting_passive(coef_vector, eq_idx):
            bf = coef_vector_to_boolean_function(coef_vector, idx_to_comb)
            store.queue_equation(bf)
            substituted.update()
            return False

        insert_fn = _insert_converting_passive

    # ---------------------------------------------------------------
    # Build the finalize closure.
    #
    # Eager stores have nothing to finalize (reduction already happened
    # per-insert).  Deferred stores need process_pending() to trigger
    # their batch computation (e.g., Gröbner basis).
    #
    # Stall detection (for deferred+filtering): if the solved-variable
    # count doesn't increase for stall_limit consecutive batches,
    # stop — further processing is unlikely to determine new variables.
    # ---------------------------------------------------------------

    if store.eager:
        def finalize_fn():
            substituted.close()
    else:
        _batch_size = 10
        has_queue = hasattr(store, 'queue')

        if store.filtering:
            def finalize_fn():
                substituted.close()
                processed = 0
                stall_count = 0
                last_solved = _get_solved()

                # the queue grows as processing adds syzygies, so there is no
                # fixed total to count toward
                processing = log.progress("Pending equations processed", unit="eq")
                while has_queue and store.queue:
                    processed += store.process_pending(batch_size=_batch_size)

                    current = _get_solved()
                    if current != last_solved:
                        stall_count = 0
                        last_solved = current
                    else:
                        stall_count += 1

                    processing.update_to(processed)
                    if processing.shown:
                        processing.set_status(f"Queue: {len(store.queue)} -- Solved: {last_solved}")

                    # No early-exit on store.is_determined here: for a
                    # deferred store (GroebnerEqStore) that only means every
                    # variable seen so far has a value -- a still-queued
                    # equation could yet reduce to a contradiction. The
                    # while condition above already stops once the queue is
                    # genuinely empty (fully drained, nothing left to check).

                    if stall_limit is not None and stall_count >= stall_limit:
                        processing.close()
                        log.info("Stalled for %d batches -- stopping (solved %d)", stall_limit, last_solved)
                        break
                processing.close()
        else:
            # Deferred but non-filtering: no early stopping possible.
            def finalize_fn():
                substituted.close()
                store.process_pending()

    return insert_fn, finalize_fn
