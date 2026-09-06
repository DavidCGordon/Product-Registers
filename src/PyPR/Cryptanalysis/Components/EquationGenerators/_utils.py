"""Shared utilities for equation generators."""
from PyPR.BooleanLogic import BooleanFunction


def normalize_output_fn(output_fn):
    """Normalize output_fn to a (list, return_list) pair.

    All three generators accept either a single BooleanFunction or
    a list of them, then yield items matching that shape. This
    function handles the check and normalization.

    :param output_fn: A single BooleanFunction or a list of them.
    :type output_fn: BooleanFunction | list[BooleanFunction]
    :return: ``(output_fn_list, return_list)`` — the list form
        and a flag indicating whether the caller passed a list.
    :rtype: tuple[list[BooleanFunction], bool]
    """
    if type(output_fn) == list:
        return output_fn, True
    elif isinstance(output_fn, BooleanFunction):
        return [output_fn], False
    else:
        raise TypeError(
            f"output_fn must be a BooleanFunction or list of Boolean functions. "
            f"Got {type(output_fn)} instead."
        )
