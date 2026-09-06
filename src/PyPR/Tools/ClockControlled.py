"""Shrinking, self-shrinking, and alternating-step generator constructions
on arbitrary binary sequences.

Each function takes raw bit lists and returns the combined output sequence.
These operate on finite prefixes — supply enough input to cover at least one
full period of each source if periodic LC analysis is intended.
"""

from __future__ import annotations


def shrinking_generator(
    control: list[int],
    data: list[int],
) -> list[int]:
    """Output data[t] at every t where control[t] == 1.

    :param control: Binary control sequence (the "shrinking sequence").
    :type control: list[int]
    :param data: Binary data sequence (the "shrunken sequence").
    :type data: list[int]
    :return: The subsequence of data selected by the control.
    :rtype: list[int]
    """
    return [d for c, d in zip(control, data) if c == 1]


def self_shrinking_generator(
    seq: list[int],
) -> list[int]:
    """Pair consecutive bits; output the second when the first is 1.

    :param seq: Binary source sequence, consumed in pairs.
    :type seq: list[int]
    :return: The self-shrunk output.
    :rtype: list[int]
    """
    return [seq[2 * t + 1] for t in range(len(seq) // 2) if seq[2 * t] == 1]


def alternating_step_generator(
    control: list[int],
    data_b: list[int],
    data_c: list[int],
) -> list[int]:
    r"""ASG: at each step, advance B when control=1 or C when control=0.

    Output is $y[t] = B[N_1(t)] \oplus C[N_0(t)]$ where $N_1(t)$ and $N_0(t)$
    count the ones and zeros of the control up to (not including) time $t$.

    :param control: Binary control sequence.
    :type control: list[int]
    :param data_b: Binary data sequence B (advances on control=1).
    :type data_b: list[int]
    :param data_c: Binary data sequence C (advances on control=0).
    :type data_c: list[int]
    :return: The ASG output sequence.
    :rtype: list[int]
    """
    out = []
    n1, n0 = 0, 0
    for c in control:
        out.append(data_b[n1] ^ data_c[n0])
        if c == 1:
            n1 += 1
        else:
            n0 += 1
    return out
