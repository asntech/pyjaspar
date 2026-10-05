"""Matrix Align numerical kernel.

The scoring of JASPAR's "Matrix Align" web tool, which runs the
``matrix_aligner`` program (Sandelin et al., Funct Integr Genomics 3:125-134,
2003; source in the ``jaspar_tools`` repository of the JASPAR team): a
semi-global variant of the Needleman-Wunsch algorithm that permits at most one
internal gap. This is an independent numpy implementation of that method.

This module works on frequency arrays and plain tuples only; the public result
types are built in ``methods.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif

_BASES = "ACGT"
# Scores whose relative difference is below this count as equal: a couple of
# units in the last place of single precision, which the reference program uses
# and which cannot tell such scores apart.
_TIE_TOLERANCE = 2e-7

# Alignment path as zero-based column boundaries: (query vertices, target vertices).
Coordinates = tuple[tuple[int, ...], tuple[int, ...]]


def _frequencies(motif: Motif) -> np.ndarray:
    """Return the motif's columns as frequencies, shape (width, 4)."""
    counts = np.array([motif.counts[b] for b in _BASES], dtype=float).T
    totals = counts.sum(axis=1, keepdims=True)
    if (totals <= 0).any():
        label = getattr(motif, "matrix_id", None) or getattr(motif, "name", None) or "motif"
        raise ValueError(f"Motif {label} has a column with no counts")
    return counts / totals


def _column_similarity(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Similarity of every column of ``a`` with every column of ``b``, in [0, 2].

    Two identical columns score 2; the squared Euclidean distance between two
    frequency vectors is at most 2.
    """
    return 2.0 - ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)


def _best_alignment(
    a: np.ndarray, b: np.ndarray, open_penalty: float, ext_penalty: float
) -> tuple[float, Coordinates]:
    """Best score and the path of the best alignment for one orientation."""
    sim = _column_similarity(a, b)
    n, m = sim.shape

    # Cumulative similarity along each diagonal, from its start (prefix) and
    # to its end (suffix). Overhanging columns are free, and every similarity
    # is >= 0, so an ungapped alignment always runs to the end of its diagonal;
    # this shortcut breaks if a similarity that can be negative is used.
    prefix = sim.copy()
    for i in range(1, n):
        prefix[i, 1:] += prefix[i - 1, :-1]
    suffix = sim.copy()
    for i in range(n - 2, -1, -1):
        suffix[i, :-1] += suffix[i + 1, 1:]

    # Ungapped: the alignment ends at the best cell of ``prefix`` and starts at
    # the beginning of its diagonal. Columns with similarity 0 at the end do not
    # lower the score, and the program counts them in the alignment, so the last
    # best cell is taken.
    flat = prefix.ravel()
    cell = int(np.flatnonzero(flat == flat.max())[-1])
    i, j = divmod(cell, m)
    best = float(flat[cell])
    start = min(i, j)
    path: Coordinates = ((i - start, i + 1), (j - start, j + 1))

    # One gap run of k columns: the alignment follows one diagonal up to a cell
    # (i, j), skips k columns of one profile, then follows the next diagonal.
    for k in range(1, max(n, m)):
        cost = open_penalty + (k - 1) * ext_penalty
        for in_second in (False, True):
            if in_second:
                if not (m - k - 1 >= 1 and n >= 2):
                    continue
                total = prefix[: n - 1, : m - k - 1] + suffix[1:, k + 1 :]
            else:
                if not (n - k - 1 >= 1 and m >= 2):
                    continue
                total = prefix[: n - k - 1, : m - 1] + suffix[k + 1 :, 1:]
            gapped = float(total.max()) - cost
            if gapped > best:
                best = gapped
                i, j = (int(x) for x in np.unravel_index(int(total.argmax()), total.shape))
                # The second diagonal starts after the gap: k columns of the second
                # profile are skipped when ``in_second``, of the first otherwise.
                # Columns with similarity 0 at its end are not part of the alignment.
                i2, j2 = (i + 1, j + k + 1) if in_second else (i + k + 1, j + 1)
                steps = np.arange(min(n - i2, m - j2))
                positive = np.flatnonzero(sim[i2 + steps, j2 + steps] > 0)
                tail = int(positive[-1]) + 1 if positive.size else 1
                start = min(i, j)
                path = (
                    (i - start, i + 1, i2, i2 + tail),
                    (j - start, j + 1, j2, j2 + tail),
                )
    return best, path


def align_frequencies(
    a: np.ndarray,
    b: np.ndarray,
    open_penalty: float,
    ext_penalty: float,
    both_strands: bool,
) -> tuple[float, bool, Coordinates]:
    """Best score, orientation of ``b`` and path; the reverse complement wins ties."""
    forward = _best_alignment(a, b, open_penalty, ext_penalty)
    if not both_strands:
        return forward[0], False, forward[1]
    reverse = _best_alignment(a, b[::-1, ::-1], open_penalty, ext_penalty)
    if forward[0] - reverse[0] > _TIE_TOLERANCE * max(forward[0], reverse[0]):
        return forward[0], False, forward[1]
    return reverse[0], True, reverse[1]
