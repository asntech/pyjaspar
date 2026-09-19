"""Offset/orientation search for aligning two motifs' PFM columns.

Finds the relative shift and strand orientation between two motifs that
maximizes a given column-similarity metric. Used by
similarity.best_correlation() and alignment.align_motifs().
"""

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif


def find_best_offset(
    motif1: Motif,
    motif2: Motif,
    metric: Callable[[Motif, Motif, int], float],
    min_overlap: int = 4,
    both_strands: bool = True,
) -> tuple[float, int, bool]:
    """Find the offset and orientation that maximizes a column-similarity metric.

    Tries every valid offset, and (if both_strands) the reverse complement
    of motif2, keeping whichever scores highest.

    `metric` must be a higher-is-better score (true for
    similarity.pearson_correlation; NOT true for
    similarity.euclidean_distance/kl_divergence, which are lower-is-better
    -- passing one of those in would silently return the worst alignment
    instead of the best, since this function always keeps the maximum).

    Args:
        motif1: First motif (reference).
        motif2: Second motif.
        metric: Column-similarity function, called as metric(motif1, motif2, offset).
        min_overlap: Minimum number of overlapping columns required.
        both_strands: If True, also try the reverse complement of motif2.

    Returns:
        Tuple of (best_score, best_offset, is_reverse_complement).
    """
    best_score = -2.0
    best_offset = 0
    best_rc = False

    len1 = motif1.length
    len2 = motif2.length

    for offset in range(-(len2 - min_overlap), len1 - min_overlap + 1):
        score = metric(motif1, motif2, offset)
        if score > best_score:
            best_score = score
            best_offset = offset
            best_rc = False

    if both_strands:
        rc_motif2 = motif2.reverse_complement()
        for offset in range(-(len2 - min_overlap), len1 - min_overlap + 1):
            score = metric(motif1, rc_motif2, offset)
            if score > best_score:
                best_score = score
                best_offset = offset
                best_rc = True

    return best_score, best_offset, best_rc
