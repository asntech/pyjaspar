"""Motif similarity and distance measures.

Provides column-wise comparison of position frequency matrices (PFMs).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "The pyjaspar.analysis module requires numpy. Install with: pip install pyjaspar[analysis]"
    ) from None

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif


def _get_pfm_columns(motif: Motif) -> list[list[float]]:
    """Extract normalized frequency columns from a motif."""
    counts = motif.counts
    length = motif.length
    columns: list[list[float]] = []
    for i in range(length):
        col = [float(counts[nt][i]) for nt in "ACGT"]
        total = sum(col)
        if total > 0:
            col = [v / total for v in col]
        columns.append(col)
    return columns


def pearson_correlation(
    motif1: Motif,
    motif2: Motif,
    offset: int = 0,
) -> float:
    """Compute average column-wise Pearson correlation between two motifs.

    Aligns motif2 at the given offset relative to motif1. Only the
    overlapping region is considered.

    Args:
        motif1: First motif.
        motif2: Second motif.
        offset: Position offset of motif2 relative to motif1.
            Positive = motif2 shifted right, negative = shifted left.

    Returns:
        Average Pearson correlation coefficient in [-1, 1].
        Returns 0.0 if there is no overlap.
    """
    cols1 = _get_pfm_columns(motif1)
    cols2 = _get_pfm_columns(motif2)

    # Determine overlapping range
    start1 = max(0, offset)
    start2 = max(0, -offset)
    overlap = min(len(cols1) - start1, len(cols2) - start2)

    if overlap <= 0:
        return 0.0

    correlations: list[float] = []
    for i in range(overlap):
        c1 = np.array(cols1[start1 + i])
        c2 = np.array(cols2[start2 + i])
        std1 = np.std(c1)
        std2 = np.std(c2)
        if std1 == 0 or std2 == 0:
            correlations.append(0.0)
        else:
            r = float(np.corrcoef(c1, c2)[0, 1])
            correlations.append(r)

    return float(np.mean(correlations))


def best_correlation(
    motif1: Motif,
    motif2: Motif,
    min_overlap: int = 4,
    both_strands: bool = True,
) -> tuple[float, int, bool]:
    """Find the best Pearson correlation across all valid offsets and orientations.

    Args:
        motif1: First motif (reference).
        motif2: Second motif.
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
        score = pearson_correlation(motif1, motif2, offset)
        if score > best_score:
            best_score = score
            best_offset = offset
            best_rc = False

    if both_strands:
        rc_motif2 = motif2.reverse_complement()
        for offset in range(-(len2 - min_overlap), len1 - min_overlap + 1):
            score = pearson_correlation(motif1, rc_motif2, offset)
            if score > best_score:
                best_score = score
                best_offset = offset
                best_rc = True

    return best_score, best_offset, best_rc


def euclidean_distance(
    motif1: Motif,
    motif2: Motif,
    offset: int = 0,
) -> float:
    """Compute average column-wise Euclidean distance between two motifs.

    Args:
        motif1: First motif.
        motif2: Second motif.
        offset: Position offset of motif2 relative to motif1.

    Returns:
        Average Euclidean distance per column. Lower = more similar.
        Returns float('inf') if there is no overlap.
    """
    cols1 = _get_pfm_columns(motif1)
    cols2 = _get_pfm_columns(motif2)

    start1 = max(0, offset)
    start2 = max(0, -offset)
    overlap = min(len(cols1) - start1, len(cols2) - start2)

    if overlap <= 0:
        return float("inf")

    total = 0.0
    for i in range(overlap):
        c1 = np.array(cols1[start1 + i])
        c2 = np.array(cols2[start2 + i])
        total += float(np.sqrt(np.sum((c1 - c2) ** 2)))

    return total / overlap


def kl_divergence(
    motif1: Motif,
    motif2: Motif,
    offset: int = 0,
    pseudocount: float = 0.001,
) -> float:
    """Compute average column-wise KL divergence from motif1 to motif2.

    Args:
        motif1: Reference motif (P distribution).
        motif2: Comparison motif (Q distribution).
        offset: Position offset of motif2 relative to motif1.
        pseudocount: Small value added to avoid log(0).

    Returns:
        Average KL divergence per column (non-negative).
        Returns float('inf') if there is no overlap.
    """
    cols1 = _get_pfm_columns(motif1)
    cols2 = _get_pfm_columns(motif2)

    start1 = max(0, offset)
    start2 = max(0, -offset)
    overlap = min(len(cols1) - start1, len(cols2) - start2)

    if overlap <= 0:
        return float("inf")

    total = 0.0
    for i in range(overlap):
        p = [v + pseudocount for v in cols1[start1 + i]]
        q = [v + pseudocount for v in cols2[start2 + i]]
        p_sum = sum(p)
        q_sum = sum(q)
        p = [v / p_sum for v in p]
        q = [v / q_sum for v in q]
        kl = sum(pi * math.log(pi / qi) for pi, qi in zip(p, q, strict=True))
        total += kl

    return total / overlap
