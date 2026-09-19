"""Motif alignment view.

Renders the best offset/orientation between two motifs as an actual
gapped, side-by-side alignment, rather than just a score.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from Bio.Align import Alignment

from ._alignment_search import find_best_offset
from .similarity import pearson_correlation

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif


@dataclass
class AlignmentResult:
    """Rendered alignment between two motifs.

    Attributes:
        motif1_id: First motif's JASPAR matrix ID.
        motif2_id: Second motif's JASPAR matrix ID.
        score: Pearson correlation at the best offset.
        offset: Position offset of motif2 relative to motif1.
        is_reverse_complement: Whether motif2 was reverse-complemented.
        alignment: The rendered Bio.Align.Alignment. ``str(alignment)``
            gives a target/query display with a match/mismatch line.
    """

    motif1_id: str
    motif2_id: str
    score: float
    offset: int
    is_reverse_complement: bool
    alignment: Alignment


def _offset_to_coordinates(len1: int, len2: int, offset: int) -> np.ndarray:
    """Convert an alignment offset into Bio.Align.Alignment coordinates.

    Positive offset means sequence 2 is shifted right relative to sequence 1.

    Args:
        len1: Length of the first sequence.
        len2: Length of the second sequence.
        offset: Position offset of sequence 2 relative to sequence 1.

    Returns:
        A (2, n) array of per-sequence coordinates suitable for
        ``Bio.Align.Alignment(sequences, coordinates)``.
    """
    shift = max(0, -offset)
    s1_start, s1_end = shift, shift + len1
    s2_start, s2_end = shift + offset, shift + offset + len2
    overlap_start = max(s1_start, s2_start)
    overlap_end = min(s1_end, s2_end)

    points = sorted({s1_start, s2_start, overlap_start, overlap_end, s1_end, s2_end})
    c1 = [max(0, min(len1, p - s1_start)) for p in points]
    c2 = [max(0, min(len2, p - s2_start)) for p in points]
    return np.array([c1, c2])


def align_motifs(
    motif1: Motif,
    motif2: Motif,
    min_overlap: int = 4,
    both_strands: bool = True,
) -> AlignmentResult:
    """Align two motifs and render the result.

    Finds the best offset and orientation between the two motifs, then
    builds the actual gapped alignment view for it.

    Args:
        motif1: First motif (reference).
        motif2: Second motif.
        min_overlap: Minimum overlapping columns required.
        both_strands: If True, also try the reverse complement of motif2.

    Returns:
        AlignmentResult with the best score/offset/orientation plus the
        rendered alignment.
    """
    score, offset, is_rc = find_best_offset(
        motif1, motif2, pearson_correlation, min_overlap, both_strands
    )
    motif2_used = motif2.reverse_complement() if is_rc else motif2

    seq1 = str(motif1.consensus)
    seq2 = str(motif2_used.consensus)
    coordinates = _offset_to_coordinates(len(seq1), len(seq2), offset)
    alignment = Alignment([seq1, seq2], coordinates)

    return AlignmentResult(
        motif1_id=motif1.matrix_id,
        motif2_id=motif2.matrix_id,
        score=score,
        offset=offset,
        is_reverse_complement=is_rc,
        alignment=alignment,
    )
