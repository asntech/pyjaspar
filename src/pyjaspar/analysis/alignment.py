"""Pairwise motif comparison and its text display.

``align_motifs`` compares two motifs with one of the methods in ``methods.py``
(``"matrix_align"`` by default, or ``"pearson"``) and returns the score and the
aligned columns. ``format_alignment`` turns such a result into a side-by-side
text alignment of the two consensus sequences.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from Bio.Align import Alignment

from .methods import AlignmentPath, AlignmentResult, _compare_pair, _resolve_options

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif


def align_motifs(
    motif1: Motif,
    motif2: Motif,
    *,
    method: str = "matrix_align",
    both_strands: bool = True,
    min_overlap: int | None = None,
    open_penalty: float | None = None,
    ext_penalty: float | None = None,
) -> AlignmentResult:
    """Compare two motifs and return the score and the aligned columns.

    Args:
        motif1: First motif (the query).
        motif2: Second motif; it is also tried reverse-complemented.
        method: ``"matrix_align"`` (default) or ``"pearson"``; see ``methods.py``.
        both_strands: If True, also try the reverse complement of ``motif2``.
        min_overlap: ``"pearson"`` only: minimum overlapping columns (default 4).
        open_penalty: ``"matrix_align"`` only: cost of a gap's first column (default 3.0).
        ext_penalty: ``"matrix_align"`` only: cost of each further gap column (default 0.01).

    Returns:
        An ``AlignmentResult``. Use ``format_alignment`` to display it.

    Raises:
        ValueError: If the method or an option is invalid, an option of the other
            method is given, a motif has invalid counts, or (``"pearson"``) a motif
            is shorter than ``min_overlap``.
    """
    options = _resolve_options(
        method,
        both_strands=both_strands,
        min_overlap=min_overlap,
        open_penalty=open_penalty,
        ext_penalty=ext_penalty,
    )
    return _compare_pair(motif1, motif2, method, options)


def _display_coordinates(path: AlignmentPath) -> np.ndarray:
    """The path extended to both full motifs, for display only.

    End columns that hang off are shown: where both motifs have unscored columns
    at an end, they are shown side by side up to the shorter flank, and the rest of
    the longer flank is shown against gaps.
    """
    query, target = path.coordinates
    n, m = path.query_length, path.target_length
    paired_start = min(query[0], target[0])
    paired_end = min(n - query[-1], m - target[-1])
    vertices = [(0, 0), (query[0] - paired_start, target[0] - paired_start)]
    vertices += list(zip(query, target, strict=True))
    vertices += [(query[-1] + paired_end, target[-1] + paired_end), (n, m)]
    kept = [vertices[0]]
    for vertex in vertices[1:]:
        if vertex != kept[-1]:
            kept.append(vertex)
    return np.array(kept).T


def format_alignment(motif1: Motif, motif2: Motif, result: AlignmentResult) -> str:
    """Text alignment of the two consensus sequences for a comparison result.

    Shows both motifs over their full length, including the columns that hang
    off either end, with the second motif in the orientation that was aligned.
    The alignment is not recomputed: the aligned columns come from ``result``.
    The match line compares consensus letters only; it is a visual aid, not the
    profile score.

    Args:
        motif1: The first motif given to ``align_motifs``.
        motif2: The second motif given to ``align_motifs``.
        result: The result of ``align_motifs(motif1, motif2, ...)``.

    Returns:
        The alignment block as text.

    Raises:
        ValueError: If the motif widths do not match the result.
    """
    path = result.alignment
    if motif1.length != path.query_length or motif2.length != path.target_length:
        raise ValueError(
            f"The motifs have {motif1.length} and {motif2.length} columns, but the result "
            f"was computed for {path.query_length} and {path.target_length}"
        )
    motif2_used = motif2.reverse_complement() if path.is_reverse_complement else motif2
    sequences = [str(motif1.consensus), str(motif2_used.consensus)]
    return str(Alignment(sequences, _display_coordinates(path)))
