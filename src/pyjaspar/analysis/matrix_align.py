"""Search a set of profiles with a query motif.

``search_profiles`` scores a query motif against candidate profiles and ranks
them. The scoring is the one of JASPAR's "Matrix Align" web tool, which runs the
``matrix_aligner`` program (Sandelin et al., Funct Integr Genomics 3:125-134,
2003; source in the ``jaspar_tools`` repository of the JASPAR team): a
semi-global variant of the Needleman-Wunsch algorithm that permits at most one
internal gap. This is an independent numpy implementation of that method; the
default gap penalties are the program's defaults.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

try:
    import numpy as np
except ImportError:
    raise ImportError(
        "The pyjaspar.analysis module requires numpy. Install with: pip install pyjaspar[analysis]"
    ) from None

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif

_BASES = "ACGT"
_MATRIX_ID = re.compile(r"[A-Z]+(\d+)\.(\d+)")
# Scores whose relative difference is below this count as equal: a couple of
# units in the last place of single precision, which the reference program uses
# and which cannot tell such scores apart.
_TIE_TOLERANCE = 2e-7


@dataclass(frozen=True)
class AlignScore:
    """Result of aligning two profiles.

    Attributes:
        score: Alignment score. Each aligned column contributes between 0 and 2,
            and a gap subtracts its penalty.
        is_reverse_complement: True if the reverse complement of the second
            profile aligned better. When both orientations score the same, this
            is True.
        gaps: Number of gap columns in the best alignment (0 if it has no gap);
            the gap is one run of that many columns.
        offset: Position in the first profile minus position in the second at
            the first pair of aligned columns (positions count from 1, in the
            orientation of the second profile that was aligned).
        alignment_length: Number of columns of the alignment, gap columns
            included and the free overhangs excluded.
    """

    score: float
    is_reverse_complement: bool
    gaps: int
    offset: int
    alignment_length: int


@dataclass(frozen=True)
class ProfileHit:
    """A candidate profile scored against a query.

    Attributes:
        matrix_id: JASPAR matrix ID of the candidate.
        name: Name of the TF the candidate belongs to.
        score: Alignment score; this is the "Score" column of the JASPAR web
            tool. It is a sum over columns, so it tends to be higher for longer
            profiles.
        percent_score: The "Percent Score" column of the JASPAR web tool,
            ``100 * score / (2 * m)``, where ``m`` is the narrowest profile seen
            so far. ``m`` starts at the query's width and is lowered by each
            candidate in matrix-ID order (the order of the web tool's table),
            and never raised again. The value therefore depends on the set of
            candidates, and it can exceed 100.
        is_reverse_complement: True if the candidate's reverse complement
            aligned better (also True when both orientations score the same).
        gaps: Number of gap columns in the best alignment (0 if it has no gap).
        width: Number of columns of the candidate.
        offset: Position in the query minus position in the candidate at the
            first pair of aligned columns (positions count from 1, in the
            orientation of the candidate that was aligned).
        alignment_length: Number of columns of the alignment, gap columns
            included and the free overhangs excluded.
    """

    matrix_id: str
    name: str
    score: float
    percent_score: float
    is_reverse_complement: bool
    gaps: int
    width: int
    offset: int
    alignment_length: int


def _frequencies(motif: Motif) -> np.ndarray:
    """Return the motif's columns as frequencies, shape (width, 4)."""
    counts = np.array([motif.counts[b] for b in _BASES], dtype=float).T
    totals = counts.sum(axis=1, keepdims=True)
    if (totals <= 0).any():
        raise ValueError(f"Motif {motif.matrix_id or motif.name} has a column with no counts")
    return counts / totals


def _column_similarity(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Similarity of every column of ``a`` with every column of ``b``, in [0, 2].

    Two identical columns score 2; the squared Euclidean distance between two
    frequency vectors is at most 2.
    """
    return 2.0 - ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)


def _best_alignment(
    a: np.ndarray, b: np.ndarray, open_penalty: float, ext_penalty: float
) -> tuple[float, int, int, int]:
    """Best score, gap columns, offset and alignment length for one orientation."""
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
    gaps, offset, length = 0, i - j, min(i, j) + 1

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
                i, j = np.unravel_index(int(total.argmax()), total.shape)
                # The second diagonal starts after the gap. Columns with similarity
                # 0 at its end are not part of the alignment.
                i2, j2 = (i + 1, j + k + 1) if in_second else (i + k + 1, j + 1)
                steps = np.arange(min(n - i2, m - j2))
                positive = np.flatnonzero(sim[i2 + steps, j2 + steps] > 0)
                tail = int(positive[-1]) + 1 if positive.size else 1
                gaps, offset = k, int(i - j)
                length = min(int(i), int(j)) + 1 + k + tail
    return best, gaps, offset, length


def align_score(
    a: Motif,
    b: Motif,
    open_penalty: float = 3.0,
    ext_penalty: float = 0.01,
) -> AlignScore:
    """Align two profiles and return the best score.

    Columns are compared as frequencies. The alignment may leave columns
    hanging off either end for free and may contain one internal gap, which
    costs ``open_penalty`` for its first column and ``ext_penalty`` for each
    further column. Both ``b`` and its reverse complement are tried; when they
    score the same (up to rounding), the reverse complement is reported.

    Args:
        a: First profile.
        b: Second profile.
        open_penalty: Penalty for the first column of a gap.
        ext_penalty: Penalty for each further column of the same gap.

    Returns:
        The best score, the orientation of ``b`` that gave it, the number of
        gap columns, the offset of the first aligned pair, and the number of
        columns of the alignment.

    Raises:
        ValueError: If a profile has a column with no counts.
    """
    fa, fb = _frequencies(a), _frequencies(b)
    forward = _best_alignment(fa, fb, open_penalty, ext_penalty)
    reverse = _best_alignment(fa, fb[::-1, ::-1], open_penalty, ext_penalty)
    if forward[0] - reverse[0] > _TIE_TOLERANCE * max(forward[0], reverse[0]):
        return AlignScore(forward[0], False, *forward[1:])
    return AlignScore(reverse[0], True, *reverse[1:])


def _matrix_id_order(candidates: list[Motif]) -> list[int]:
    """Indices of the candidates in matrix-ID order, or as given if an ID is not standard."""
    keys = []
    for candidate in candidates:
        match = _MATRIX_ID.fullmatch(candidate.matrix_id or "")
        if match is None:
            return list(range(len(candidates)))
        keys.append((int(match.group(1)), int(match.group(2))))
    return sorted(range(len(candidates)), key=keys.__getitem__)


def search_profiles(
    query: Motif,
    candidates: Iterable[Motif],
    open_penalty: float = 3.0,
    ext_penalty: float = 0.01,
    sort_by: str = "score",
    top: int | None = None,
) -> list[ProfileHit]:
    """Score a query motif against candidate profiles and rank them.

    Build the candidates with ``JasparDB.fetch_motifs``, which offers the same
    filters as the web tool (collection, taxonomic group, latest or all
    versions).

    Args:
        query: The motif to search with.
        candidates: Profiles to compare the query with.
        open_penalty: Penalty for the first column of a gap.
        ext_penalty: Penalty for each further column of the same gap.
        sort_by: ``"score"`` or ``"percent_score"``; best first.
        top: Return only this many hits; all of them if None.

    Returns:
        The scored candidates, best first.

    Raises:
        ValueError: If ``sort_by`` is unknown or a profile has a column with
            no counts.
    """
    if sort_by not in ("score", "percent_score"):
        raise ValueError("sort_by must be 'score' or 'percent_score'")

    candidates = list(candidates)
    results = [align_score(query, c, open_penalty, ext_penalty) for c in candidates]

    # The web tool divides by the narrowest profile seen so far, going through
    # its table in matrix-ID order; see ProfileHit.percent_score.
    narrowest = query.length
    percent = {}
    for i in _matrix_id_order(candidates):
        narrowest = min(narrowest, candidates[i].length)
        percent[i] = 100 * results[i].score / (2 * narrowest)

    hits = [
        ProfileHit(
            matrix_id=c.matrix_id,
            name=c.name,
            score=r.score,
            percent_score=percent[i],
            is_reverse_complement=r.is_reverse_complement,
            gaps=r.gaps,
            width=c.length,
            offset=r.offset,
            alignment_length=r.alignment_length,
        )
        for i, (c, r) in enumerate(zip(candidates, results, strict=True))
    ]
    hits.sort(key=lambda h: (-getattr(h, sort_by), h.matrix_id))
    return hits if top is None else hits[:top]
