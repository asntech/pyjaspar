"""Search a set of profiles with a query motif.

Scores a query motif against candidate profiles with the alignment used by
JASPAR's "Matrix Align" web tool: a semi-global variant of the
Needleman-Wunsch algorithm that permits at most one internal gap
(Sandelin et al., Funct Integr Genomics 3:125-134, 2003, as documented for
TFBS::Matrix::Alignment; the default gap penalties are the documented ones).
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


@dataclass(frozen=True)
class AlignScore:
    """Result of aligning two profiles.

    Attributes:
        score: Alignment score. Each aligned column contributes between 0 and 2,
            and a gap subtracts its penalty.
        is_reverse_complement: True if the reverse complement of the second
            profile aligned better.
        gaps: Number of internal gap runs in the best alignment (0 or 1).
    """

    score: float
    is_reverse_complement: bool
    gaps: int


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
            aligned better.
        gaps: Number of internal gap runs in the best alignment (0 or 1).
        width: Number of columns of the candidate.
    """

    matrix_id: str
    name: str
    score: float
    percent_score: float
    is_reverse_complement: bool
    gaps: int
    width: int


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
) -> tuple[float, int]:
    """Best score and number of gap runs for one orientation."""
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

    best = float(prefix.max())
    gaps = 0

    # One gap run of length k: the alignment follows one diagonal up to a cell,
    # skips k columns of one profile, then follows the next diagonal.
    for k in range(1, max(n, m)):
        cost = open_penalty + (k - 1) * ext_penalty
        if n - k - 1 >= 1 and m >= 2:
            gapped = float((prefix[: n - k - 1, : m - 1] + suffix[k + 1 :, 1:]).max()) - cost
            if gapped > best:
                best, gaps = gapped, 1
        if m - k - 1 >= 1 and n >= 2:
            gapped = float((prefix[: n - 1, : m - k - 1] + suffix[1:, k + 1 :]).max()) - cost
            if gapped > best:
                best, gaps = gapped, 1
    return best, gaps


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
    further column. Both ``b`` and its reverse complement are tried.

    Args:
        a: First profile.
        b: Second profile.
        open_penalty: Penalty for the first column of a gap.
        ext_penalty: Penalty for each further column of the same gap.

    Returns:
        The best score, the orientation of ``b`` that gave it, and the number
        of gap runs used.

    Raises:
        ValueError: If a profile has a column with no counts.
    """
    fa, fb = _frequencies(a), _frequencies(b)
    forward, forward_gaps = _best_alignment(fa, fb, open_penalty, ext_penalty)
    reverse, reverse_gaps = _best_alignment(fa, fb[::-1, ::-1], open_penalty, ext_penalty)
    if reverse > forward:
        return AlignScore(reverse, True, reverse_gaps)
    return AlignScore(forward, False, forward_gaps)


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
        )
        for i, (c, r) in enumerate(zip(candidates, results, strict=True))
    ]
    hits.sort(key=lambda h: (-getattr(h, sort_by), h.matrix_id))
    return hits if top is None else hits[:top]
