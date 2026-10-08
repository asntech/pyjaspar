"""Search a set of profiles with a query motif.

``search_profiles`` compares a query motif with candidate profiles using one of
the methods in ``methods.py`` and returns the candidates best first. With the
default ``"matrix_align"`` this reproduces JASPAR's "Matrix Align" web tool,
including its "Percent Score" column.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable, Sequence
from typing import TYPE_CHECKING, Any

from .methods import (
    AlignmentResult,
    ProfileHit,
    _check_method,
    _compare,
    _matrix_id,
    _resolve_options,
    _validate_profile,
)

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif


_MATRIX_ID = re.compile(r"[A-Z]+(\d+)\.(\d+)")


def search_profiles(
    query: Motif,
    candidates: Iterable[Motif],
    *,
    method: str = "matrix_align",
    both_strands: bool = True,
    min_overlap: int | None = None,
    open_penalty: float | None = None,
    ext_penalty: float | None = None,
    sort_by: str | None = None,
    top: int | None = None,
) -> list[ProfileHit]:
    """Compare a query motif with candidate profiles and rank them, best first.

    Build the candidates with ``JasparDB.fetch_motifs``, which offers the same
    filters as the web tool (collection, taxonomic group, latest or all versions).

    Args:
        query: The motif to search with.
        candidates: Profiles to compare the query with.
        method: ``"matrix_align"`` (default) or ``"pearson"``; see ``methods.py``.
        both_strands: If True, also try each candidate's reverse complement.
        min_overlap: ``"pearson"`` only: minimum overlapping columns (default 4).
        open_penalty: ``"matrix_align"`` only: cost of a gap's first column (default 3.0).
        ext_penalty: ``"matrix_align"`` only: cost of each further gap column (default 0.01).
        sort_by: ``"score"`` (default), or ``"percent_score"`` for ``"matrix_align"``.
        top: Return only this many hits; all of them if None.

    Returns:
        ``ProfileHit`` objects, best first. Ties are broken by
        matrix ID, then by input order.

    Raises:
        ValueError: If the method, an option, ``sort_by`` or ``top`` is invalid, or a
            profile is invalid (``"pearson"``: also shorter than ``min_overlap``).
    """
    key = _sort_key(method, sort_by)
    options = _resolve_options(
        method,
        both_strands=both_strands,
        min_overlap=min_overlap,
        open_penalty=open_penalty,
        ext_penalty=ext_penalty,
    )
    if top is not None and (isinstance(top, bool) or not isinstance(top, int) or top < 0):
        raise ValueError(f"top must be None or an integer >= 0, not {top!r}")
    items = list(candidates)
    if not items:
        return []
    _validate_profile(query)
    for candidate in items:
        _validate_profile(candidate)
    results = [_compare(query, c, method, options) for c in items]
    if method == "matrix_align":
        percents: list[float | None] = list(_percent_scores(query, items, results))
    else:
        percents = [None] * len(items)
    hits = [
        ProfileHit(_matrix_id(c), getattr(c, "name", None), r, p)
        for c, r, p in zip(items, results, percents, strict=True)
    ]
    hits.sort(key=lambda hit: (key(hit), hit.matrix_id or ""))
    return hits if top is None else hits[:top]


def _percent_scores(
    query: Motif, candidates: Sequence[Motif], results: Sequence[AlignmentResult]
) -> list[float]:
    """The web tool's Percent Score, computed over the full set before ranking or truncating.

    The web tool divides by the narrowest profile seen so far, going through its
    table in matrix-ID order; see ``ProfileHit.percent_score``.
    """
    narrowest = query.length
    percents = [0.0] * len(candidates)
    for i in _matrix_id_order(candidates):
        narrowest = min(narrowest, candidates[i].length)
        percents[i] = 100 * results[i].score / (2 * narrowest)
    return percents


def _matrix_id_order(candidates: Sequence[Motif]) -> list[int]:
    """Indices in matrix-ID order, or input order if any ID is not standard."""
    keys = []
    for candidate in candidates:
        match = _MATRIX_ID.fullmatch(_matrix_id(candidate) or "")
        if match is None:
            return list(range(len(candidates)))
        keys.append((int(match.group(1)), int(match.group(2))))
    return sorted(range(len(candidates)), key=keys.__getitem__)


def _sort_key(method: object, sort_by: str | None) -> Callable[[Any], float]:
    """Validate sorting before options or candidates; return a best-first key."""
    method = _check_method(method)
    if sort_by in (None, "score"):
        return lambda hit: -hit.score
    if method == "pearson":
        raise ValueError(f"sort_by must be 'score' for method='pearson', not {sort_by!r}")
    if sort_by == "percent_score":
        return lambda hit: -hit.percent_score
    raise ValueError(
        f"sort_by must be 'score' or 'percent_score' for method='matrix_align', not {sort_by!r}"
    )
