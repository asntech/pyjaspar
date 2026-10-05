"""Search a set of profiles with a query motif.

``search_profiles`` compares a query motif with candidate profiles using one of
the methods in ``methods.py`` and returns the candidates best first. With the
default ``"matrix_align"`` this reproduces JASPAR's "Matrix Align" web tool,
including its "Percent Score" column.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

from .methods import ProfileHit, _compare_many, _resolve_options, _sort_key

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif


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
    hits = _compare_many(query, items, method, options)
    hits.sort(key=lambda hit: (key(hit), hit.matrix_id or ""))
    return hits if top is None else hits[:top]
