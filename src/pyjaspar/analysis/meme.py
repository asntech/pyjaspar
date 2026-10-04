"""Motif comparison with p-values, through memesuite-lite.

``tomtom`` ranks candidate profiles by how significant their similarity to a query
motif is, with the Tomtom algorithm of the MEME suite as implemented by
memesuite-lite (https://github.com/jmschrei/memesuite-lite). memesuite-lite is not
installed with pyjaspar: ``pip install pyjaspar[meme]``.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif

_BASES = "ACGT"


@dataclass(frozen=True)
class TomtomHit:
    """A candidate profile compared with a query by Tomtom.

    Attributes:
        matrix_id: JASPAR matrix ID of the candidate.
        name: Name of the TF the candidate belongs to.
        pvalue: p-value of the best alignment between the query and the candidate.
        offset: Position in the query minus position in the candidate at the first
            aligned column, in the orientation of the candidate that was aligned (the
            sign of the offset of ``best_correlation``).
        overlap: Number of columns that the query and the candidate share.
        is_reverse_complement: True if the reverse complement of the candidate
            aligned better.
    """

    matrix_id: str
    name: str
    pvalue: float
    offset: int
    overlap: int
    is_reverse_complement: bool


def _probabilities(motif: Motif) -> np.ndarray:
    """Return the motif's columns as probabilities, shape (4, width)."""
    counts = np.array([motif.counts[b] for b in _BASES], dtype=float)
    totals = counts.sum(axis=0)
    if (totals <= 0).any():
        raise ValueError(f"Motif {motif.matrix_id or motif.name} has a column with no counts")
    return counts / totals


def tomtom(
    query: Motif,
    candidates: Iterable[Motif],
    top: int | None = None,
) -> list[TomtomHit]:
    """Compare a query motif with candidate profiles and rank them by p-value.

    Each column is turned into probabilities and the comparison is done by
    memesuite-lite. The p-values are computed against a background made of the
    columns of the candidates, so they depend on which candidates are given: use a
    large set, such as a whole collection (memesuite-lite warns below 25
    candidates). The first call on a machine compiles code and takes about half a
    minute; the compiled code is cached, so later calls take about a second.

    Build the candidates with ``JasparDB.fetch_motifs``.

    Args:
        query: The motif to search with.
        candidates: Profiles to compare the query with.
        top: Return only this many hits; all of them if None.

    Returns:
        The candidates, smallest p-value first.

    Raises:
        ImportError: If memesuite-lite is not installed.
        ValueError: If a profile has a column with no counts.
    """
    try:
        from memelite import tomtom as _tomtom
    except ImportError:
        raise ImportError(
            "tomtom() requires memesuite-lite. Install with: pip install pyjaspar[meme]"
        ) from None

    candidates = list(candidates)
    if not candidates:
        return []

    pvalues, _, offsets, overlaps, strands = _tomtom(
        [_probabilities(query)], [_probabilities(c) for c in candidates]
    )

    # memesuite-lite counts the offset as the position in the candidate minus the one in
    # the query; it is reported here the other way round, with the sign of the offset of
    # best_correlation.
    hits = [
        TomtomHit(
            matrix_id=c.matrix_id,
            name=c.name,
            pvalue=float(pvalues[0, i]),
            offset=-int(offsets[0, i]),
            overlap=int(overlaps[0, i]),
            is_reverse_complement=bool(strands[0, i]),
        )
        for i, c in enumerate(candidates)
    ]
    hits.sort(key=lambda h: (h.pvalue, h.matrix_id))
    return hits if top is None else hits[:top]
