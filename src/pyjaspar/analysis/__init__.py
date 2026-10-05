"""Motif analysis tools for pyjaspar.

Requires the ``analysis`` extra: ``pip install pyjaspar[analysis]``
"""

from .alignment import AlignmentResult, align_motifs
from .enrichment import EnrichmentResult, motif_enrichment
from .matrix_align import AlignScore, ProfileHit, align_score, search_profiles
from .scanning import ScanHit, scan_sequence
from .similarity import (
    best_correlation,
    euclidean_distance,
    kl_divergence,
    pearson_correlation,
)

__all__ = [
    "best_correlation",
    "euclidean_distance",
    "kl_divergence",
    "pearson_correlation",
    "ScanHit",
    "scan_sequence",
    "AlignScore",
    "ProfileHit",
    "align_score",
    "search_profiles",
    "EnrichmentResult",
    "motif_enrichment",
    "AlignmentResult",
    "align_motifs",
]
