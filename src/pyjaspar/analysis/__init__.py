"""Motif analysis tools for pyjaspar.

Requires the ``analysis`` extra: ``pip install pyjaspar[analysis]``
"""

from .alignment import align_motifs, format_alignment
from .enrichment import EnrichmentResult, motif_enrichment
from .methods import (
    AlignmentPath,
    AlignmentResult,
    MatrixAlignHit,
    MatrixAlignResult,
    PearsonHit,
    PearsonResult,
    ProfileHit,
)
from .profile_search import search_profiles
from .scanning import ScanHit, scan_sequence
from .similarity import (
    euclidean_distance,
    kl_divergence,
    pearson_correlation,
)

__all__ = [
    "euclidean_distance",
    "kl_divergence",
    "pearson_correlation",
    "ScanHit",
    "scan_sequence",
    "align_motifs",
    "format_alignment",
    "search_profiles",
    "AlignmentPath",
    "AlignmentResult",
    "MatrixAlignResult",
    "PearsonResult",
    "ProfileHit",
    "MatrixAlignHit",
    "PearsonHit",
    "EnrichmentResult",
    "motif_enrichment",
]
