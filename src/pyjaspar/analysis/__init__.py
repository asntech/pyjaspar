"""Motif analysis tools for pyjaspar.

Requires the ``analysis`` extra: ``pip install pyjaspar[analysis]``
"""

from .alignment import align_motifs, format_alignment
from .enrichment import EnrichmentResult, motif_enrichment
from .inference import InferenceHit, infer_profiles
from .methods import (
    AlignmentPath,
    AlignmentResult,
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
    "InferenceHit",
    "infer_profiles",
    "ScanHit",
    "scan_sequence",
    "align_motifs",
    "format_alignment",
    "search_profiles",
    "AlignmentPath",
    "AlignmentResult",
    "ProfileHit",
    "EnrichmentResult",
    "motif_enrichment",
]
