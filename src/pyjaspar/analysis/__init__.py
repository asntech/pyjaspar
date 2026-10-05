"""Motif analysis tools for pyjaspar.

Requires the ``analysis`` extra: ``pip install pyjaspar[analysis]``
"""

from .enrichment import EnrichmentResult, motif_enrichment
from .inference import InferenceHit, infer_profiles
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
    "InferenceHit",
    "infer_profiles",
    "ScanHit",
    "scan_sequence",
    "EnrichmentResult",
    "motif_enrichment",
]
