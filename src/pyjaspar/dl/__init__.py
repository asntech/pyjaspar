"""JASPAR Deep Learning (DL) collection support.

The DL collection was introduced in JASPAR 2026 and contains motif patterns
derived from BPNet deep learning models trained on ENCODE TF ChIP-seq data.

Requires the ``dl`` extra: ``pip install pyjaspar[dl]``
"""

from ._client import DLClient
from ._models import (
    DLModel,
    DLModelSummary,
    DLMotifPattern,
    DLProfile,
    DLProfileSummary,
    Matrix,
)

__all__ = [
    "DLClient",
    "DLModel",
    "DLModelSummary",
    "DLMotifPattern",
    "DLProfile",
    "DLProfileSummary",
    "Matrix",
]
