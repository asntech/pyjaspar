"""pyJASPAR: A Pythonic interface to JASPAR transcription factor motifs.

Example::

    from pyjaspar import JasparDB

    jdb = JasparDB(release='JASPAR2026')
    motif = jdb.fetch_motif_by_id('MA0095.2')
    print(motif)
"""

from ._constants import (
    JASPAR_DFLT_COLLECTION,
    JASPAR_LATEST_RELEASE,
    __version__,
    jaspar_releases,
)
from ._exceptions import (
    DatabaseError,
    DLNotAvailableError,
    MotifNotFoundError,
    PyJasparError,
    ReleaseNotFoundError,
)
from .db import JasparDB

# Backward compatibility alias — existing code using ``jaspardb`` continues to work.
jaspardb = JasparDB

__all__ = [
    "__version__",
    "JasparDB",
    "jaspardb",
    "JASPAR_DFLT_COLLECTION",
    "JASPAR_LATEST_RELEASE",
    "jaspar_releases",
    "DLNotAvailableError",
    "DatabaseError",
    "MotifNotFoundError",
    "PyJasparError",
    "ReleaseNotFoundError",
]
