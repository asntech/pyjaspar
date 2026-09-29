"""Profile inference: predict TF binding profiles from a protein sequence.

The search runs on a JASPAR server, so it needs network access and takes
several seconds per sequence.
"""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass

# The JASPAR release hosts that answer /api/v1/infer/. The current release
# host returns an empty result for every sequence, so it is not listed.
INFER_HOSTS = {
    "JASPAR2024": "https://jaspar2024.elixir.no",
    "JASPAR2022": "https://jaspar2022.genereg.net",
    "JASPAR2020": "https://jaspar2020.genereg.net",
}
DEFAULT_INFER_RELEASE = "JASPAR2024"


@dataclass
class InferenceHit:
    """A JASPAR profile predicted for a protein sequence.

    Attributes:
        matrix_id: JASPAR matrix ID, e.g. ``MA0162.2``.
        name: Name of the TF the profile belongs to.
        evalue: E-value of the DNA-binding domain match; lower is more significant.
        dbd_identity: Share of identical amino acids in the DNA-binding domain,
            as a fraction (0-1).
        logo_url: URL of the profile's sequence logo (SVG).
        release: JASPAR release the matrix ID belongs to.
    """

    matrix_id: str
    name: str
    evalue: float
    dbd_identity: float
    logo_url: str
    release: str


def infer_profiles(
    sequence: str,
    release: str = DEFAULT_INFER_RELEASE,
    timeout: float = 120.0,
) -> list[InferenceHit]:
    """Predict the JASPAR profiles bound by a protein.

    The protein is compared with the DNA-binding domains of TFs that already
    have a JASPAR profile, so the hits are its closest relatives.

    Args:
        sequence: Amino acid sequence, full length or the DNA-binding domain.
        release: JASPAR release to search; one of ``INFER_HOSTS``.
        timeout: Seconds to wait for the server.

    Returns:
        Predicted profiles, best match first. Empty if nothing matches.

    Raises:
        ValueError: If the sequence is not a protein sequence or the release
            has no inference service.
        urllib.error.URLError: If the server cannot be reached.
    """
    host = INFER_HOSTS.get(release)
    if host is None:
        raise ValueError(
            f"Inference is not available for {release}; use one of {list(INFER_HOSTS)}"
        )

    sequence = "".join(sequence.split()).rstrip("*").upper()
    if not re.fullmatch(r"[A-Z]+", sequence):
        raise ValueError("The sequence must contain only amino acid letters")

    url = f"{host}/api/v1/infer/{urllib.parse.quote(sequence)}/?format=json"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = json.load(response)

    return [
        InferenceHit(
            matrix_id=r["matrix_id"],
            name=r["name"],
            evalue=r["evalue"],
            dbd_identity=r["dbd"],
            logo_url=r["sequence_logo"],
            release=release,
        )
        for r in data["results"] or []
    ]
