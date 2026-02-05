"""Constants for the pyjaspar package."""

from __future__ import annotations

__version__ = "5.0.0"

jaspar_releases: dict[str, str] = {
    "JASPAR2026": "JASPAR2026.sqlite",
    "JASPAR2024": "JASPAR2024.sqlite",
    "JASPAR2022": "JASPAR2022.sqlite",
    "JASPAR2020": "JASPAR2020.sqlite",
    "JASPAR2018": "JASPAR2018.sqlite",
    "JASPAR2016": "JASPAR2016.sqlite",
    "JASPAR2014": "JASPAR2014.sqlite",
}

JASPAR_LATEST_RELEASE: str = "JASPAR2026"

JASPAR_DFLT_COLLECTION: str = "CORE"

JASPAR_DL_RELEASES: frozenset[str] = frozenset({"JASPAR2026"})
