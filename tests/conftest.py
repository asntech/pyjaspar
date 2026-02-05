"""Shared fixtures for pyjaspar tests."""

from __future__ import annotations

import pytest

from pyjaspar import jaspardb


@pytest.fixture(scope="session")
def jdb_2026():
    """A jaspardb connection to JASPAR2026 reused across all tests."""
    db = jaspardb(release="JASPAR2026")
    yield db


@pytest.fixture(scope="session")
def jdb_2020():
    """A jaspardb connection to JASPAR2020 for cross-release tests."""
    db = jaspardb(release="JASPAR2020")
    yield db


@pytest.fixture(scope="session")
def sample_motif(jdb_2026):
    """A known motif (AGL3) for deterministic assertions."""
    return jdb_2026.fetch_motif_by_id("MA0001.1")
