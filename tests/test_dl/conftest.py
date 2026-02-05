"""Fixtures for DL collection tests."""

from __future__ import annotations

import pytest

from pyjaspar import JasparDB


@pytest.fixture(scope="session")
def dl_client():
    """A DL sub-client connected to JASPAR2026."""
    db = JasparDB(release="JASPAR2026")
    yield db.dl
    db.close()
