"""Tests for DL availability guards."""

from __future__ import annotations

import pytest

from pyjaspar import DLNotAvailableError, JasparDB


def test_dl_not_available_on_old_release():
    db = JasparDB(release="JASPAR2020")
    with pytest.raises(DLNotAvailableError, match="not available"):
        _ = db.dl
    db.close()


def test_dl_available_on_2026():
    db = JasparDB(release="JASPAR2026")
    assert db.dl is not None
    db.close()


def test_dl_client_is_cached():
    db = JasparDB(release="JASPAR2026")
    dl1 = db.dl
    dl2 = db.dl
    assert dl1 is dl2
    db.close()


def test_dl_not_available_error_is_pyjaspar_error():
    """DLNotAvailableError should be a subclass of PyJasparError."""
    from pyjaspar import PyJasparError

    assert issubclass(DLNotAvailableError, PyJasparError)
