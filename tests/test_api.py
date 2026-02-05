"""Smoke tests for the pyjaspar public API."""

from __future__ import annotations

from Bio.motifs.jaspar import Motif, Record


def test_fetch_motif_by_id(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001.1")
    assert motif is not None
    assert isinstance(motif, Motif)
    assert motif.matrix_id == "MA0001.1"
    assert motif.name == "AGL3"


def test_fetch_motif_by_id_base_only(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001")
    assert motif is not None
    assert motif.base_id == "MA0001"


def test_fetch_motifs_by_name(jdb_2026):
    motifs = jdb_2026.fetch_motifs_by_name("CTCF")
    assert isinstance(motifs, Record)
    assert len(motifs) > 0
    for m in motifs:
        assert m.name == "CTCF"


def test_fetch_motifs_core_vertebrates(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group="vertebrates")
    assert isinstance(motifs, Record)
    assert len(motifs) > 0


def test_get_releases(jdb_2026):
    releases = jdb_2026.get_releases()
    assert isinstance(releases, list)
    assert "JASPAR2026" in releases
    assert "JASPAR2020" in releases


def test_constants_importable():
    from pyjaspar import (
        JASPAR_DFLT_COLLECTION,
        JASPAR_LATEST_RELEASE,
        __version__,
        jaspar_releases,
    )

    assert __version__ == "5.0.0"
    assert JASPAR_LATEST_RELEASE == "JASPAR2026"
    assert JASPAR_DFLT_COLLECTION == "CORE"
    assert isinstance(jaspar_releases, dict)
    assert len(jaspar_releases) == 7


def test_cross_release(jdb_2020):
    motif = jdb_2020.fetch_motif_by_id("MA0001.1")
    assert motif is not None
    assert motif.matrix_id == "MA0001.1"
