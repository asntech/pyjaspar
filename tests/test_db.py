"""Tests for the JasparDB database layer."""

from __future__ import annotations

import warnings

import pytest
from Bio.motifs.jaspar import Motif, Record

from pyjaspar import JasparDB
from pyjaspar._exceptions import ReleaseNotFoundError


def test_constructor_default():
    db = JasparDB()
    assert db.release == "JASPAR2026"
    db.close()


def test_constructor_specific_release():
    db = JasparDB(release="JASPAR2020")
    assert db.release == "JASPAR2020"
    db.close()


def test_constructor_invalid_release():
    with pytest.raises(ReleaseNotFoundError, match="JASPAR9999"):
        JasparDB(release="JASPAR9999")


def test_constructor_invalid_path():
    with pytest.raises(FileNotFoundError):
        JasparDB(sqlite_db_path="/nonexistent/path.sqlite")


def test_context_manager():
    with JasparDB() as jdb:
        motif = jdb.fetch_motif_by_id("MA0001.1")
        assert motif is not None


def test_conn_property():
    """Backward compat: .conn attribute should still work."""
    db = JasparDB()
    assert db.conn is not None
    db.close()


def test_str():
    db = JasparDB()
    s = str(db)
    assert "JASPAR" in s
    db.close()


def test_fetch_motif_by_id_returns_motif(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001.1")
    assert isinstance(motif, Motif)


def test_fetch_motif_by_id_resolves_species_name(sample_motif):
    assert sample_motif.species == ["3702"]
    assert sample_motif.species_name == ["Arabidopsis thaliana"]


def test_fetch_motif_by_id_nonexistent_returns_none(jdb_2026):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        motif = jdb_2026.fetch_motif_by_id("MA9999.1")
    assert motif is None


def test_fetch_motifs_returns_record(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group="vertebrates")
    assert isinstance(motifs, Record)
    assert len(motifs) > 0


def test_fetch_motifs_data_type(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", data_type="ChIP-seq")
    assert len(motifs) > 0


def test_fetch_motifs_species(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", species=9606)
    assert len(motifs) > 0


def test_fetch_motifs_min_sites(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group="vertebrates", min_sites=100)
    for m in motifs:
        num_sites = sum(m.counts[nt][0] for nt in "ACGT")
        assert num_sites >= 100


def test_get_releases():
    db = JasparDB()
    releases = db.get_releases()
    assert len(releases) == 7
    assert releases[0] == "JASPAR2026"
    db.close()
