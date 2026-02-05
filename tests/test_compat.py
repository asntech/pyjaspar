"""Backward compatibility tests.

These tests verify the existing public API contract. They must pass
after every refactoring commit to ensure nothing breaks for users.
"""

from __future__ import annotations

import warnings

from Bio.motifs.jaspar import Motif, Record

# --- Import compatibility ---


def test_import_jaspardb_class():
    """The lowercase jaspardb class must remain importable."""
    from pyjaspar import jaspardb as cls

    assert callable(cls)


def test_import_constants():
    from pyjaspar import (
        JASPAR_DFLT_COLLECTION,
        JASPAR_LATEST_RELEASE,
        __version__,
        jaspar_releases,
    )

    assert isinstance(__version__, str)
    assert isinstance(JASPAR_LATEST_RELEASE, str)
    assert isinstance(JASPAR_DFLT_COLLECTION, str)
    assert isinstance(jaspar_releases, dict)


def test_jaspar_releases_content():
    from pyjaspar import jaspar_releases

    expected = {
        "JASPAR2026",
        "JASPAR2024",
        "JASPAR2022",
        "JASPAR2020",
        "JASPAR2018",
        "JASPAR2016",
        "JASPAR2014",
    }
    assert set(jaspar_releases.keys()) == expected
    for _key, val in jaspar_releases.items():
        assert val.endswith(".sqlite")


# --- fetch_motif_by_id ---


def test_fetch_motif_by_id_full_id(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001.1")
    assert isinstance(motif, Motif)
    assert motif.matrix_id == "MA0001.1"
    assert motif.base_id == "MA0001"
    assert int(motif.version) == 1
    assert motif.name == "AGL3"


def test_fetch_motif_by_id_base_only_returns_latest(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001")
    assert isinstance(motif, Motif)
    assert motif.base_id == "MA0001"
    # Should be the latest version (>= 1)
    assert int(motif.version) >= 1


def test_fetch_motif_by_id_nonexistent(jdb_2026):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        motif = jdb_2026.fetch_motif_by_id("MA9999.1")
    assert motif is None


def test_motif_has_counts(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001.1")
    assert motif is not None
    counts = motif.counts
    assert "A" in counts
    assert "C" in counts
    assert "G" in counts
    assert "T" in counts
    assert len(counts["A"]) > 0
    # All rows same length
    assert len(counts["A"]) == len(counts["C"]) == len(counts["G"]) == len(counts["T"])


def test_motif_attributes(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0095.2")
    assert motif is not None
    assert motif.name == "YY1"
    assert motif.collection == "CORE"
    assert isinstance(motif.species, list)
    assert isinstance(motif.acc, list)
    assert isinstance(motif.tf_class, list)
    assert isinstance(motif.tf_family, list)


# --- fetch_motifs_by_name ---


def test_fetch_motifs_by_name_single(jdb_2026):
    motifs = jdb_2026.fetch_motifs_by_name("CTCF")
    assert isinstance(motifs, Record)
    assert len(motifs) > 0
    for m in motifs:
        assert m.name == "CTCF"


def test_fetch_motifs_by_name_returns_record(jdb_2026):
    motifs = jdb_2026.fetch_motifs_by_name("AGL3")
    assert isinstance(motifs, Record)


# --- fetch_motifs ---


def test_fetch_motifs_default_core(jdb_2026):
    """Default collection is CORE."""
    motifs = jdb_2026.fetch_motifs()
    assert isinstance(motifs, Record)
    assert len(motifs) > 0
    for m in motifs:
        assert m.collection == "CORE"


def test_fetch_motifs_tax_group(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group="vertebrates")
    assert len(motifs) > 0


def test_fetch_motifs_multiple_tax_groups(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group=["vertebrates", "insects"])
    assert len(motifs) > 0


def test_fetch_motifs_tf_class(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tf_class="C2H2 zinc finger factors")
    assert len(motifs) > 0


def test_fetch_motifs_min_ic(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group="vertebrates", min_ic=12)
    assert len(motifs) > 0
    for m in motifs:
        assert m.pssm.mean() >= 12


def test_fetch_motifs_min_length(jdb_2026):
    motifs = jdb_2026.fetch_motifs(collection="CORE", tax_group="vertebrates", min_length=10)
    assert len(motifs) > 0
    for m in motifs:
        assert m.length >= 10


def test_fetch_motifs_all_versions(jdb_2026):
    motifs_latest = jdb_2026.fetch_motifs(collection="CORE", tf_name="CTCF")
    motifs_all = jdb_2026.fetch_motifs(collection="CORE", tf_name="CTCF", all_versions=True)
    assert len(motifs_all) >= len(motifs_latest)


def test_fetch_motifs_all_flag(jdb_2026):
    motifs = jdb_2026.fetch_motifs(all=True)
    assert len(motifs) > 100  # Should be many motifs


# --- get_releases ---


def test_get_releases(jdb_2026):
    releases = jdb_2026.get_releases()
    assert isinstance(releases, list)
    assert len(releases) == 7
    assert "JASPAR2026" in releases
    assert "JASPAR2014" in releases


# --- Multi-release ---


def test_multi_release_same_motif(jdb_2020, jdb_2026):
    m_2020 = jdb_2020.fetch_motif_by_id("MA0001.1")
    m_2026 = jdb_2026.fetch_motif_by_id("MA0001.1")
    assert m_2020 is not None
    assert m_2026 is not None
    assert m_2020.matrix_id == m_2026.matrix_id
    assert m_2020.name == m_2026.name


# --- Return types ---


def test_return_type_motif(jdb_2026):
    motif = jdb_2026.fetch_motif_by_id("MA0001.1")
    assert isinstance(motif, Motif)


def test_return_type_record(jdb_2026):
    motifs = jdb_2026.fetch_motifs_by_name("CTCF")
    assert isinstance(motifs, Record)
