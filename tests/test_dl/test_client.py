"""Integration tests for DL client against real JASPAR2026 database."""

from __future__ import annotations

import pytest

from pyjaspar._exceptions import MotifNotFoundError
from pyjaspar.dl._models import (
    DLModel,
    DLModelSummary,
    DLMotifPattern,
    DLProfile,
    DLProfileSummary,
    Matrix,
)


def test_fetch_profile(dl_client):
    profile = dl_client.fetch_profile("DL0001.1")
    assert isinstance(profile, DLProfile)
    assert profile.profile_id == "DL0001.1"
    assert profile.tf_name == "ARNT2"
    assert profile.tax_group == "vertebrates"


def test_fetch_profile_primary_motif(dl_client):
    profile = dl_client.fetch_profile("DL0001.1")
    assert profile.primary_motif is not None
    assert profile.primary_motif.motif_type == "pmp"
    assert "PFM" in profile.primary_motif.matrices
    assert profile.primary_motif.pfm.length > 0


def test_fetch_profile_latest_version(dl_client):
    profile = dl_client.fetch_profile("DL0001")
    assert profile.base_id == "DL0001"
    assert profile.version >= 1


def test_fetch_profile_not_found(dl_client):
    with pytest.raises(MotifNotFoundError):
        dl_client.fetch_profile("DL9999.1")


def test_fetch_profile_has_alt_motifs(dl_client):
    """DL0002 (ATF2) has PMP + 2 AMPs."""
    profile = dl_client.fetch_profile("DL0002.1")
    assert len(profile.alt_motifs) == 2
    for m in profile.alt_motifs:
        assert m.motif_type == "amp"


def test_fetch_profile_jaspar_matches(dl_client):
    profile = dl_client.fetch_profile("DL0001.1")
    assert len(profile.jaspar_matches) >= 1
    for _cluster_id, match_id in profile.jaspar_matches.items():
        assert match_id.startswith("MA")


def test_fetch_profile_linked_models(dl_client):
    profile = dl_client.fetch_profile("DL0001.1")
    assert len(profile.linked_model_ids) >= 1
    for mid in profile.linked_model_ids:
        assert mid.startswith("BP")


def test_fetch_model(dl_client):
    model = dl_client.fetch_model("BP000001.1")
    assert isinstance(model, DLModel)
    assert model.model_id == "BP000001.1"
    assert model.tf_name == "REST"
    assert model.cell_line == "K562"
    assert model.tax_id == 9606
    assert model.data_type == "ChIP-seq"


def test_fetch_model_motifs(dl_client):
    model = dl_client.fetch_model("BP000001.1")
    assert len(model.motifs) > 0
    assert model.primary_motif is not None
    assert model.primary_motif.motif_type == "pmp"


def test_fetch_model_latest_version(dl_client):
    model = dl_client.fetch_model("BP000001")
    assert model.base_id == "BP000001"
    assert model.version >= 1


def test_fetch_model_not_found(dl_client):
    with pytest.raises(MotifNotFoundError):
        dl_client.fetch_model("BP999999.1")


def test_fetch_motif_pattern(dl_client):
    motif = dl_client.fetch_motif_pattern("MO000001.1")
    assert isinstance(motif, DLMotifPattern)
    assert motif.motif_id == "MO000001.1"
    assert motif.tf_name == "ARNT2"
    assert motif.motif_type == "pmp"


def test_fetch_motif_pattern_matrices(dl_client):
    motif = dl_client.fetch_motif_pattern("MO000001.1")
    assert motif.pfm is not None
    assert isinstance(motif.pfm, Matrix)
    assert motif.pfm.kind == "PFM"
    assert motif.pfm.length > 0
    assert motif.cwm is not None
    assert motif.cwm.kind == "CWM"


def test_fetch_motif_pattern_jaspar_match(dl_client):
    motif = dl_client.fetch_motif_pattern("MO000001.1")
    assert motif.jaspar_match == "MA1464.2"


def test_fetch_motif_pattern_not_found(dl_client):
    with pytest.raises(MotifNotFoundError):
        dl_client.fetch_motif_pattern("MO999999.1")


def test_motif_to_biopython(dl_client):
    motif = dl_client.fetch_motif_pattern("MO000001.1")
    bp_motif = motif.to_biopython()
    assert bp_motif.matrix_id == "MO000001.1"
    assert bp_motif.name == "ARNT2"
    assert "A" in bp_motif.counts


def test_search_profiles_by_name(dl_client):
    results = dl_client.search_profiles(tf_name="ARNT2")
    assert len(results) >= 1
    assert all(isinstance(r, DLProfileSummary) for r in results)
    assert all(r.tf_name == "ARNT2" for r in results)


def test_search_profiles_by_tax_group(dl_client):
    results = dl_client.search_profiles(tax_group="vertebrates")
    assert len(results) == 240


def test_search_profiles_no_results(dl_client):
    results = dl_client.search_profiles(tf_name="NONEXISTENT_TF_XYZ")
    assert len(results) == 0


def test_search_models_by_name(dl_client):
    results = dl_client.search_models(tf_name="REST")
    assert len(results) >= 1
    assert all(isinstance(r, DLModelSummary) for r in results)
    assert all(r.tf_name == "REST" for r in results)


def test_search_models_by_cell_line(dl_client):
    results = dl_client.search_models(cell_line="K562")
    assert len(results) > 0
    assert all(r.cell_line == "K562" for r in results)


def test_search_models_combined(dl_client):
    results = dl_client.search_models(tf_name="REST", cell_line="K562")
    assert len(results) >= 1
    assert all(r.tf_name == "REST" and r.cell_line == "K562" for r in results)
