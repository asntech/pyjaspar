"""CLI smoke tests for DL commands."""

from __future__ import annotations

import json

from click.testing import CliRunner

from pyjaspar.cli.main import cli


def test_dl_help():
    result = CliRunner().invoke(cli, ["dl", "--help"])
    assert result.exit_code == 0
    assert "profile" in result.output
    assert "model" in result.output


def test_dl_profile_json():
    result = CliRunner().invoke(cli, ["dl", "profile", "DL0001.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["profile_id"] == "DL0001.1"
    assert data["tf_name"] == "ARNT2"
    assert data["primary_motif"] is not None


def test_dl_profile_tsv():
    result = CliRunner().invoke(cli, ["dl", "profile", "DL0001.1", "--format", "tsv"])
    assert result.exit_code == 0
    assert "ARNT2" in result.output


def test_dl_profile_not_found():
    result = CliRunner().invoke(cli, ["dl", "profile", "DL9999.1"])
    assert result.exit_code == 1


def test_dl_model_json():
    result = CliRunner().invoke(cli, ["dl", "model", "BP000001.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["model_id"] == "BP000001.1"
    assert data["tf_name"] == "REST"


def test_dl_model_tsv():
    result = CliRunner().invoke(cli, ["dl", "model", "BP000001.1", "--format", "tsv"])
    assert result.exit_code == 0
    assert "REST" in result.output


def test_dl_motif_json():
    result = CliRunner().invoke(cli, ["dl", "motif", "MO000001.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["motif_id"] == "MO000001.1"
    assert "PFM" in data["matrices"]


def test_dl_motif_tsv():
    result = CliRunner().invoke(cli, ["dl", "motif", "MO000001.1", "--format", "tsv"])
    assert result.exit_code == 0
    assert "MO000001.1" in result.output


def test_dl_profiles_search():
    result = CliRunner().invoke(cli, ["dl", "profiles", "--tf-name", "ARNT2"])
    assert result.exit_code == 0
    assert "ARNT2" in result.output


def test_dl_profiles_search_json():
    result = CliRunner().invoke(cli, ["dl", "profiles", "--tf-name", "ARNT2", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["count"] >= 1


def test_dl_models_search():
    result = CliRunner().invoke(cli, ["dl", "models", "--tf-name", "REST"])
    assert result.exit_code == 0
    assert "REST" in result.output


def test_dl_models_search_json():
    result = CliRunner().invoke(cli, ["dl", "models", "--tf-name", "REST", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["count"] >= 1


def test_dl_models_search_cell_line():
    result = CliRunner().invoke(cli, ["dl", "models", "--cell-line", "K562"])
    assert result.exit_code == 0
    assert "K562" in result.output
