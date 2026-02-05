"""Tests for the pyjaspar CLI."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from click.testing import CliRunner

from pyjaspar.cli.main import cli


@staticmethod
def runner():
    return CliRunner()


def test_help():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "pyJASPAR" in result.output


def test_version():
    result = CliRunner().invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert "pyJASPAR" in result.output


# --- releases ---


def test_releases():
    result = CliRunner().invoke(cli, ["releases"])
    assert result.exit_code == 0
    assert "JASPAR2026" in result.output
    assert "JASPAR2014" in result.output


def test_releases_latest():
    result = CliRunner().invoke(cli, ["releases", "--latest"])
    assert result.exit_code == 0
    assert "JASPAR2026" in result.output


# --- motif-by-id ---


def test_motif_by_id():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA0001.1"])
    assert result.exit_code == 0
    assert "MA0001.1" in result.output
    assert "AGL3" in result.output


def test_motif_by_id_jaspar_format():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA0001.1", "--motif-format", "jaspar"])
    assert result.exit_code == 0
    assert "MA0001.1" in result.output


def test_motif_by_id_transfac_format():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA0001.1", "--motif-format", "transfac"])
    assert result.exit_code == 0


def test_motif_by_id_metadata_tsv():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA0001.1", "--metadata", "--format", "tsv"])
    assert result.exit_code == 0
    assert "matrix_id" in result.output
    assert "MA0001.1" in result.output


def test_motif_by_id_metadata_json():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA0001.1", "--metadata", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["count"] == 1
    assert data["results"][0]["matrix_id"] == "MA0001.1"


def test_motif_by_id_nonexistent():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA9999.1"])
    assert result.exit_code != 0


def test_motif_by_id_release():
    result = CliRunner().invoke(cli, ["motif-by-id", "MA0001.1", "-r", "2020"])
    assert result.exit_code == 0
    assert "MA0001.1" in result.output


# --- motifs-by-name ---


def test_motifs_by_name():
    result = CliRunner().invoke(cli, ["motifs-by-name", "CTCF"])
    assert result.exit_code == 0
    assert "CTCF" in result.output


def test_motifs_by_name_metadata():
    result = CliRunner().invoke(cli, ["motifs-by-name", "CTCF", "--metadata"])
    assert result.exit_code == 0
    assert "matrix_id" in result.output


# --- motifs ---


def test_motifs_default():
    result = CliRunner().invoke(cli, ["motifs", "--tax-group", "Vertebrates", "--min-length", "10"])
    assert result.exit_code == 0
    assert len(result.output) > 0


# --- metadata ---


def test_metadata():
    result = CliRunner().invoke(cli, ["metadata", "--tax-group", "Vertebrates"])
    assert result.exit_code == 0
    assert "matrix_id" in result.output


def test_metadata_json():
    result = CliRunner().invoke(cli, ["metadata", "--tax-group", "Vertebrates", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "count" in data
    assert "results" in data


def test_metadata_output_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = str(Path(tmpdir) / "output.tsv")
        result = CliRunner().invoke(cli, ["metadata", "--tax-group", "Vertebrates", "-o", outpath])
        assert result.exit_code == 0
        content = Path(outpath).read_text()
        assert "matrix_id" in content


# --- cite ---


def test_cite():
    result = CliRunner().invoke(cli, ["cite"])
    assert result.exit_code == 0
    assert "JASPAR 2024" in result.output
