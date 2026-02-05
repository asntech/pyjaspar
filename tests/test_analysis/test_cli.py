"""CLI smoke tests for analysis commands."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from click.testing import CliRunner

from pyjaspar.cli.main import cli


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    """Write a simple FASTA file."""
    with open(path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")


def test_analysis_commands_in_help():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "scan" in result.output
    assert "similarity" in result.output
    assert "enrichment" in result.output


# --- scan tests ---


def test_scan_literal():
    result = CliRunner().invoke(
        cli,
        [
            "scan",
            "ACGTACGTACGTACGT" * 5,
            "--motif-id",
            "MA0001.1",
            "--threshold",
            "0.3",
        ],
    )
    assert result.exit_code == 0


def test_scan_json():
    result = CliRunner().invoke(
        cli,
        [
            "scan",
            "ACGTACGTACGTACGT" * 5,
            "--motif-id",
            "MA0001.1",
            "--threshold",
            "0.3",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "results" in data
    assert "count" in data
    assert data["motif_id"] == "MA0001.1"


def test_scan_fasta_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "test.fasta"
        _write_fasta(
            fasta_path,
            {
                "seq1": "ACGTACGTACGTACGT" * 5,
                "seq2": "TGCATGCATGCATGCA" * 5,
            },
        )
        result = CliRunner().invoke(
            cli,
            [
                "scan",
                str(fasta_path),
                "--motif-id",
                "MA0001.1",
                "--threshold",
                "0.3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "results" in data


def test_scan_invalid_motif():
    result = CliRunner().invoke(
        cli,
        ["scan", "ACGTACGT", "--motif-id", "MA9999.1"],
    )
    assert result.exit_code == 1


def test_scan_no_reverse():
    result = CliRunner().invoke(
        cli,
        [
            "scan",
            "ACGTACGTACGTACGT" * 5,
            "--motif-id",
            "MA0001.1",
            "--threshold",
            "0.3",
            "--no-reverse",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0


# --- similarity tests ---


def test_similarity_pearson():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1"])
    assert result.exit_code == 0


def test_similarity_euclidean():
    result = CliRunner().invoke(
        cli, ["similarity", "MA0001.1", "MA0002.1", "--metric", "euclidean"]
    )
    assert result.exit_code == 0


def test_similarity_kl():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1", "--metric", "kl"])
    assert result.exit_code == 0


def test_similarity_best():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1", "--best"])
    assert result.exit_code == 0


def test_similarity_json():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "score" in data
    assert data["motif1"] == "MA0001.1"
    assert data["motif2"] == "MA0002.1"


def test_similarity_best_json():
    result = CliRunner().invoke(
        cli, ["similarity", "MA0001.1", "MA0002.1", "--best", "--format", "json"]
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "best_score" in data
    assert "best_offset" in data
    assert "is_reverse_complement" in data


def test_similarity_invalid_motif():
    result = CliRunner().invoke(cli, ["similarity", "MA9999.1", "MA0002.1"])
    assert result.exit_code == 1


def test_similarity_with_offset():
    result = CliRunner().invoke(
        cli, ["similarity", "MA0001.1", "MA0002.1", "--offset", "2", "--format", "json"]
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["offset"] == 2


# --- enrichment tests ---


def test_enrichment_basic():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {f"fg{i}": "ACGTACGTACGTACGT" * 5 for i in range(5)})
        _write_fasta(bg_path, {f"bg{i}": "TGCATGCATGCATGCA" * 5 for i in range(5)})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA0001.1",
                "--threshold",
                "0.3",
            ],
        )
        assert result.exit_code == 0


def test_enrichment_json():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {f"fg{i}": "ACGTACGTACGTACGT" * 5 for i in range(5)})
        _write_fasta(bg_path, {f"bg{i}": "TGCATGCATGCATGCA" * 5 for i in range(5)})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA0001.1",
                "--threshold",
                "0.3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "results" in data
        assert "count" in data


def test_enrichment_multiple_motifs():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {f"fg{i}": "ACGTACGTACGTACGT" * 5 for i in range(3)})
        _write_fasta(bg_path, {f"bg{i}": "TGCATGCATGCATGCA" * 5 for i in range(3)})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA0001.1,MA0002.1",
                "--threshold",
                "0.3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["count"] == 2


def test_enrichment_invalid_motif_warns():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {"fg1": "ACGTACGTACGTACGT" * 5})
        _write_fasta(bg_path, {"bg1": "TGCATGCATGCATGCA" * 5})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA9999.1",
            ],
        )
        assert result.exit_code == 1
