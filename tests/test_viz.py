"""Tests for motif visualization."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from pyjaspar import JasparDB
from pyjaspar.viz import plot_logo


@pytest.fixture(scope="module")
def motif_ctcf():
    return JasparDB().fetch_motif_by_id("MA0139.1")


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_plot_logo_returns_logo_with_one_row_per_position(motif_ctcf):
    logo = plot_logo(motif_ctcf)
    assert len(logo.df) == motif_ctcf.length
    assert list(logo.df.columns) == ["A", "C", "G", "T"]


def test_plot_logo_positions_are_one_based(motif_ctcf):
    logo = plot_logo(motif_ctcf)
    assert logo.df.index[0] == 1
    assert logo.df.index[-1] == motif_ctcf.length


def test_plot_logo_heights_within_two_bits(motif_ctcf):
    logo = plot_logo(motif_ctcf)
    assert logo.df.sum(axis=1).max() <= 2.0 + 1e-9


def test_plot_logo_axes_labels(motif_ctcf):
    logo = plot_logo(motif_ctcf)
    assert logo.ax.get_xlabel() == "Position"
    assert logo.ax.get_ylabel() == "Information content (bits)"
    assert logo.ax.get_ylim() == (0, 2)


def test_plot_logo_draws_on_given_axes(motif_ctcf):
    _, ax = plt.subplots()
    logo = plot_logo(motif_ctcf, ax=ax)
    assert logo.ax is ax


# --- logo CLI tests ---


def test_logo_cli_writes_png(tmp_path):
    from click.testing import CliRunner

    from pyjaspar.cli.main import cli

    out = tmp_path / "ctcf.png"
    result = CliRunner().invoke(cli, ["logo", "MA0139.1", "-o", str(out)])
    assert result.exit_code == 0
    assert out.read_bytes().startswith(b"\x89PNG")


def test_logo_cli_writes_svg(tmp_path):
    from click.testing import CliRunner

    from pyjaspar.cli.main import cli

    out = tmp_path / "ctcf.svg"
    result = CliRunner().invoke(cli, ["logo", "MA0139.1", "-o", str(out)])
    assert result.exit_code == 0
    assert "<svg" in out.read_text()


def test_logo_cli_invalid_motif(tmp_path):
    from click.testing import CliRunner

    from pyjaspar.cli.main import cli

    out = tmp_path / "x.png"
    result = CliRunner().invoke(cli, ["logo", "MA9999.9", "-o", str(out)])
    assert result.exit_code == 1
    assert not out.exists()


def test_logo_in_help():
    from click.testing import CliRunner

    from pyjaspar.cli.main import cli

    result = CliRunner().invoke(cli, ["--help"])
    assert "logo" in result.output
