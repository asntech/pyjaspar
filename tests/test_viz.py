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
