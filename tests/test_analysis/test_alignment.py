"""Tests for motif alignment."""

from __future__ import annotations

import numpy as np
import pytest
from Bio.Align import Alignment

from pyjaspar import JasparDB
from pyjaspar.analysis.alignment import _offset_to_coordinates, align_motifs
from pyjaspar.analysis.similarity import pearson_correlation


@pytest.fixture(scope="module")
def jdb():
    return JasparDB()


@pytest.fixture(scope="module")
def motif_agl3(jdb):
    return jdb.fetch_motif_by_id("MA0001.1")


@pytest.fixture(scope="module")
def motif_runx1(jdb):
    return jdb.fetch_motif_by_id("MA0002.1")


# --- _offset_to_coordinates ---


def test_offset_to_coordinates_negative():
    """offset < 0: sequence 2 leads, sequence 1 trails."""
    coordinates = _offset_to_coordinates(len1=10, len2=11, offset=-6)
    assert coordinates.tolist() == [[0, 0, 5, 10], [0, 6, 11, 11]]
    alignment = Alignment(["A" * 10, "C" * 11], coordinates)
    assert str(alignment)  # renders without error


def test_offset_to_coordinates_zero():
    """offset == 0: no leading gap, only a trailing gap for the shorter sequence."""
    coordinates = _offset_to_coordinates(len1=10, len2=11, offset=0)
    alignment = Alignment(["A" * 10, "C" * 11], coordinates)
    assert str(alignment)


def test_offset_to_coordinates_positive():
    """offset > 0: sequence 1 leads, sequence 2 trails."""
    coordinates = _offset_to_coordinates(len1=10, len2=11, offset=3)
    alignment = Alignment(["A" * 10, "C" * 11], coordinates)
    assert str(alignment)


# --- align_motifs ---


def test_align_motifs_self(motif_agl3):
    """Self-alignment should have score 1.0, offset 0, no reverse complement."""
    result = align_motifs(motif_agl3, motif_agl3)
    assert result.score > 0.99
    assert result.offset == 0
    assert result.is_reverse_complement is False
    assert isinstance(result.alignment, Alignment)


def test_align_motifs_different(motif_agl3, motif_runx1):
    result = align_motifs(motif_agl3, motif_runx1)
    assert isinstance(result.score, float)
    assert isinstance(result.offset, int)
    assert isinstance(result.is_reverse_complement, bool)
    assert result.motif1_id == motif_agl3.matrix_id
    assert result.motif2_id == motif_runx1.matrix_id


def test_align_motifs_score_matches_best_correlation(motif_agl3, motif_runx1):
    """The reported score must be reproducible from pearson_correlation directly.

    Guards against the offset/coordinate math silently drifting from the
    scorer it's built on top of.
    """
    result = align_motifs(motif_agl3, motif_runx1)
    motif2_used = motif_runx1.reverse_complement() if result.is_reverse_complement else motif_runx1
    recomputed = pearson_correlation(motif_agl3, motif2_used, result.offset)
    assert recomputed == pytest.approx(result.score)


def test_align_motifs_overlap_width(motif_agl3, motif_runx1):
    """The aligned overlap width must match the offset-derived expectation."""
    result = align_motifs(motif_agl3, motif_runx1)
    motif2_used = motif_runx1.reverse_complement() if result.is_reverse_complement else motif_runx1
    len1, len2 = motif_agl3.length, motif2_used.length
    expected_overlap = min(len1 - max(0, result.offset), len2 - max(0, -result.offset))

    coordinates = np.asarray(result.alignment.coordinates)
    row1, row2 = coordinates
    overlap = sum(
        (row1[i + 1] - row1[i]) == (row2[i + 1] - row2[i]) and (row1[i + 1] - row1[i]) > 0
        for i in range(len(row1) - 1)
    )
    aligned_overlap = sum(
        (row1[i + 1] - row1[i])
        for i in range(len(row1) - 1)
        if (row1[i + 1] - row1[i]) == (row2[i + 1] - row2[i]) and (row1[i + 1] - row1[i]) > 0
    )
    assert overlap >= 1
    assert aligned_overlap == expected_overlap
