"""Tests for pairwise motif comparison and its text display."""

from __future__ import annotations

import numpy as np
import pytest
from Bio.Align import Alignment
from Bio.motifs.jaspar import Motif

from pyjaspar import JasparDB
from pyjaspar.analysis import (
    AlignmentPath,
    MatrixAlignResult,
    PearsonResult,
    align_motifs,
    format_alignment,
)
from pyjaspar.analysis.alignment import _display_coordinates
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


def one_hot(sequence, matrix_id="X"):
    """A motif whose columns each have all counts on one base."""
    return Motif(
        matrix_id=matrix_id,
        name=matrix_id,
        counts={b: [10 if c == b else 0 for c in sequence] for b in "ACGT"},
    )


# --- align_motifs ---


def test_default_method_is_matrix_align(motif_agl3, motif_runx1):
    result = align_motifs(motif_agl3, motif_runx1)
    assert isinstance(result, MatrixAlignResult)
    assert result.method == "matrix_align"
    assert result.parameters == {"open_penalty": 3.0, "ext_penalty": 0.01, "both_strands": True}


def test_pearson_self(motif_agl3):
    """Self-alignment should have correlation 1.0, offset 0, no reverse complement."""
    result = align_motifs(motif_agl3, motif_agl3, method="pearson")
    assert isinstance(result, PearsonResult)
    assert result.correlation > 0.99
    assert result.score == result.correlation
    assert result.offset == 0
    assert result.is_reverse_complement is False
    assert isinstance(result.alignment, AlignmentPath)
    assert result.gaps == 0


def test_pearson_different(motif_agl3, motif_runx1):
    result = align_motifs(motif_agl3, motif_runx1, method="pearson")
    assert isinstance(result.correlation, float)
    assert isinstance(result.offset, int)
    assert isinstance(result.is_reverse_complement, bool)
    assert result.motif1_id == motif_agl3.matrix_id
    assert result.motif2_id == motif_runx1.matrix_id
    assert result.parameters == {"min_overlap": 4, "both_strands": True}


def test_pearson_score_is_reproducible_from_pearson_correlation(motif_agl3, motif_runx1):
    """The reported correlation is pearson_correlation at the reported offset and orientation."""
    result = align_motifs(motif_agl3, motif_runx1, method="pearson")
    motif2_used = motif_runx1.reverse_complement() if result.is_reverse_complement else motif_runx1
    recomputed = pearson_correlation(motif_agl3, motif2_used, result.offset)
    assert recomputed == pytest.approx(result.correlation)


def test_pearson_overlap_width(motif_agl3, motif_runx1):
    """The aligned overlap matches the width implied by the offset."""
    result = align_motifs(motif_agl3, motif_runx1, method="pearson")
    len1, len2 = motif_agl3.length, motif_runx1.length
    expected = min(len1 - max(0, result.offset), len2 - max(0, -result.offset))
    assert result.overlap == expected
    assert result.alignment_length == expected


# --- display ---


@pytest.mark.parametrize(
    ("offset", "expected"),
    [
        # the coordinates #11 rendered for two motifs of 10 and 11 columns
        (-6, [[0, 0, 5, 10], [0, 6, 11, 11]]),
        (0, [[0, 10, 10], [0, 10, 11]]),
        (3, [[0, 3, 10, 10], [0, 0, 7, 11]]),
    ],
)
def test_display_coordinates_for_an_ungapped_path(offset, expected):
    start1, start2 = max(0, offset), max(0, -offset)
    overlap = min(10 - start1, 11 - start2)
    path = AlignmentPath(10, 11, ((start1, start1 + overlap), (start2, start2 + overlap)), False)
    assert _display_coordinates(path).tolist() == expected


def test_pearson_display_matches_the_offset_display(motif_agl3, motif_runx1):
    """format_alignment shows the same block #11 rendered from the offset."""
    result = align_motifs(motif_agl3, motif_runx1, method="pearson")
    motif2_used = motif_runx1.reverse_complement() if result.is_reverse_complement else motif_runx1
    sequences = [str(motif_agl3.consensus), str(motif2_used.consensus)]
    expected = str(Alignment(sequences, np.array([[0, 0, 5, 10], [0, 6, 11, 11]])))
    assert result.offset == -6
    assert format_alignment(motif_agl3, motif_runx1, result) == expected


def test_display_shows_the_internal_gap_and_both_full_motifs():
    query, target = one_hot("AAAACCGGGGT"), one_hot("AAAAGGGGA")
    result = align_motifs(query, target)
    # the stored core stops before the final column, whose similarity is 0
    assert result.alignment.coordinates == ((0, 4, 6, 10), (0, 4, 4, 8))
    # the display adds it back, unscored, and covers both motifs completely
    assert _display_coordinates(result.alignment).tolist() == [[0, 4, 6, 10, 11], [0, 4, 4, 8, 9]]
    text = format_alignment(query, target, result)
    assert "AAAACCGGGGT" in text
    assert "AAAA--GGGGA" in text


def test_display_uses_the_stored_path_without_recomputing():
    query, target = one_hot("ACGTACGT"), one_hot("ACGTACGT")
    shifted = AlignmentPath(8, 8, ((2, 6), (0, 4)), False)
    result = PearsonResult("X", "X", 0.0, shifted, {"min_overlap": 4, "both_strands": True})
    expected = str(
        Alignment(["ACGTACGT", "ACGTACGT"], np.array([[0, 2, 6, 8, 8], [0, 0, 4, 6, 8]]))
    )
    assert format_alignment(query, target, result) == expected
    assert result.alignment.coordinates == ((2, 6), (0, 4))


def test_display_rejects_motifs_that_do_not_match_the_result(motif_agl3, motif_runx1):
    result = align_motifs(motif_agl3, motif_runx1)
    with pytest.raises(ValueError, match="columns"):
        format_alignment(motif_runx1, motif_agl3, result)
