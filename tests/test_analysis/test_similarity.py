"""Tests for motif similarity measures."""

from __future__ import annotations

import pytest

from pyjaspar import JasparDB
from pyjaspar.analysis.similarity import (
    best_correlation,
    euclidean_distance,
    kl_divergence,
    pearson_correlation,
)


@pytest.fixture(scope="module")
def jdb():
    return JasparDB()


@pytest.fixture(scope="module")
def motif_agl3(jdb):
    return jdb.fetch_motif_by_id("MA0001.1")


@pytest.fixture(scope="module")
def motif_runx1(jdb):
    return jdb.fetch_motif_by_id("MA0002.1")


def test_pearson_self(motif_agl3):
    """Self-correlation should be close to 1.0."""
    score = pearson_correlation(motif_agl3, motif_agl3, offset=0)
    assert score > 0.99


def test_pearson_different(motif_agl3, motif_runx1):
    """Different motifs should have lower correlation."""
    score = pearson_correlation(motif_agl3, motif_runx1, offset=0)
    assert -1 <= score <= 1


def test_pearson_no_overlap(motif_agl3, motif_runx1):
    """Very large offset should give 0 (no overlap)."""
    score = pearson_correlation(motif_agl3, motif_runx1, offset=100)
    assert score == 0.0


def test_best_correlation(motif_agl3, motif_runx1):
    score, offset, is_rc = best_correlation(motif_agl3, motif_runx1, min_overlap=4)
    assert isinstance(score, float)
    assert isinstance(offset, int)
    assert isinstance(is_rc, bool)
    assert -1 <= score <= 1


def test_best_correlation_self(motif_agl3):
    score, offset, is_rc = best_correlation(motif_agl3, motif_agl3, min_overlap=4)
    assert score > 0.9
    assert offset == 0
    assert is_rc is False


def test_euclidean_self(motif_agl3):
    """Self-distance should be 0."""
    dist = euclidean_distance(motif_agl3, motif_agl3, offset=0)
    assert dist < 0.01


def test_euclidean_different(motif_agl3, motif_runx1):
    dist = euclidean_distance(motif_agl3, motif_runx1, offset=0)
    assert dist > 0


def test_euclidean_no_overlap(motif_agl3, motif_runx1):
    dist = euclidean_distance(motif_agl3, motif_runx1, offset=100)
    assert dist == float("inf")


def test_kl_self(motif_agl3):
    """Self KL-divergence should be very close to 0."""
    kl = kl_divergence(motif_agl3, motif_agl3, offset=0)
    assert kl < 0.01


def test_kl_nonnegative(motif_agl3, motif_runx1):
    kl = kl_divergence(motif_agl3, motif_runx1, offset=0)
    assert kl >= 0


def test_kl_no_overlap(motif_agl3, motif_runx1):
    kl = kl_divergence(motif_agl3, motif_runx1, offset=100)
    assert kl == float("inf")
