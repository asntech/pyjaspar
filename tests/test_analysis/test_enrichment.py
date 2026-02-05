"""Tests for motif enrichment analysis."""

from __future__ import annotations

import pytest

from pyjaspar import JasparDB
from pyjaspar.analysis.enrichment import EnrichmentResult, motif_enrichment


@pytest.fixture(scope="module")
def jdb():
    return JasparDB()


@pytest.fixture(scope="module")
def test_motifs(jdb):
    """A small set of motifs for enrichment testing."""
    m1 = jdb.fetch_motif_by_id("MA0001.1")
    m2 = jdb.fetch_motif_by_id("MA0002.1")
    return [m1, m2]


def test_enrichment_returns_results(test_motifs):
    """Basic enrichment test with synthetic data."""
    foreground = [
        "ACGTACGTACGTACGT",
        "TGCATGCATGCATGCA",
        "AAAACCCCGGGGTTTT",
    ]
    background = [
        "AAAAAAAAAAAAAAAA",
        "TTTTTTTTTTTTTTTT",
        "CCCCCCCCCCCCCCCC",
    ]
    results = motif_enrichment(foreground, background, test_motifs, threshold=0.3)
    assert isinstance(results, list)
    assert len(results) == 2  # one per motif
    for r in results:
        assert isinstance(r, EnrichmentResult)
        assert 0 <= r.pvalue <= 1
        assert r.fg_total == 3
        assert r.bg_total == 3
        assert r.qvalue is not None


def test_enrichment_sorted_by_pvalue(test_motifs):
    """Results should be sorted by p-value ascending."""
    foreground = ["ACGTACGTACGTACGT"] * 5
    background = ["AAAAAAAAAAAAAAAA"] * 5
    results = motif_enrichment(foreground, background, test_motifs, threshold=0.3)
    for i in range(len(results) - 1):
        assert results[i].pvalue <= results[i + 1].pvalue


def test_enrichment_empty_foreground(test_motifs):
    """Empty foreground should not crash."""
    results = motif_enrichment([], ["ACGTACGTACGTACGTACGT"], test_motifs, threshold=0.3)
    assert len(results) == 2
    for r in results:
        assert r.fg_total == 0
        assert r.fg_hits == 0


def test_enrichment_empty_background(test_motifs):
    """Empty background should not crash."""
    results = motif_enrichment(["ACGTACGTACGTACGTACGT"], [], test_motifs, threshold=0.3)
    assert len(results) == 2
    for r in results:
        assert r.bg_total == 0


def test_enrichment_result_fields(test_motifs):
    """Check all fields are populated."""
    foreground = ["ACGTACGTACGTACGT"]
    background = ["TTTTTTTTTTTTTTTT"]
    results = motif_enrichment(foreground, background, test_motifs, threshold=0.3)
    r = results[0]
    assert isinstance(r.motif_id, str)
    assert isinstance(r.motif_name, str)
    assert isinstance(r.fold_enrichment, float)
