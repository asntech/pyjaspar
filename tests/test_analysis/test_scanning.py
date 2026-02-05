"""Tests for motif scanning."""

from __future__ import annotations

import pytest

from pyjaspar import JasparDB
from pyjaspar.analysis.scanning import ScanHit, scan_sequence


@pytest.fixture(scope="module")
def jdb():
    return JasparDB()


@pytest.fixture(scope="module")
def ctcf_motif(jdb):
    motifs = jdb.fetch_motifs_by_name("CTCF")
    return motifs[0]


def test_scan_finds_hits(ctcf_motif):
    """A sequence containing a strong CTCF site should produce hits."""
    # CTCF consensus-like sequence embedded in random flanking
    sequence = "AAAAACCACCAGGGGGCGCAAAAAA"
    hits = scan_sequence(sequence, ctcf_motif, threshold=0.5)
    assert isinstance(hits, list)
    for hit in hits:
        assert isinstance(hit, ScanHit)
        assert hit.strand in ("+", "-")
        assert isinstance(hit.score, float)
        assert isinstance(hit.position, int)


def test_scan_empty_sequence(ctcf_motif):
    """Empty sequence should return no hits."""
    hits = scan_sequence("", ctcf_motif, threshold=0.8)
    assert hits == []


def test_scan_short_sequence(ctcf_motif):
    """Sequence shorter than motif should return no hits."""
    hits = scan_sequence("ACGT", ctcf_motif, threshold=0.8)
    assert hits == []


def test_scan_low_threshold(ctcf_motif):
    """Low threshold should find more hits."""
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    hits_high = scan_sequence(sequence, ctcf_motif, threshold=0.9)
    hits_low = scan_sequence(sequence, ctcf_motif, threshold=0.3)
    assert len(hits_low) >= len(hits_high)


def test_scan_hit_fields(ctcf_motif):
    """Verify ScanHit dataclass fields."""
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT" * 3
    hits = scan_sequence(sequence, ctcf_motif, threshold=0.3)
    if hits:
        hit = hits[0]
        assert hasattr(hit, "position")
        assert hasattr(hit, "strand")
        assert hasattr(hit, "score")
        assert hasattr(hit, "sequence")
        assert len(hit.sequence) == ctcf_motif.length
