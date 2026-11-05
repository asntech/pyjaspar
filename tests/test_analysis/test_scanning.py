"""Tests for motif scanning."""

from __future__ import annotations

import math
import random

import pytest
from Bio.motifs.jaspar import Motif
from Bio.Seq import Seq

from pyjaspar import JasparDB
from pyjaspar.analysis.scanning import ScanHit, scan_sequence


@pytest.fixture
def acg_motif():
    return Motif(
        matrix_id="TEST",
        name="ACG",
        counts={base: [10 if letter == base else 0 for letter in "ACG"] for base in "ACGT"},
    )


@pytest.mark.parametrize("sequence", ["CGTAAAAAAA", "AAAACGTAAA", "AAAAAAACGT"])
def test_reverse_hit_uses_original_sequence_coordinates(acg_motif, sequence):
    hits = scan_sequence(sequence, acg_motif, threshold=0.99)
    reverse_hits = [hit for hit in hits if hit.strand == "-"]
    assert len(reverse_hits) == 1
    hit = reverse_hits[0]
    assert hit.position == sequence.index("CGT")
    assert hit.sequence == "CGT"
    assert hit.sequence == sequence[hit.position : hit.position + acg_motif.length]
    pssm = acg_motif.counts.normalize(0.001).log_odds()
    assert hit.score == pytest.approx(pssm.calculate(Seq(hit.sequence).reverse_complement()))


def test_mixed_strand_hits_are_sorted_in_original_coordinates(acg_motif):
    sequence = "AAAACGTAAA"
    hits = scan_sequence(sequence, acg_motif, threshold=0.99)
    assert [(hit.position, hit.strand, hit.sequence) for hit in hits] == [
        (3, "+", "ACG"),
        (4, "-", "CGT"),
    ]
    assert scan_sequence(sequence, acg_motif, threshold=0.99, both_strands=False) == hits[:1]


def test_reverse_only_site_is_excluded_when_both_strands_is_false(acg_motif):
    assert scan_sequence("CGTAAAAAAA", acg_motif, threshold=0.99, both_strands=False) == []


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


def _brute_force_hits(sequence, motif, threshold, pseudocount=0.001):
    """Score every window on both strands without Biopython's search.

    The log-odds follow what ``scan_sequence`` asks Biopython for: the pseudocount
    is added to each base's count (``counts.normalize``), the background is uniform,
    and scores are log2 ratios (``log_odds``). Only the scoring rule is shared;
    positions and strands come from plain slicing here, so a coordinate error in
    ``scan_sequence`` cannot be hidden by the same error in the reference.
    """
    complement = str.maketrans("ACGT", "TGCA")
    columns = []
    for i in range(motif.length):
        counts = {base: motif.counts[base][i] for base in "ACGT"}
        total = sum(counts.values()) + 4 * pseudocount  # added to each base
        columns.append(
            {base: math.log2((counts[base] + pseudocount) / total / 0.25) for base in "ACGT"}
        )
    best = sum(max(column.values()) for column in columns)
    worst = sum(min(column.values()) for column in columns)
    cutoff = worst + threshold * (best - worst)

    def score(window):
        return sum(column[base] for column, base in zip(columns, window, strict=True))

    hits = {}
    for start in range(len(sequence) - motif.length + 1):
        window = sequence[start : start + motif.length]
        for strand, read in (("+", window), ("-", window.translate(complement)[::-1])):
            value = score(read)
            if value >= cutoff:
                hits[(start, strand)] = value
    return hits


@pytest.mark.parametrize("seed", range(5))
def test_scan_matches_brute_force_scan_on_both_strands(ctcf_motif, seed):
    """Every hit, on both strands, is where an independent scan of all windows finds it.

    This checks pyjaspar's handling of Biopython's search results (in particular the
    negative positions it returns for the minus strand) against a direct calculation
    that does not use Biopython at all.
    """
    rng = random.Random(seed)
    for _ in range(20):
        sequence = "".join(
            rng.choice("ACGT") for _ in range(rng.randint(ctcf_motif.length + 1, 120))
        )
        threshold = rng.choice([0.5, 0.6, 0.7])
        hits = scan_sequence(sequence, ctcf_motif, threshold=threshold)
        expected = _brute_force_hits(sequence, ctcf_motif, threshold)
        assert {(hit.position, hit.strand) for hit in hits} == set(expected)
        for hit in hits:
            # Biopython scores in single precision
            assert hit.score == pytest.approx(expected[(hit.position, hit.strand)], abs=1e-4)
            assert hit.sequence == sequence[hit.position : hit.position + ctcf_motif.length]
