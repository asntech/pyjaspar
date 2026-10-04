"""Tests for the memesuite-lite wrapper.

The tests that call memesuite-lite are skipped when it is not installed
(``pip install pyjaspar[meme]``); the first of them takes about half a minute
because memesuite-lite compiles its code on first use.
"""

from __future__ import annotations

import sys

import pytest
from Bio.motifs.jaspar import Motif

from pyjaspar import JasparDB
from pyjaspar.analysis import best_correlation
from pyjaspar.analysis.meme import TomtomHit, tomtom


def test_missing_memesuite_lite_says_how_to_install(monkeypatch, jdb):
    monkeypatch.setitem(sys.modules, "memelite", None)
    query = jdb.fetch_motif_by_id("MA0139.2")
    with pytest.raises(ImportError, match=r"pip install pyjaspar\[meme\]"):
        tomtom(query, [query])


@pytest.fixture(scope="module")
def jdb():
    return JasparDB("JASPAR2026")


@pytest.fixture(scope="module")
def ctcf(jdb):
    return jdb.fetch_motif_by_id("MA0139.2")


@pytest.fixture(scope="module")
def others(jdb, ctcf):
    """Enough other profiles for memesuite-lite to give stable p-values (it warns below 25)."""
    candidates = jdb.fetch_motifs(
        collection=["CORE"], tax_group=["Vertebrates"], all_versions=False
    )
    return [c for c in candidates if c.matrix_id != ctcf.matrix_id][:60]


def renamed(motif, matrix_id, counts=None):
    counts = counts or {b: list(motif.counts[b]) for b in "ACGT"}
    return Motif(matrix_id=matrix_id, name=matrix_id, counts=counts)


@pytest.fixture(scope="module")
def memelite():
    return pytest.importorskip("memelite")


def test_query_is_its_own_best_hit(memelite, ctcf, others):
    hits = tomtom(ctcf, [*others, ctcf])
    best = hits[0]
    assert isinstance(best, TomtomHit)
    assert best.matrix_id == "MA0139.2"
    assert best.name == "CTCF"
    assert best.offset == 0
    assert best.overlap == ctcf.length
    assert best.is_reverse_complement is False
    assert 0 <= best.pvalue < 1e-10


def test_reverse_complement_candidate_is_found_in_reverse(memelite, ctcf, others):
    twin = renamed(ctcf, "RC", counts=ctcf.reverse_complement().counts)
    hit = next(h for h in tomtom(ctcf, [*others, twin]) if h.matrix_id == "RC")
    assert hit.is_reverse_complement is True
    assert hit.offset == 0
    assert hit.overlap == ctcf.length


def test_offset_is_the_query_position_minus_the_candidate_position(memelite, ctcf, others):
    # four flat columns in front of the query: its first column meets the candidate's fifth
    counts = {b: [25.0] * 4 + list(ctcf.counts[b]) for b in "ACGT"}
    shifted = renamed(ctcf, "SHIFTED", counts=counts)
    hit = next(h for h in tomtom(ctcf, [*others, shifted]) if h.matrix_id == "SHIFTED")
    assert hit.offset == -4
    assert hit.overlap == ctcf.length
    # the same sign as the offset of best_correlation
    assert best_correlation(ctcf, shifted)[1] == hit.offset


def test_hits_are_sorted_by_pvalue_and_top_limits_them(memelite, ctcf, others):
    hits = tomtom(ctcf, others)
    assert len(hits) == len(others)
    assert [h.pvalue for h in hits] == sorted(h.pvalue for h in hits)
    assert tomtom(ctcf, others, top=3) == hits[:3]


def test_no_candidates_gives_no_hits(memelite, ctcf):
    assert tomtom(ctcf, []) == []


def test_column_without_counts_is_rejected(memelite, ctcf, others):
    counts = {b: list(ctcf.counts[b]) for b in "ACGT"}
    for b in "ACGT":
        counts[b][3] = 0
    with pytest.raises(ValueError, match="no counts"):
        tomtom(renamed(ctcf, "EMPTY", counts=counts), others)
