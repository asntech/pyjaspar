"""Tests for searching profiles with a query motif.

Expected scores come from two public sources, recorded next to each value:
the TFBSTools manual example, and the "Score" column of the JASPAR web tool
(https://jaspar.elixir.no/align/, queried on 2026-09-30 with collection CORE,
taxonomic group vertebrates, latest versions).
"""

from __future__ import annotations

import random

import pytest
from Bio.motifs.jaspar import Motif

from pyjaspar import JasparDB
from pyjaspar.analysis import ProfileHit, align_motifs, search_profiles


@pytest.fixture(scope="module")
def jdb():
    return JasparDB("JASPAR2026")


@pytest.fixture(scope="module")
def egr1(jdb):
    return jdb.fetch_motif_by_id("MA0162.2")


@pytest.fixture(scope="module")
def ctcf(jdb):
    return jdb.fetch_motif_by_id("MA0139.2")


@pytest.fixture(scope="module")
def core_vertebrates(jdb):
    """The profiles the web form searched: CORE, vertebrates, latest versions."""
    return jdb.fetch_motifs(collection=["CORE"], tax_group=["Vertebrates"], all_versions=False)


def with_flat_columns(motif, position, how_many):
    """Copy of ``motif`` with uninformative columns inserted at ``position``."""
    counts = {b: list(motif.counts[b]) for b in "ACGT"}
    total = sum(counts[b][position] for b in "ACGT")
    for b in "ACGT":
        counts[b][position:position] = [total / 4] * how_many
    return Motif(matrix_id="", name="flat", counts=counts)


def test_tfbstools_manual_example():
    # TFBSTools manual, PFMSimilarity(MA0003.2, MA0004.1) on JASPAR2014 data: score 7.294736
    db = JasparDB("JASPAR2014")
    a = db.fetch_motif_by_id("MA0003.2")
    b = db.fetch_motif_by_id("MA0004.1")
    assert align_motifs(a, b).score == pytest.approx(7.294736, abs=1e-4)
    assert align_motifs(b, a).score == pytest.approx(7.294736, abs=1e-4)


@pytest.mark.parametrize(
    ("matrix_id", "web_score"),
    [
        # JASPAR web tool, query = human EGR1 MA0162.2 (JASPAR2026 counts)
        ("MA0002.3", 13.1023),
        ("MA0003.5", 13.3877),
        ("MA0004.1", 8.64621),
        ("MA0006.2", 8.67515),
        ("MA0007.4", 17.2406),
        ("MA1723.2", 24.6726),
        ("MA0073.2", 23.9691),
        ("MA1929.2", 23.8812),
        ("MA0162.5", 18.9821),
        # the best alignments of these two contain a gap
        ("MA2457.1", 16.9322),
        ("MA1654.2", 20.3021),
    ],
)
def test_matches_web_tool_for_egr1_query(jdb, egr1, matrix_id, web_score):
    candidate = jdb.fetch_motif_by_id(matrix_id)
    assert align_motifs(egr1, candidate).score == pytest.approx(web_score, abs=1e-3)


@pytest.mark.parametrize(
    ("query_id", "matrix_id", "is_reverse_complement", "gaps", "offset", "alignment_length"),
    [
        # Output of the matrix_aligner program (jaspar_tools, commit 54aa156, built with g++)
        # on JASPAR2026 counts, with the default penalties; checked on 2026-10-04.
        ("MA0162.2", "MA0002.3", True, 0, 0, 9),
        ("MA0162.2", "MA0004.1", True, 0, 8, 6),
        ("MA0162.2", "MA2457.1", True, 8, 0, 22),
        ("MA0162.2", "MA1654.2", True, 6, 0, 20),
        ("MA0162.2", "MA2331.1", False, 4, -1, 18),
        ("MA0162.2", "MA1723.2", True, 0, -3, 14),
        ("MA0139.2", "MA0139.2", False, 0, 0, 15),
        ("MA0139.2", "MA1930.2", False, 0, -18, 15),
        ("MA0139.2", "MA0107.1", True, 0, 4, 10),
        ("MA0139.2", "MA0149.1", False, 2, 0, 17),
        ("MA0139.2", "MA2100.1", False, 1, -4, 15),
        # MA0854.2 is its own reverse complement: both orientations score the same
        ("MA0139.2", "MA0854.2", True, 0, 0, 8),
    ],
)
def test_matches_matrix_aligner_alignment(
    jdb, query_id, matrix_id, is_reverse_complement, gaps, offset, alignment_length
):
    result = align_motifs(jdb.fetch_motif_by_id(query_id), jdb.fetch_motif_by_id(matrix_id))
    assert result.is_reverse_complement is is_reverse_complement
    assert result.gaps == gaps
    assert result.offset == offset
    assert result.alignment_length == alignment_length


def test_inserted_columns_match_matrix_aligner_alignment(jdb, ctcf):
    # matrix_aligner, query = CTCF MA0139.2 with two flat columns inserted after column 7
    query = with_flat_columns(ctcf, position=7, how_many=2)
    own_row = align_motifs(query, jdb.fetch_motif_by_id("MA0139.2"))
    assert (own_row.gaps, own_row.offset, own_row.alignment_length) == (2, 0, 17)
    other_row = align_motifs(query, jdb.fetch_motif_by_id("MA1930.2"))
    assert (other_row.gaps, other_row.offset, other_row.alignment_length) == (9, -7, 26)


@pytest.mark.parametrize(
    ("matrix_id", "web_score"),
    [
        # JASPAR web tool, query = CTCF MA0139.2: its own row is 2 * 15 columns
        ("MA0139.2", 30.0),
        ("MA1930.2", 29.888),
        ("MA1929.2", 27.8202),
        ("MA2510.1", 24.2936),
    ],
)
def test_matches_web_tool_for_ctcf_query(jdb, ctcf, matrix_id, web_score):
    candidate = jdb.fetch_motif_by_id(matrix_id)
    assert align_motifs(ctcf, candidate).score == pytest.approx(web_score, abs=1e-3)


@pytest.mark.parametrize(
    ("matrix_id", "web_score"),
    [
        # JASPAR web tool, query = CTCF MA0139.2 with two flat columns inserted after column 7;
        # its own row is 30 - (3 + 0.01): one gap of two columns
        ("MA0139.2", 26.99),
        ("MA1930.2", 28.9098),
        ("MA1929.2", 27.8626),
        ("MA1987.2", 27.3149),
    ],
)
def test_matches_web_tool_for_query_with_inserted_columns(jdb, ctcf, matrix_id, web_score):
    query = with_flat_columns(ctcf, position=7, how_many=2)
    candidate = jdb.fetch_motif_by_id(matrix_id)
    assert align_motifs(query, candidate).score == pytest.approx(web_score, abs=1e-3)


def test_self_alignment_scores_two_per_column(ctcf):
    result = align_motifs(ctcf, ctcf)
    assert result.score == pytest.approx(2 * ctcf.length)
    assert result.gaps == 0
    assert result.is_reverse_complement is False


def test_reverse_complement_does_not_change_the_score(jdb, egr1):
    candidate = jdb.fetch_motif_by_id("MA0002.3")
    forward = align_motifs(egr1, candidate)
    reverse = align_motifs(egr1, candidate.reverse_complement())
    assert forward.score == pytest.approx(reverse.score)
    assert forward.is_reverse_complement != reverse.is_reverse_complement


def test_inserted_columns_cost_the_gap_penalty(ctcf):
    one = align_motifs(with_flat_columns(ctcf, 7, 1), ctcf)
    two = align_motifs(with_flat_columns(ctcf, 7, 2), ctcf)
    assert (one.gaps, two.gaps) == (1, 2)
    assert one.score == pytest.approx(30.0 - 3.0)
    assert two.score == pytest.approx(30.0 - 3.0 - 0.01)


def test_gap_penalties_are_parameters(ctcf):
    query = with_flat_columns(ctcf, 7, 2)
    cheap = align_motifs(query, ctcf, open_penalty=1.0, ext_penalty=0.5)
    assert cheap.score == pytest.approx(30.0 - 1.0 - 0.5)


def test_column_without_counts_is_rejected(ctcf):
    counts = {b: list(ctcf.counts[b]) for b in "ACGT"}
    for b in "ACGT":
        counts[b][3] = 0
    broken = Motif(matrix_id="broken", name="broken", counts=counts)
    with pytest.raises(ValueError, match="no counts"):
        align_motifs(ctcf, broken)


def test_search_returns_hits_best_first(jdb, ctcf):
    candidates = [jdb.fetch_motif_by_id(i) for i in ("MA1929.2", "MA0139.2", "MA1930.2")]
    hits = search_profiles(ctcf, candidates)
    assert [h.matrix_id for h in hits] == ["MA0139.2", "MA1930.2", "MA1929.2"]
    assert all(isinstance(h, ProfileHit) for h in hits)
    assert hits[0].name == "CTCF"
    assert hits[0].score == pytest.approx(30.0)
    assert hits[0].width == 15


@pytest.mark.parametrize(
    ("matrix_id", "web_percent"),
    [
        # JASPAR web tool "Percent Score", query = human EGR1 MA0162.2
        ("MA0002.3", 72.79055555555556),  # first row of the web table: m = 9
        ("MA0003.5", 74.37611111111111),
        ("MA0004.1", 72.05175),  # m = 6
        ("MA0006.2", 86.7515),  # m = 5
        ("MA0007.4", 172.406),  # above 100: m is still 5
        ("MA1723.2", 246.72599999999997),
        ("MA2557.1", 72.544375),  # width 4: m = 4 from here on
        ("MA2558.1", 144.775),
        ("MA2587.1", 271.37375000000003),
        ("MA0162.5", 189.82099999999997),
    ],
)
def test_percent_score_matches_web_tool_for_egr1_query(
    core_vertebrates, egr1, matrix_id, web_percent
):
    hits = {h.matrix_id: h for h in search_profiles(egr1, core_vertebrates)}
    assert hits[matrix_id].percent_score == pytest.approx(web_percent, abs=1e-2)


@pytest.mark.parametrize(
    ("matrix_id", "web_percent"),
    [
        # JASPAR web tool "Percent Score", query = CTCF MA0139.2
        ("MA0139.2", 300.0),
        ("MA1930.2", 298.88),
        ("MA0002.3", 70.58),
        ("MA2557.1", 69.513125),
    ],
)
def test_percent_score_matches_web_tool_for_ctcf_query(
    core_vertebrates, ctcf, matrix_id, web_percent
):
    hits = {h.matrix_id: h for h in search_profiles(ctcf, core_vertebrates)}
    assert hits[matrix_id].percent_score == pytest.approx(web_percent, abs=1e-2)


@pytest.mark.parametrize(
    ("matrix_id", "web_percent"),
    [
        # JASPAR web tool "Percent Score", query = CTCF MA0139.2 with two flat columns inserted
        ("MA0139.2", 269.9),
        ("MA1930.2", 289.098),
        ("MA0002.3", 75.495),
        ("MA2557.1", 76.967625),
    ],
)
def test_percent_score_matches_web_tool_for_inserted_columns_query(
    core_vertebrates, ctcf, matrix_id, web_percent
):
    query = with_flat_columns(ctcf, position=7, how_many=2)
    hits = {h.matrix_id: h for h in search_profiles(query, core_vertebrates)}
    assert hits[matrix_id].percent_score == pytest.approx(web_percent, abs=1e-2)


def test_percent_score_does_not_depend_on_the_order_of_the_candidates(core_vertebrates, egr1):
    shuffled = list(core_vertebrates)
    random.Random(0).shuffle(shuffled)
    expected = {h.matrix_id: h.percent_score for h in search_profiles(egr1, core_vertebrates)}
    actual = {h.matrix_id: h.percent_score for h in search_profiles(egr1, shuffled)}
    assert actual == pytest.approx(expected)


def test_percent_score_follows_the_given_order_without_standard_ids(ctcf):
    # matrix_id "" is not a JASPAR ID, so the candidates are taken in the order given
    wide = with_flat_columns(ctcf, 7, 2)  # 17 columns
    narrow = Motif(
        matrix_id="",
        name="narrow",
        counts={b: list(ctcf.counts[b])[:10] for b in "ACGT"},
    )
    hits = {h.name: h for h in search_profiles(ctcf, [wide, narrow])}
    # m = min(15, 17) = 15 for the first candidate, then min(15, 10) = 10 for the second
    assert hits["flat"].percent_score == pytest.approx(100 * hits["flat"].score / (2 * 15))
    assert hits["narrow"].percent_score == pytest.approx(100 * hits["narrow"].score / (2 * 10))


def test_search_sort_by_percent_score(core_vertebrates, egr1):
    by_score = search_profiles(egr1, core_vertebrates)
    by_percent = search_profiles(egr1, core_vertebrates, sort_by="percent_score")
    assert [h.score for h in by_score] == sorted((h.score for h in by_score), reverse=True)
    assert [h.percent_score for h in by_percent] == sorted(
        (h.percent_score for h in by_percent), reverse=True
    )


def test_search_top_truncates(jdb, ctcf):
    candidates = [jdb.fetch_motif_by_id(i) for i in ("MA1929.2", "MA0139.2", "MA1930.2")]
    assert len(search_profiles(ctcf, candidates, top=2)) == 2


def test_search_with_no_candidates(ctcf):
    assert search_profiles(ctcf, []) == []


def test_search_rejects_unknown_sort_key(ctcf):
    with pytest.raises(ValueError, match="sort_by"):
        search_profiles(ctcf, [], sort_by="pvalue")
