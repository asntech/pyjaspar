"""Tests for method selection, option validation, alignment paths, and search semantics."""

from __future__ import annotations

import math
import random

import numpy as np
import pytest
from Bio.motifs.jaspar import Motif

from pyjaspar import JasparDB
from pyjaspar.analysis import (
    MatrixAlignHit,
    PearsonHit,
    align_motifs,
    search_profiles,
)


@pytest.fixture(scope="module")
def jdb():
    return JasparDB("JASPAR2026")


@pytest.fixture(scope="module")
def ctcf(jdb):
    return jdb.fetch_motif_by_id("MA0139.2")


@pytest.fixture(scope="module")
def sample(jdb):
    """A spread of widths (6 to 33 columns), including palindromes and wide profiles."""
    ids = [
        "MA0139.2",
        "MA1930.2",
        "MA0162.2",
        "MA0079.5",
        "MA0004.1",
        "MA0854.2",
        "MA0114.4",
        "MA0107.1",
        "MA0149.1",
        "MA2100.1",
        "MA0001.1",
        "MA0002.1",
    ]
    return [jdb.fetch_motif_by_id(i) for i in ids]


def one_hot(sequence, matrix_id="X"):
    return Motif(
        matrix_id=matrix_id,
        name=matrix_id,
        counts={b: [10 if c == b else 0 for c in sequence] for b in "ACGT"},
    )


def with_flat_columns(motif, position, how_many):
    counts = {b: list(motif.counts[b]) for b in "ACGT"}
    total = sum(counts[b][position] for b in "ACGT")
    for b in "ACGT":
        counts[b][position:position] = [total / 4] * how_many
    return Motif(matrix_id="", name="flat", counts=counts)


def frequencies(motif):
    counts = np.array([motif.counts[b] for b in "ACGT"], dtype=float).T
    return counts / counts.sum(axis=1, keepdims=True)


def segments(path):
    query, target = path.coordinates
    return list(zip(query, query[1:], target, target[1:], strict=False))


def rescore(a, b, result):
    """Score recomputed from the stored path alone, independent of the search code."""
    fa = frequencies(a)
    fb = frequencies(b)
    if result.is_reverse_complement:
        fb = fb[::-1, ::-1]
    paired = [
        (q0 + k, t0 + k)
        for q0, q1, t0, t1 in segments(result.alignment)
        if q1 > q0 and t1 > t0
        for k in range(q1 - q0)
    ]
    if result.method == "pearson":
        correlations = []
        for i, j in paired:
            x, y = fa[i], fb[j]
            if np.std(x) == 0 or np.std(y) == 0:
                correlations.append(0.0)
            else:
                correlations.append(float(np.corrcoef(x, y)[0, 1]))
        return float(np.mean(correlations))
    score = sum(2.0 - float(((fa[i] - fb[j]) ** 2).sum()) for i, j in paired)
    if result.gaps:
        options = result.parameters
        score -= options["open_penalty"] + (result.gaps - 1) * options["ext_penalty"]
    return score


# --- method selection and options ---


@pytest.mark.parametrize("method", ["Pearson", "matrix-align", "", None, 3, object()])
def test_unknown_methods_are_rejected(ctcf, method):
    with pytest.raises(ValueError, match="method"):
        align_motifs(ctcf, ctcf, method=method)
    with pytest.raises(ValueError, match="method"):
        search_profiles(ctcf, [ctcf], method=method)


@pytest.mark.parametrize(
    ("method", "options"),
    [
        ("pearson", {"open_penalty": 3.0}),
        ("pearson", {"ext_penalty": 0.01}),
        ("matrix_align", {"min_overlap": 4}),
    ],
)
def test_options_of_the_other_method_are_rejected(ctcf, method, options):
    with pytest.raises(ValueError, match="only to method="):
        align_motifs(ctcf, ctcf, method=method, **options)
    with pytest.raises(ValueError, match="only to method="):
        search_profiles(ctcf, [], method=method, **options)


@pytest.mark.parametrize("value", [0, -1, True, 2.5, "4"])
def test_invalid_min_overlap_is_rejected(ctcf, value):
    with pytest.raises(ValueError, match="min_overlap"):
        align_motifs(ctcf, ctcf, method="pearson", min_overlap=value)


@pytest.mark.parametrize("value", [-1.0, math.nan, math.inf, True, "3"])
def test_invalid_penalties_are_rejected(ctcf, value):
    with pytest.raises(ValueError, match="open_penalty"):
        align_motifs(ctcf, ctcf, open_penalty=value)
    with pytest.raises(ValueError, match="ext_penalty"):
        align_motifs(ctcf, ctcf, ext_penalty=value)


@pytest.mark.parametrize("value", [1, None, "yes"])
def test_both_strands_must_be_a_bool(ctcf, value):
    with pytest.raises(ValueError, match="both_strands"):
        align_motifs(ctcf, ctcf, both_strands=value)


def test_parameters_record_the_effective_options(ctcf):
    result = align_motifs(ctcf, ctcf, open_penalty=2, both_strands=False)
    assert result.parameters == {"open_penalty": 2.0, "ext_penalty": 0.01, "both_strands": False}
    assert align_motifs(ctcf, ctcf, method="pearson", min_overlap=6).parameters == {
        "min_overlap": 6,
        "both_strands": True,
    }


# --- invalid profiles ---


def test_profile_with_a_zero_column_is_rejected_by_both_methods(ctcf):
    broken = one_hot("ACGT")
    for b in "ACGT":
        broken.counts[b][2] = 0
    for method in ("matrix_align", "pearson"):
        with pytest.raises(ValueError, match="no counts"):
            align_motifs(ctcf, broken, method=method)


def test_negative_counts_are_rejected(ctcf):
    broken = Motif(matrix_id="NEG", name="neg", counts={b: [5, 5, 5, 5] for b in "ACGT"})
    broken.counts["A"][1] = -1
    with pytest.raises(ValueError, match="negative"):
        align_motifs(ctcf, broken)


def test_pearson_rejects_a_motif_shorter_than_min_overlap(ctcf):
    short = one_hot("ACG", "SHORT")
    with pytest.raises(ValueError, match="fewer than min_overlap=4"):
        align_motifs(ctcf, short, method="pearson")
    # Matrix Align has no minimum overlap
    assert align_motifs(ctcf, short).overlap >= 1


def test_pearson_search_rejects_a_candidate_shorter_than_min_overlap(ctcf):
    with pytest.raises(ValueError, match="SHORT"):
        search_profiles(ctcf, [ctcf, one_hot("ACG", "SHORT")], method="pearson")


# --- alignment paths ---


@pytest.mark.parametrize("method", ["matrix_align", "pearson"])
def test_paths_are_well_formed_and_rescore_to_the_reported_score(sample, method):
    for a in sample:
        for b in sample:
            result = align_motifs(a, b, method=method)
            path = result.alignment
            query, target = path.coordinates
            assert (path.query_length, path.target_length) == (a.length, b.length)
            assert len(query) == len(target) >= 2
            assert all(0 <= x <= a.length for x in query) and all(
                0 <= y <= b.length for y in target
            )
            parts = segments(path)
            assert all(q1 >= q0 and t1 >= t0 and (q1 > q0 or t1 > t0) for q0, q1, t0, t1 in parts)
            # paired at both ends, at most one gap run (none for Pearson)
            for q0, q1, t0, t1 in (parts[0], parts[-1]):
                assert q1 - q0 == t1 - t0 > 0
            gap_runs = [p for p in parts if (p[1] > p[0]) != (p[3] > p[2])]
            assert len(gap_runs) <= (0 if method == "pearson" else 1)
            assert result.alignment_length == result.overlap + result.gaps
            assert result.offset == query[0] - target[0]
            assert rescore(a, b, result) == pytest.approx(result.score, abs=1e-9)


def test_a_query_insertion_is_a_gap_in_the_target_row(ctcf):
    query = with_flat_columns(ctcf, position=7, how_many=2)
    result = align_motifs(query, ctcf)
    assert result.alignment.coordinates == ((0, 7, 9, 17), (0, 7, 7, 15))
    assert (result.gaps, result.overlap, result.alignment_length, result.offset) == (2, 15, 17, 0)


def test_a_target_insertion_is_a_gap_in_the_query_row(ctcf):
    target = with_flat_columns(ctcf, position=7, how_many=2)
    result = align_motifs(ctcf, target)
    assert result.alignment.coordinates == ((0, 7, 7, 15), (0, 7, 9, 17))


def test_reverse_complement_coordinates_refer_to_the_oriented_target(ctcf):
    result = align_motifs(ctcf, ctcf.reverse_complement())
    assert result.is_reverse_complement is True
    assert result.alignment.coordinates == ((0, 15), (0, 15))
    assert result.score == pytest.approx(30.0)
    assert result.motif2_id is None  # a reverse-complemented Biopython motif has no matrix ID


def test_trailing_zero_similarity_columns(ctcf):
    # ungapped: a final column with similarity 0 is still part of the alignment
    ungapped = align_motifs(one_hot("ACGA"), one_hot("ACGC"), both_strands=False)
    assert ungapped.alignment.coordinates == ((0, 4), (0, 4))
    assert ungapped.score == pytest.approx(6.0)
    # gapped: after the gap it is not
    gapped = align_motifs(one_hot("AAAACCGGGGT"), one_hot("AAAAGGGGA"))
    assert gapped.alignment.coordinates == ((0, 4, 6, 10), (0, 4, 4, 8))
    assert gapped.score == pytest.approx(16.0 - 3.01)


def test_forward_only(ctcf):
    result = align_motifs(ctcf, ctcf.reverse_complement(), both_strands=False)
    assert result.is_reverse_complement is False
    assert result.score < 30.0


# --- search ---


def test_search_default_is_matrix_align(ctcf, sample):
    hits = search_profiles(ctcf, sample)
    assert all(isinstance(h, MatrixAlignHit) for h in hits)
    assert hits[0].matrix_id == "MA0139.2"
    assert hits[0].score == pytest.approx(30.0)


def test_percent_score_uses_the_narrowest_profile_seen_in_matrix_id_order(ctcf, sample):
    hits = {h.matrix_id: h for h in search_profiles(ctcf, sample)}
    narrowest = ctcf.length
    for motif in sorted(sample, key=lambda m: tuple(int(x) for x in m.matrix_id[2:].split("."))):
        narrowest = min(narrowest, motif.length)
        hit = hits[motif.matrix_id]
        assert hit.percent_score == pytest.approx(100 * hit.score / (2 * narrowest))
    # in this set the 6-column MA0004.1 comes before CTCF, so CTCF's own row is 30 / 12
    assert hits["MA0139.2"].percent_score == pytest.approx(250.0)


def test_pearson_search_returns_pearson_hits_best_first(ctcf, sample):
    hits = search_profiles(ctcf, sample, method="pearson")
    assert all(isinstance(h, PearsonHit) for h in hits)
    correlations = [h.correlation for h in hits]
    assert correlations == sorted(correlations, reverse=True)
    assert not hasattr(hits[0], "percent_score")
    assert search_profiles(ctcf, sample, method="pearson", sort_by="score") == hits


@pytest.mark.parametrize("method", ["matrix_align", "pearson"])
def test_every_hit_equals_the_direct_pair_comparison(ctcf, sample, method):
    for hit in search_profiles(ctcf, sample, method=method, both_strands=False):
        candidate = next(c for c in sample if c.matrix_id == hit.matrix_id)
        assert hit.result == align_motifs(ctcf, candidate, method=method, both_strands=False)
        assert hit.width == candidate.length


@pytest.mark.parametrize(
    ("method", "sort_by"),
    [("pearson", "percent_score"), ("pearson", "pvalue"), ("matrix_align", "correlation")],
)
def test_inapplicable_sort_fields_are_rejected_before_scoring(ctcf, method, sort_by):
    with pytest.raises(ValueError, match="sort_by"):
        search_profiles(ctcf, [], method=method, sort_by=sort_by)


@pytest.mark.parametrize("top", [-1, True, 2.5, "3"])
def test_invalid_top_is_rejected(ctcf, top):
    with pytest.raises(ValueError, match="top"):
        search_profiles(ctcf, [], top=top)


def test_top_zero_and_truncation_keep_percent_scores(ctcf, sample):
    assert search_profiles(ctcf, sample, top=0) == []
    full = search_profiles(ctcf, sample)
    assert search_profiles(ctcf, sample, top=3) == full[:3]


def test_candidates_can_be_a_generator(ctcf, sample):
    assert search_profiles(ctcf, (c for c in sample)) == search_profiles(ctcf, sample)


def test_duplicates_and_missing_ids_are_kept(ctcf):
    no_id = ctcf.reverse_complement()  # no matrix_id
    hits = search_profiles(ctcf, [ctcf, ctcf, no_id])
    assert len(hits) == 3
    assert [h.matrix_id for h in hits].count("MA0139.2") == 2
    assert any(h.matrix_id is None for h in hits)


def test_percent_score_does_not_depend_on_input_order(ctcf, sample):
    shuffled = list(sample)
    random.Random(1).shuffle(shuffled)
    expected = {h.matrix_id: h.percent_score for h in search_profiles(ctcf, sample)}
    actual = {h.matrix_id: h.percent_score for h in search_profiles(ctcf, shuffled)}
    assert actual == pytest.approx(expected)
