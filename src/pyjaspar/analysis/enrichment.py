"""Motif enrichment analysis.

Tests whether motifs are enriched in a foreground set of sequences
compared to a background set using Fisher's exact test.
"""

from __future__ import annotations

from dataclasses import dataclass

try:
    from scipy import stats as scipy_stats
except ImportError:
    raise ImportError(
        "The pyjaspar.analysis.enrichment module requires scipy. "
        "Install with: pip install pyjaspar[analysis]"
    ) from None

from Bio.motifs.jaspar import Motif
from Bio.Seq import Seq

from .scanning import scan_sequence


@dataclass
class EnrichmentResult:
    """Result of enrichment analysis for a single motif.

    Attributes:
        motif_id: JASPAR matrix ID.
        motif_name: Transcription factor name.
        fg_hits: Number of foreground sequences with at least one hit.
        bg_hits: Number of background sequences with at least one hit.
        fg_total: Total foreground sequences.
        bg_total: Total background sequences.
        fold_enrichment: Ratio of foreground hit rate to background hit rate.
        pvalue: Fisher's exact test p-value.
        qvalue: Benjamini-Hochberg adjusted p-value (if computed).
    """

    motif_id: str
    motif_name: str
    fg_hits: int
    bg_hits: int
    fg_total: int
    bg_total: int
    fold_enrichment: float
    pvalue: float
    qvalue: float | None = None


def _count_sequences_with_hits(
    sequences: list[str | Seq],
    motif: Motif,
    threshold: float,
) -> int:
    """Count how many sequences have at least one motif hit."""
    count = 0
    for seq in sequences:
        hits = scan_sequence(seq, motif, threshold=threshold, both_strands=True)
        if len(hits) > 0:
            count += 1
    return count


def _bh_correction(pvalues: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg correction to a list of p-values."""
    n = len(pvalues)
    if n == 0:
        return []

    # Sort indices by p-value
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    qvalues = [0.0] * n

    cummin = 1.0
    for rank_minus_1 in range(n - 1, -1, -1):
        idx, pval = indexed[rank_minus_1]
        rank = rank_minus_1 + 1
        adjusted = pval * n / rank
        cummin = min(cummin, adjusted)
        qvalues[idx] = min(cummin, 1.0)

    return qvalues


def motif_enrichment(
    foreground: list[str | Seq],
    background: list[str | Seq],
    motifs: list[Motif],
    threshold: float = 0.8,
) -> list[EnrichmentResult]:
    """Test motif enrichment in foreground vs background sequences.

    For each motif, counts the number of sequences with at least one
    PSSM hit above the threshold, then applies a one-sided Fisher's
    exact test (greater). Results are sorted by p-value with
    Benjamini-Hochberg correction applied.

    Args:
        foreground: Target sequences to test for enrichment.
        background: Control sequences.
        motifs: JASPAR motifs to test.
        threshold: PSSM score threshold as fraction of max score (0.0-1.0).

    Returns:
        List of EnrichmentResult, sorted by p-value (ascending).
    """
    fg_total = len(foreground)
    bg_total = len(background)
    results: list[EnrichmentResult] = []

    for motif in motifs:
        fg_hits = _count_sequences_with_hits(foreground, motif, threshold)
        bg_hits = _count_sequences_with_hits(background, motif, threshold)

        # Fisher's exact test (2x2 contingency table)
        # [[fg_hits, fg_no_hit], [bg_hits, bg_no_hit]]
        fg_no_hit = fg_total - fg_hits
        bg_no_hit = bg_total - bg_hits

        table = [[fg_hits, fg_no_hit], [bg_hits, bg_no_hit]]
        _, pvalue = scipy_stats.fisher_exact(table, alternative="greater")

        # Fold enrichment
        fg_rate = fg_hits / fg_total if fg_total > 0 else 0
        bg_rate = bg_hits / bg_total if bg_total > 0 else 0
        fold = fg_rate / bg_rate if bg_rate > 0 else float("inf") if fg_rate > 0 else 1.0

        results.append(
            EnrichmentResult(
                motif_id=motif.matrix_id,
                motif_name=motif.name,
                fg_hits=fg_hits,
                bg_hits=bg_hits,
                fg_total=fg_total,
                bg_total=bg_total,
                fold_enrichment=fold,
                pvalue=pvalue,
            )
        )

    # Apply BH correction
    pvalues = [r.pvalue for r in results]
    qvalues = _bh_correction(pvalues)
    for r, q in zip(results, qvalues, strict=True):
        r.qvalue = q

    results.sort(key=lambda r: r.pvalue)
    return results
