"""Motif scanning against DNA sequences.

Wraps BioPython's PSSM search with a cleaner API and structured results.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio.motifs.jaspar import Motif
from Bio.Seq import Seq


@dataclass
class ScanHit:
    """A single motif hit in a sequence.

    Attributes:
        position: 0-based start position in the sequence.
        strand: ``"+"`` or ``"-"``.
        score: PSSM score at this position.
        sequence: The matched subsequence.
    """

    position: int
    strand: str
    score: float
    sequence: str


def scan_sequence(
    sequence: str | Seq,
    motif: Motif,
    threshold: float = 0.8,
    both_strands: bool = True,
    pseudocount: float = 0.001,
    background: dict[str, float] | None = None,
) -> list[ScanHit]:
    """Scan a DNA sequence for motif occurrences using the PSSM.

    The threshold is specified as a fraction of the maximum possible
    PSSM score.

    Args:
        sequence: DNA sequence to scan (string or Bio.Seq.Seq).
        motif: JASPAR motif to scan with.
        threshold: Score threshold as fraction of max score (0.0-1.0).
        both_strands: If True, also scan the reverse complement strand.
        pseudocount: Pseudocount for PSSM calculation.
        background: Background nucleotide frequencies.
            Defaults to uniform (0.25 each).

    Returns:
        List of ScanHit objects, sorted by position.
    """
    if background is None:
        background = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

    seq_str = str(sequence).upper()
    motif_length = motif.length

    # Sequence must be at least as long as the motif
    if len(seq_str) < motif_length:
        return []

    # Build PSSM via BioPython
    pwm = motif.counts.normalize(pseudocount)
    pssm = pwm.log_odds(background)

    max_score = pssm.max
    min_score = pssm.min
    score_threshold = min_score + threshold * (max_score - min_score)

    hits: list[ScanHit] = []

    # Scan forward strand
    for pos, score in pssm.search(Seq(seq_str), threshold=score_threshold):
        pos_int = int(pos)
        score_float = float(score)
        if pos_int >= 0:
            subseq = seq_str[pos_int : pos_int + motif_length]
            hits.append(ScanHit(position=pos_int, strand="+", score=score_float, sequence=subseq))
        elif both_strands:
            rc_pos = -(pos_int + motif_length)
            subseq = seq_str[rc_pos : rc_pos + motif_length]
            hits.append(ScanHit(position=rc_pos, strand="-", score=score_float, sequence=subseq))

    hits.sort(key=lambda h: h.position)
    return hits
