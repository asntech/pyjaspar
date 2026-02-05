"""Analysis CLI subcommands."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import rich_click as click

from pyjaspar import JASPAR_LATEST_RELEASE, JasparDB, jaspar_releases
from pyjaspar.utils import dict_list_to_tsv

_RELEASE_YEARS = [key.replace("JASPAR", "") for key in jaspar_releases]
_LATEST_YEAR = JASPAR_LATEST_RELEASE.replace("JASPAR", "").strip()

_SCAN_TSV_KEYS = ["seq_id", "position", "strand", "score", "sequence"]
_SIMILARITY_TSV_KEYS = ["motif1", "motif2", "metric", "offset", "score"]
_SIMILARITY_BEST_TSV_KEYS = [
    "motif1",
    "motif2",
    "metric",
    "best_score",
    "best_offset",
    "is_reverse_complement",
]
_ENRICHMENT_TSV_KEYS = [
    "motif_id",
    "motif_name",
    "fg_hits",
    "bg_hits",
    "fg_total",
    "bg_total",
    "fold_enrichment",
    "pvalue",
    "qvalue",
]


def _print_output(data: Any, fmt: str, tsv_keys: list[str] | None = None) -> None:
    """Print data in the requested format."""
    if fmt == "json":
        output = json.dumps(data, indent=2)
    else:
        keys = tsv_keys or []
        items = data if isinstance(data, list) else [data]
        filtered = [{k: row.get(k, "") for k in keys} for row in items]
        output = dict_list_to_tsv(filtered, keys)
    click.echo(click.style(output, fg="green"))


def _is_fasta_file(value: str) -> bool:
    """Check if the value looks like a FASTA file path."""
    p = Path(value)
    return p.exists() and p.suffix.lower() in (".fa", ".fasta", ".fna", ".fas")


def _read_fasta(path: str) -> list[tuple[str, str]]:
    """Read a FASTA file and return list of (header, sequence) tuples."""
    from Bio import SeqIO

    return [(record.id, str(record.seq)) for record in SeqIO.parse(path, "fasta")]


def _require_analysis() -> None:
    """Check that the analysis extra is installed."""
    try:
        import pyjaspar.analysis  # noqa: F401
    except ImportError:
        click.echo(
            click.style(
                "Analysis commands require the 'analysis' extra. "
                "Install with: pip install pyjaspar[analysis]",
                fg="red",
            ),
            err=True,
        )
        sys.exit(1)


@click.command("scan")
@click.argument("sequence")
@click.option("--motif-id", required=True, help="JASPAR motif ID (e.g. MA0095.2)")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "--threshold",
    default=0.8,
    show_default=True,
    type=float,
    help="Score threshold (0.0-1.0)",
)
@click.option("--no-reverse", is_flag=True, default=False, help="Only scan forward strand")
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def scan(
    sequence: str,
    motif_id: str,
    release: str,
    threshold: float,
    no_reverse: bool,
    fmt: str,
) -> None:
    """Scan a DNA sequence for motif occurrences.

    SEQUENCE can be a literal DNA string or a path to a FASTA file
    (.fa, .fasta, .fna, .fas).
    """
    _require_analysis()
    from pyjaspar.analysis import scan_sequence

    jdb = JasparDB(f"JASPAR{release}")
    motif = jdb.fetch_motif_by_id(motif_id)
    if motif is None:
        click.echo(click.style(f"No motif found with ID: {motif_id}", fg="red"), err=True)
        sys.exit(1)

    # Determine if input is a FASTA file or literal DNA sequence
    sequences = _read_fasta(sequence) if _is_fasta_file(sequence) else [("input", sequence)]

    all_hits: list[dict[str, Any]] = []
    for seq_id, seq_str in sequences:
        hits = scan_sequence(seq_str, motif, threshold=threshold, both_strands=not no_reverse)
        for hit in hits:
            all_hits.append(
                {
                    "seq_id": seq_id,
                    "position": hit.position,
                    "strand": hit.strand,
                    "score": round(hit.score, 4),
                    "sequence": hit.sequence,
                }
            )

    if fmt == "json":
        data = {
            "motif_id": motif_id,
            "threshold": threshold,
            "count": len(all_hits),
            "results": all_hits,
        }
        _print_output(data, "json")
    else:
        _print_output(all_hits, "tsv", _SCAN_TSV_KEYS)


@click.command("similarity")
@click.argument("id1")
@click.argument("id2")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "--metric",
    default="pearson",
    show_default=True,
    type=click.Choice(["pearson", "euclidean", "kl"], case_sensitive=False),
    help="Similarity metric",
)
@click.option(
    "--best",
    "find_best",
    is_flag=True,
    default=False,
    help="Find best alignment (Pearson only)",
)
@click.option(
    "--offset",
    default=0,
    show_default=True,
    type=int,
    help="Alignment offset (ignored with --best)",
)
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def similarity(
    id1: str,
    id2: str,
    release: str,
    metric: str,
    find_best: bool,
    offset: int,
    fmt: str,
) -> None:
    """Compare two motifs using similarity metrics.

    Computes similarity between two JASPAR motifs (ID1, ID2).
    """
    _require_analysis()
    from pyjaspar.analysis import (
        best_correlation,
        euclidean_distance,
        kl_divergence,
        pearson_correlation,
    )

    jdb = JasparDB(f"JASPAR{release}")
    motif1 = jdb.fetch_motif_by_id(id1)
    if motif1 is None:
        click.echo(click.style(f"No motif found with ID: {id1}", fg="red"), err=True)
        sys.exit(1)

    motif2 = jdb.fetch_motif_by_id(id2)
    if motif2 is None:
        click.echo(click.style(f"No motif found with ID: {id2}", fg="red"), err=True)
        sys.exit(1)

    if find_best:
        score, best_offset, is_rc = best_correlation(motif1, motif2)
        result = {
            "motif1": id1,
            "motif2": id2,
            "metric": "pearson",
            "best_score": round(score, 6),
            "best_offset": best_offset,
            "is_reverse_complement": is_rc,
        }
        _print_output(result, fmt, _SIMILARITY_BEST_TSV_KEYS)
    else:
        metric_funcs = {
            "pearson": pearson_correlation,
            "euclidean": euclidean_distance,
            "kl": kl_divergence,
        }
        func = metric_funcs[metric]
        score = func(motif1, motif2, offset=offset)
        result = {
            "motif1": id1,
            "motif2": id2,
            "metric": metric,
            "offset": offset,
            "score": round(score, 6),
        }
        _print_output(result, fmt, _SIMILARITY_TSV_KEYS)


@click.command("enrichment")
@click.option(
    "--foreground",
    "-fg",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of foreground sequences",
)
@click.option(
    "--background",
    "-bg",
    required=True,
    type=click.Path(exists=True),
    help="FASTA file of background sequences",
)
@click.option("--motif-ids", required=True, help="Comma-separated JASPAR motif IDs")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "--threshold",
    default=0.8,
    show_default=True,
    type=float,
    help="Score threshold (0.0-1.0)",
)
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def enrichment(
    foreground: str,
    background: str,
    motif_ids: str,
    release: str,
    threshold: float,
    fmt: str,
) -> None:
    """Test motif enrichment in foreground vs background sequences.

    Reads FASTA files and tests whether motifs are enriched in
    foreground sequences compared to background using Fisher's exact test.
    """
    _require_analysis()
    from dataclasses import asdict

    from pyjaspar.analysis import motif_enrichment

    jdb = JasparDB(f"JASPAR{release}")

    # Parse motif IDs and fetch motifs
    ids = [mid.strip() for mid in motif_ids.split(",")]
    motifs = []
    for mid in ids:
        motif = jdb.fetch_motif_by_id(mid)
        if motif is None:
            click.echo(
                click.style(f"Warning: motif {mid} not found, skipping", fg="yellow"),
                err=True,
            )
        else:
            motifs.append(motif)

    if not motifs:
        click.echo(click.style("No valid motifs found", fg="red"), err=True)
        sys.exit(1)

    # Read sequences
    fg_seqs = [seq for _, seq in _read_fasta(foreground)]
    bg_seqs = [seq for _, seq in _read_fasta(background)]

    results = motif_enrichment(fg_seqs, bg_seqs, motifs, threshold=threshold)

    rows = []
    for r in results:
        d = asdict(r)
        d["pvalue"] = round(d["pvalue"], 6)
        d["fold_enrichment"] = round(d["fold_enrichment"], 4)
        if d["qvalue"] is not None:
            d["qvalue"] = round(d["qvalue"], 6)
        rows.append(d)

    if fmt == "json":
        data = {"count": len(rows), "threshold": threshold, "results": rows}
        _print_output(data, "json")
    else:
        _print_output(rows, "tsv", _ENRICHMENT_TSV_KEYS)
