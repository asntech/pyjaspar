"""DL collection CLI subcommands."""

from __future__ import annotations

import json
import sys
from dataclasses import asdict
from typing import Any

import rich_click as click

from pyjaspar import JasparDB
from pyjaspar.utils import dict_list_to_tsv


def _numpy_serializer(obj: Any) -> Any:
    """Custom JSON serializer for numpy arrays."""
    if hasattr(obj, "tolist"):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _print_output(data: Any, fmt: str, tsv_keys: list[str] | None = None) -> None:
    """Print data in the requested format."""
    if fmt == "json":
        output = json.dumps(data, indent=2, default=_numpy_serializer)
    else:
        keys = tsv_keys or []
        items = data if isinstance(data, list) else [data]
        filtered = [{k: row.get(k, "") for k in keys} for row in items]
        output = dict_list_to_tsv(filtered, keys)
    click.echo(click.style(output, fg="green"))


_PROFILE_SUMMARY_KEYS = ["profile_id", "tf_name", "tax_group"]
_MODEL_SUMMARY_KEYS = [
    "model_id",
    "tf_name",
    "cell_line",
    "tax_id",
    "data_type",
    "model_name",
    "source",
]


@click.group("dl")
def dl_group() -> None:
    """JASPAR Deep Learning (DL) collection commands.

    Access DL profiles, models, and motif patterns from JASPAR 2026.
    """


@dl_group.command("profile")
@click.argument("profile_id")
@click.option(
    "--format",
    "fmt",
    default="json",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def fetch_profile(profile_id: str, fmt: str) -> None:
    """Fetch a DL profile by ID (e.g. DL0001.1)."""
    jdb = JasparDB()
    try:
        profile = jdb.dl.fetch_profile(profile_id)
    except Exception as e:
        click.echo(click.style(str(e), fg="red"), err=True)
        sys.exit(1)

    data = asdict(profile)
    if fmt == "tsv":
        summary = {
            "profile_id": profile.profile_id,
            "tf_name": profile.tf_name,
            "tax_group": profile.tax_group,
            "primary_motif": profile.primary_motif.motif_id if profile.primary_motif else "",
            "alt_motifs": ",".join(m.motif_id for m in profile.alt_motifs),
            "jaspar_matches": ",".join(f"{k}={v}" for k, v in profile.jaspar_matches.items()),
            "linked_models": ",".join(profile.linked_model_ids),
        }
        _print_output(summary, "tsv", list(summary.keys()))
    else:
        _print_output(data, "json")


@dl_group.command("model")
@click.argument("model_id")
@click.option(
    "--format",
    "fmt",
    default="json",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def fetch_model(model_id: str, fmt: str) -> None:
    """Fetch a DL model by ID (e.g. BP000001.1)."""
    jdb = JasparDB()
    try:
        model = jdb.dl.fetch_model(model_id)
    except Exception as e:
        click.echo(click.style(str(e), fg="red"), err=True)
        sys.exit(1)

    data = asdict(model)
    if fmt == "tsv":
        summary = {
            "model_id": model.model_id,
            "tf_name": model.tf_name,
            "cell_line": model.cell_line,
            "tax_id": model.tax_id,
            "data_type": model.data_type,
            "model_name": model.model_name,
            "source": model.source,
            "source_id": model.source_id,
            "source_url": model.source_url,
            "num_motifs": len(model.motifs),
        }
        _print_output(summary, "tsv", list(summary.keys()))
    else:
        _print_output(data, "json")


@dl_group.command("motif")
@click.argument("motif_id")
@click.option(
    "--format",
    "fmt",
    default="json",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def fetch_motif(motif_id: str, fmt: str) -> None:
    """Fetch a DL motif pattern by MO* ID (e.g. MO000001.1)."""
    jdb = JasparDB()
    try:
        motif = jdb.dl.fetch_motif_pattern(motif_id)
    except Exception as e:
        click.echo(click.style(str(e), fg="red"), err=True)
        sys.exit(1)

    data = asdict(motif)
    if fmt == "tsv":
        summary = {
            "motif_id": motif.motif_id,
            "tf_name": motif.tf_name,
            "motif_type": motif.motif_type,
            "name": motif.name,
            "num_seqlets": motif.num_seqlets,
            "num_hits": motif.num_hits,
            "jaspar_match": motif.jaspar_match or "",
            "has_pfm": "PFM" in motif.matrices,
            "has_cwm": "CWM" in motif.matrices,
        }
        _print_output(summary, "tsv", list(summary.keys()))
    else:
        _print_output(data, "json")


@dl_group.command("profiles")
@click.option("-n", "--tf-name", default=None, help="TF name filter")
@click.option("--tax-group", default=None, help="Taxonomic group filter")
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def search_profiles(tf_name: str | None, tax_group: str | None, fmt: str) -> None:
    """Search DL profiles."""
    jdb = JasparDB()
    results = jdb.dl.search_profiles(tf_name=tf_name, tax_group=tax_group)

    if fmt == "json":
        data = {"count": len(results), "results": [asdict(r) for r in results]}
        _print_output(data, "json")
    else:
        rows = [asdict(r) for r in results]
        _print_output(rows, "tsv", _PROFILE_SUMMARY_KEYS)


@dl_group.command("models")
@click.option("-n", "--tf-name", default=None, help="TF name filter")
@click.option("--cell-line", default=None, help="Cell line filter")
@click.option("--tax-id", type=int, default=None, help="Taxonomy ID filter")
@click.option("--data-type", default=None, help="Data type filter")
@click.option("--model-name", default=None, help="Model name filter")
@click.option("--source", default=None, help="Source filter")
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
def search_models(
    tf_name: str | None,
    cell_line: str | None,
    tax_id: int | None,
    data_type: str | None,
    model_name: str | None,
    source: str | None,
    fmt: str,
) -> None:
    """Search DL models."""
    jdb = JasparDB()
    results = jdb.dl.search_models(
        tf_name=tf_name,
        cell_line=cell_line,
        tax_id=tax_id,
        data_type=data_type,
        model_name=model_name,
        source=source,
    )

    if fmt == "json":
        data = {"count": len(results), "results": [asdict(r) for r in results]}
        _print_output(data, "json")
    else:
        rows = [asdict(r) for r in results]
        _print_output(rows, "tsv", _MODEL_SUMMARY_KEYS)
