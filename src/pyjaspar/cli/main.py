"""pyJASPAR command-line interface.

Built with rich-click for enhanced terminal output.
This is a thin wrapper around the pyjaspar core API.
"""

from __future__ import annotations

import json
import sys
from typing import Any

import rich_click as click

from pyjaspar import (
    JASPAR_DFLT_COLLECTION,
    JASPAR_LATEST_RELEASE,
    JasparDB,
    jaspar_releases,
)
from pyjaspar.utils import dict_list_to_tsv

# --- rich-click configuration ---

click.rich_click.COMMAND_GROUPS = {
    "pyjaspar": [
        {
            "name": "Motif retrieval commands",
            "commands": ["motifs", "motif-by-id", "motifs-by-name", "metadata"],
        },
        {
            "name": "Analysis commands",
            "commands": ["scan", "similarity", "enrichment"],
        },
        {
            "name": "Deep Learning (DL) collection",
            "commands": ["dl"],
        },
        {
            "name": "General commands",
            "commands": ["releases", "collections", "cite"],
        },
    ]
}

TAXONOMIC_GROUPS = ["Fungi", "Insects", "Nematodes", "Plants", "Urochordates", "Vertebrates"]

_RELEASE_YEARS = [key.replace("JASPAR", "") for key in jaspar_releases]
_LATEST_YEAR = JASPAR_LATEST_RELEASE.replace("JASPAR", "").strip()

_METADATA_KEYS = [
    "matrix_id",
    "name",
    "collection",
    "length",
    "tf_class",
    "tf_family",
    "species",
    "tax_group",
    "acc",
    "data_type",
    "medline",
    "pazar_id",
    "comment",
]


# --- Helper functions ---


def _print_matrix_data(
    motifs: Any,
    output_file: str | None = None,
    motif_format: str = "jaspar",
) -> None:
    """Format and output motif matrices."""
    motifs_output = ""
    for motif in motifs:
        motifs_output += str(motif.format(motif_format))
    if output_file:
        with open(output_file, "w") as f:
            f.write(motifs_output)
    else:
        click.echo(click.style(motifs_output, fg="green"))


def _print_metadata(
    motifs: Any,
    output_file: str | None = None,
    fmt: str = "tsv",
) -> None:
    """Format and output motif metadata."""
    motifs_list: list[dict[str, Any]] = []
    for motif in motifs:
        motif_dict = vars(motif)
        motifs_list.append({k: motif_dict[k] for k in _METADATA_KEYS if k in motif_dict})

    if fmt == "json":
        motifs_dict = {"count": len(motifs_list), "results": motifs_list}
        output = json.dumps(motifs_dict, indent=2)
    else:
        output = dict_list_to_tsv(motifs_list, _METADATA_KEYS)

    if output_file:
        with open(output_file, "w") as f:
            f.write(output)
    else:
        click.echo(click.style(output, fg="green"))


# --- CLI commands ---


@click.group()
@click.version_option(package_name="pyjaspar", prog_name="pyJASPAR")
def cli() -> None:
    """pyJASPAR: A Pythonic interface to JASPAR transcription factor motifs.

    https://github.com/asntech/pyjaspar
    """


@cli.command("motif-by-id")
@click.argument("id")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "--motif-format",
    default="jaspar",
    show_default=True,
    type=click.Choice(["pfm", "jaspar", "transfac"], case_sensitive=False),
    help="Motif output format",
)
@click.option("--metadata", is_flag=True, default=False, help="Return metadata instead of matrix")
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Metadata output format",
)
def motif_by_id(id: str, release: str, motif_format: str, metadata: bool, fmt: str) -> None:
    """Get motif matrix by JASPAR ID."""
    jdb = JasparDB(f"JASPAR{release}")
    motif = jdb.fetch_motif_by_id(id)
    if motif is None:
        click.echo(click.style(f"No motif found with ID: {id}", fg="red"), err=True)
        sys.exit(1)

    if metadata:
        _print_metadata([motif], fmt=fmt)
    else:
        _print_matrix_data([motif], motif_format=motif_format)


@cli.command("motifs-by-name")
@click.argument("tf_name")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "--motif-format",
    default="jaspar",
    show_default=True,
    type=click.Choice(["pfm", "jaspar", "transfac"], case_sensitive=False),
    help="Motif output format",
)
@click.option("--metadata", is_flag=True, default=False, help="Return metadata instead of matrix")
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Metadata output format",
)
def motifs_by_name(tf_name: str, release: str, motif_format: str, metadata: bool, fmt: str) -> None:
    """Get motif(s) by transcription factor name.

    TF_NAME is the transcription factor name (e.g. CTCF, YY1).
    """
    jdb = JasparDB(f"JASPAR{release}")
    motifs = jdb.fetch_motifs_by_name(tf_name)

    if metadata:
        _print_metadata(motifs, fmt=fmt)
    else:
        _print_matrix_data(motifs, motif_format=motif_format)


@cli.command("motifs")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "-c",
    "--collection",
    default=JASPAR_DFLT_COLLECTION,
    show_default=True,
    help="JASPAR collection(s)",
)
@click.option("-n", "--tf-name", default=None, help="TF name")
@click.option("--tf-class", default=None, help="TF class(es)")
@click.option("--tf-family", default=None, help="TF family(ies)")
@click.option("--matrix-id", default=None, help="JASPAR matrix ID(s)")
@click.option(
    "--tax-group",
    default=None,
    type=click.Choice(TAXONOMIC_GROUPS, case_sensitive=False),
    help="Taxonomic supergroup",
)
@click.option("--species", default=None, help="Species (taxonomy IDs)")
@click.option("--pazar-id", default=None, help="PAZAR TF ID(s)")
@click.option("--data-type", default=None, help="Data type(s)")
@click.option("--medline", default=None, help="Medline (PubMed IDs)")
@click.option("--min-ic", type=float, default=None, help="Minimum information content")
@click.option("--min-length", type=int, default=None, help="Minimum motif length")
@click.option("--min-sites", type=int, default=None, help="Minimum binding sites")
@click.option("-o", "--output-file", default=None, help="Save output to file")
@click.option(
    "--motif-format",
    default="jaspar",
    show_default=True,
    type=click.Choice(["pfm", "jaspar", "transfac"], case_sensitive=False),
    help="Motif output format",
)
@click.option("--redundant", is_flag=True, default=False, help="Return all versions of motifs")
def motifs(
    release: str,
    collection: str,
    tf_name: str | None,
    tf_class: str | None,
    tf_family: str | None,
    matrix_id: str | None,
    tax_group: str | None,
    species: str | None,
    pazar_id: str | None,
    data_type: str | None,
    medline: str | None,
    min_ic: float | None,
    min_length: int | None,
    min_sites: int | None,
    redundant: bool,
    motif_format: str,
    output_file: str | None,
) -> None:
    """Get JASPAR motifs in different formats."""
    jdb = JasparDB(f"JASPAR{release}")
    result = jdb.fetch_motifs(
        collection=collection,
        tf_name=tf_name,
        tf_class=tf_class,
        tf_family=tf_family,
        matrix_id=matrix_id,
        tax_group=tax_group,
        species=species,
        pazar_id=pazar_id,
        data_type=data_type,
        medline=medline,
        min_ic=min_ic or 0,
        min_length=min_length or 0,
        min_sites=min_sites or 0,
        all_versions=redundant,
    )
    _print_matrix_data(result, output_file=output_file, motif_format=motif_format)


@cli.command("metadata")
@click.option(
    "-c",
    "--collection",
    default=JASPAR_DFLT_COLLECTION,
    show_default=True,
    help="JASPAR collection",
)
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
@click.option(
    "--tax-group",
    default=None,
    type=click.Choice(TAXONOMIC_GROUPS, case_sensitive=False),
    help="Taxonomic supergroup",
)
@click.option("--species", default=None, help="Species (taxonomy IDs)")
@click.option(
    "--format",
    "fmt",
    default="tsv",
    show_default=True,
    type=click.Choice(["json", "tsv"], case_sensitive=False),
    help="Output format",
)
@click.option("-o", "--output-file", default=None, help="Save output to file")
@click.option("--redundant", is_flag=True, default=False, help="Return all versions of motifs")
def metadata(
    collection: str,
    release: str,
    tax_group: str | None,
    species: str | None,
    fmt: str,
    output_file: str | None,
    redundant: bool,
) -> None:
    """Get metadata for motif(s)."""
    jdb = JasparDB(f"JASPAR{release}")
    result = jdb.fetch_motifs(
        collection=collection,
        tax_group=tax_group,
        species=species,
        all_versions=redundant,
    )
    _print_metadata(result, output_file=output_file, fmt=fmt)


@cli.command("releases")
@click.option("-l", "--latest", is_flag=True, default=False, help="Show only the latest release")
def releases(latest: bool) -> None:
    """Get available JASPAR releases."""
    if latest:
        click.echo(JASPAR_LATEST_RELEASE)
    else:
        for release in jaspar_releases:
            click.echo(release)


@cli.command("collections")
@click.option(
    "-r",
    "--release",
    default=_LATEST_YEAR,
    show_default=True,
    type=click.Choice(_RELEASE_YEARS, case_sensitive=False),
    help="JASPAR release year",
)
def collections(release: str) -> None:
    """Get JASPAR collections available in a release.

    CORE profiles are curated with orthogonal supporting evidence.
    UNVALIDATED profiles are computationally sound but not yet
    independently validated.
    """
    with JasparDB(f"JASPAR{release}") as jdb:
        for collection in jdb.get_collections():
            click.echo(collection)


@cli.command("cite")
def cite() -> None:
    """Show citation for JASPAR."""
    click.echo(
        click.style(
            "Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, Castro-Mondragon JA, "
            "Ferenc K, Kumar V, Lemma RB, Lucas J, Chèneby J, Baranasic D, "
            "Khan A, Fornes O, Gundersen S, Johansen M, Hovig E, Lenhard B, "
            "Sandelin A, Wasserman WW, Parcy F, Mathelier A. "
            "JASPAR 2024: 20th anniversary of the open-access database of "
            "transcription factor binding profiles. "
            "Nucleic Acids Res. 2024 Jan 5;52(D1):D174-D182.",
            fg="green",
        )
    )


from .analysis import enrichment, scan, similarity  # noqa: E402
from .dl import dl_group  # noqa: E402

#cli.add_command(scan)
#cli.add_command(similarity)
#cli.add_command(enrichment)

cli.add_command(dl_group)

if __name__ == "__main__":
    cli()
