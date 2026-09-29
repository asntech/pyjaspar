"""Visualization CLI subcommands."""

from __future__ import annotations

import sys

import rich_click as click

from pyjaspar import JASPAR_LATEST_RELEASE, JasparDB, jaspar_releases

_RELEASE_YEARS = [key.replace("JASPAR", "") for key in jaspar_releases]
_LATEST_YEAR = JASPAR_LATEST_RELEASE.replace("JASPAR", "").strip()


def _require_viz() -> None:
    """Check that the viz extra is installed."""
    try:
        import pyjaspar.viz  # noqa: F401
    except ImportError:
        click.echo(
            click.style(
                "The logo command requires the 'viz' extra. "
                "Install with: pip install pyjaspar[viz]",
                fg="red",
            ),
            err=True,
        )
        sys.exit(1)


@click.command("logo")
@click.argument("motif_id")
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Output image file; the format follows the extension (.png, .svg, .pdf)",
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
    "--dpi", default=150, show_default=True, type=int, help="Resolution for raster output"
)
def logo(motif_id: str, output: str, release: str, dpi: int) -> None:
    """Draw the sequence logo of a motif to an image file."""
    _require_viz()
    import matplotlib

    matplotlib.use("Agg")
    from pyjaspar.viz import plot_logo

    with JasparDB(f"JASPAR{release}") as jdb:
        motif = jdb.fetch_motif_by_id(motif_id)
    if motif is None:
        click.echo(click.style(f"No motif found with ID: {motif_id}", fg="red"), err=True)
        sys.exit(1)

    drawn = plot_logo(motif)
    drawn.ax.set_title(f"{motif.matrix_id} {motif.name}")
    drawn.fig.savefig(output, dpi=dpi, bbox_inches="tight")
    click.echo(click.style(f"Saved {output}", fg="green"))
