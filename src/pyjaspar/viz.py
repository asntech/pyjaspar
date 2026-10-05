"""Motif visualization.

Requires the ``viz`` extra: ``pip install pyjaspar[viz]``
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

try:
    import logomaker
    import pandas as pd
    from matplotlib.ticker import MaxNLocator
except ImportError:
    raise ImportError(
        "The pyjaspar.viz module requires logomaker and matplotlib. "
        "Install with: pip install pyjaspar[viz]"
    ) from None

if TYPE_CHECKING:
    from Bio.motifs.jaspar import Motif
    from matplotlib.axes import Axes


def plot_logo(
    motif: Motif,
    ax: Axes | None = None,
    color_scheme: str = "classic",
    **logo_kwargs: Any,
) -> logomaker.Logo:
    """Draw the sequence logo of a motif.

    Letter heights are information content in bits, with the most frequent
    letter of each position on top. Positions are numbered from 1.

    Args:
        motif: Motif whose PFM counts are drawn.
        ax: Matplotlib axes to draw on. A new figure is created if None.
        color_scheme: Any logomaker color scheme name.
        **logo_kwargs: Passed through to ``logomaker.Logo``.

    Returns:
        The logomaker ``Logo``; its ``ax`` and ``fig`` attributes give access
        to the matplotlib objects.
    """
    counts = pd.DataFrame(motif.counts)
    counts.index = range(1, motif.length + 1)
    info = logomaker.transform_matrix(counts, from_type="counts", to_type="information")

    logo = logomaker.Logo(info, ax=ax, color_scheme=color_scheme, **logo_kwargs)
    logo.ax.set_ylim(0, 2)
    logo.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    logo.ax.set_xlabel("Position")
    logo.ax.set_ylabel("Information content (bits)")
    logo.style_spines(visible=False)
    return logo
