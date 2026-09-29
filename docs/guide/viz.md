# Visualization

Install with `pip install pyjaspar[viz]` (adds logomaker, matplotlib, pandas).

## Sequence logo

`plot_logo()` draws the sequence logo of a motif. Letter heights are
information content in bits, with the most frequent letter of each position
on top.

```python
from pyjaspar import JasparDB
from pyjaspar.viz import plot_logo

jdb = JasparDB()
motif = jdb.fetch_motif_by_id("MA0139.2")  # CTCF, human

logo = plot_logo(motif)
logo.ax.set_title(f"{motif.matrix_id} {motif.name}")
logo.fig.savefig("ctcf_logo.png", dpi=150, bbox_inches="tight")
```

Pass `ax=` to draw several logos on one figure:

```python
import matplotlib.pyplot as plt

# CTCF in vertebrates vs insects
motifs = [jdb.fetch_motif_by_id(i) for i in ("MA0139.2", "MA0531.2")]

fig, axes = plt.subplots(2, 1, sharex=True)
for ax, motif in zip(axes, motifs):
    plot_logo(motif, ax=ax)
    ax.set_title(f"{motif.matrix_id} {motif.name} ({motif.tax_group})")
```
