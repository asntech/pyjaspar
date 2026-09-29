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
motif = jdb.fetch_motif_by_id("MA0139.1")

logo = plot_logo(motif)
logo.ax.set_title(f"{motif.matrix_id} {motif.name}")
logo.fig.savefig("ctcf_logo.png", dpi=150, bbox_inches="tight")
```

Pass `ax=` to draw several logos on one figure:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 1, sharex=True)
for ax, motif in zip(axes, jdb.fetch_motifs_by_name("ATF3")[:2]):
    plot_logo(motif, ax=ax)
    ax.set_title(f"{motif.matrix_id} {motif.name}")
```
