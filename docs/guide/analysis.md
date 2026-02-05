# Analysis Tools

Install with `pip install pyjaspar[analysis]` (adds numpy, scipy).

## Sequence scanning

Scan DNA sequences for motif occurrences using PSSM scoring.

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import scan_sequence

jdb = JasparDB()
motif = jdb.fetch_motifs_by_name("CTCF")[0]

hits = scan_sequence(
    "CCACCAGGGGGCGCAC" * 10,
    motif,
    threshold=0.7,        # fraction of max PSSM score (0.0-1.0)
    both_strands=True,    # also scan reverse complement
)

for hit in hits:
    print(f"Position {hit.position} ({hit.strand}): score={hit.score:.2f} {hit.sequence}")
```

### ScanHit fields

| Field | Type | Description |
|-------|------|-------------|
| `position` | int | 0-based start position |
| `strand` | str | `"+"` or `"-"` |
| `score` | float | PSSM score |
| `sequence` | str | Matched subsequence |

## Motif similarity

Compare two motifs using column-wise metrics.

### Pearson correlation

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import pearson_correlation, best_correlation

jdb = JasparDB()
m1 = jdb.fetch_motif_by_id("MA0001.1")
m2 = jdb.fetch_motif_by_id("MA0002.1")

# At a specific offset
score = pearson_correlation(m1, m2, offset=0)
print(f"Pearson r = {score:.4f}")

# Best alignment across all offsets and orientations
score, offset, is_rc = best_correlation(m1, m2, min_overlap=4)
print(f"Best: r={score:.4f}, offset={offset}, reverse_complement={is_rc}")
```

### Euclidean distance

```python
from pyjaspar.analysis import euclidean_distance

dist = euclidean_distance(m1, m2, offset=0)
print(f"Distance = {dist:.4f}")  # lower = more similar
```

### KL divergence

```python
from pyjaspar.analysis import kl_divergence

kl = kl_divergence(m1, m2, offset=0)
print(f"KL divergence = {kl:.4f}")  # non-negative, 0 = identical
```

### Metric comparison

| Metric | Range | Interpretation |
|--------|-------|----------------|
| Pearson correlation | [-1, 1] | 1 = identical, 0 = uncorrelated, -1 = anti-correlated |
| Euclidean distance | [0, inf) | 0 = identical, lower = more similar |
| KL divergence | [0, inf) | 0 = identical distributions |

All metrics return 0.0 (Pearson) or `float('inf')` (Euclidean, KL) when there is no overlap between motifs.

## Enrichment analysis

Test whether motifs are enriched in a foreground set of sequences compared to a background set.

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import motif_enrichment

jdb = JasparDB()
motifs = jdb.fetch_motifs(collection="CORE", tax_group="vertebrates", min_ic=15)

results = motif_enrichment(
    foreground=["ACGT..." * 50],   # your sequences of interest
    background=["TGCA..." * 50],   # control sequences
    motifs=motifs[:10],
    threshold=0.8,
)

for r in results:
    print(f"{r.motif_id} ({r.motif_name}): "
          f"p={r.pvalue:.4f}, q={r.qvalue:.4f}, "
          f"fold={r.fold_enrichment:.2f}")
```

Results are sorted by p-value with Benjamini-Hochberg correction applied.

### EnrichmentResult fields

| Field | Type | Description |
|-------|------|-------------|
| `motif_id` | str | JASPAR matrix ID |
| `motif_name` | str | TF name |
| `fg_hits` | int | Foreground sequences with at least one hit |
| `bg_hits` | int | Background sequences with at least one hit |
| `fg_total` | int | Total foreground sequences |
| `bg_total` | int | Total background sequences |
| `fold_enrichment` | float | Foreground hit rate / background hit rate |
| `pvalue` | float | Fisher's exact test p-value |
| `qvalue` | float | Benjamini-Hochberg adjusted p-value |
