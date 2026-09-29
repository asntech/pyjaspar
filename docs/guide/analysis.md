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

## Profile inference

Predict which JASPAR profiles a protein binds from its amino acid sequence
(full length or just the DNA-binding domain). The protein is compared with the
DNA-binding domains of TFs that already have a JASPAR profile, so the hits are
its closest *relatives*. The search runs on a JASPAR server, so it needs
network access and takes several seconds. The method is described in the
[JASPAR documentation](https://jaspar.elixir.no/docs/) and the
[JASPAR 2016 paper](https://doi.org/10.1093/nar/gkv1176); the tool itself is
[on GitHub](https://github.com/wassermanlab/JASPAR-inference-tool).

```python
from pyjaspar.analysis import infer_profiles

dbd = "ACPVETCDRRFSRSDELTRHIRIHTGQKPFQCRICMRNFSRSDHLTTHIRTHTGEKPFACEICGRKFARSDERKRHTKIHMRQKDKKAEKGA"  # a zinc-finger domain

hits = infer_profiles(dbd)

for hit in hits:
    print(f"{hit.matrix_id}  {hit.name:<5}  E-value={hit.evalue:<9.2e}  DBD identity={hit.dbd_identity:.0%}")
# MA0162.2  EGR1   E-value=1.90e-63   DBD identity=94%
# MA0732.1  EGR3   E-value=2.34e-59   DBD identity=90%
# MA0472.2  EGR2   E-value=1.50e-54   DBD identity=93%
# MA0733.1  EGR4   E-value=9.15e-51   DBD identity=77%
```

The server returns matrix IDs and logo URLs, not the matrices themselves.
`JASPAR2024` is searched by default; pass `release=` to search another. The
matrix IDs belong to that release, so fetch the profile from the same one:

```python
motif = JasparDB("JASPAR2024").fetch_motif_by_id(hits[0].matrix_id)
print(motif.consensus)  # CCCCCGCCCCCGCC
```

An empty list means the server found no matching profile.

### InferenceHit fields

| Field | Type | Description |
|-------|------|-------------|
| `matrix_id` | str | JASPAR matrix ID |
| `name` | str | Name of the TF the profile belongs to (a relative of your protein) |
| `evalue` | float | E-value of the DNA-binding domain match; lower is more significant |
| `dbd_identity` | float | Share of identical amino acids in the DNA-binding domain, as a fraction (0-1) |
| `logo_url` | str | URL of the profile's sequence logo (SVG) |
| `release` | str | JASPAR release the matrix ID belongs to |
