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

## Profile search

`search_profiles()` scores a query motif against a set of profiles and ranks them,
like the **Matrix Align** tool on the JASPAR web site. Build the set with
`fetch_motifs`, which has the same filters as the web form (collection, taxonomic
group, latest or all versions).

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import search_profiles

jdb = JasparDB()
query = jdb.fetch_motif_by_id("MA0139.2")  # CTCF

candidates = jdb.fetch_motifs(
    collection=["CORE"], tax_group=["Vertebrates"], all_versions=False
)
for hit in search_profiles(query, candidates, top=3):
    print(f"{hit.matrix_id}  {hit.name:<8}  score={hit.score:.4f}")
# MA0139.2  CTCF      score=30.0000
# MA1930.2  CTCF      score=29.8880
# MA1929.2  CTCF      score=27.8202
```

`score` is the web tool's "Score" column. Every aligned column contributes between
0 and 2, so a profile aligned with itself scores twice its width, and a longer
profile tends to score higher. The alignment may leave columns hanging off either
end for free and may contain one internal gap (`open_penalty=3.0` for its first
column, `ext_penalty=0.01` for each further column); both the candidate and its
reverse complement are tried.

`percent_score` reproduces the web tool's "Percent Score" column: `100 * score / (2 * m)`,
where `m` is the narrowest profile seen so far. `m` starts at the query's width and is
lowered by each candidate in matrix-ID order, the order of the web tool's table, and is
never raised again. The value therefore depends on the set of candidates, not only on the
pair, and it can exceed 100. Use `sort_by="percent_score"` to rank by it.

The method is a semi-global variant of the Needleman-Wunsch algorithm that permits
one internal gap, as documented for
[`TFBS::Matrix::Alignment`](https://manpages.ubuntu.com/manpages/noble/man3/TFBS::Matrix::Alignment.3pm.html)
(Sandelin et al., [Funct Integr Genomics 3:125-134, 2003](https://doi.org/10.1007/s10142-003-0086-6)).

### ProfileHit fields

| Field | Type | Description |
|-------|------|-------------|
| `matrix_id` | str | JASPAR matrix ID of the candidate |
| `name` | str | Name of the TF the candidate belongs to |
| `score` | float | Alignment score (the web tool's "Score") |
| `percent_score` | float | The web tool's "Percent Score" (depends on the set of candidates, see above) |
| `is_reverse_complement` | bool | The candidate's reverse complement aligned better |
| `gaps` | int | Internal gap runs in the best alignment (0 or 1) |
| `width` | int | Number of columns of the candidate |

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
