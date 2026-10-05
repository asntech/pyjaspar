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

`search_profiles()` scores a query motif against a set of profiles and ranks them. Build the set with `fetch_motifs`, which has the same filters as the web form of the **Matrix Align** tool on the JASPAR web site (collection, taxonomic group, latest or all versions).

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

### Scoring

With the default `method="matrix_align"`, the scoring implemented here is the one of the Matrix Align web tool (`method="pearson"` ranks by the best Pearson correlation instead, see [Motif alignment](#motif-alignment)). `score` is its "Score" column. Every aligned column contributes between 0 and 2, so a profile aligned with itself scores twice its width, and a longer profile tends to score higher. The alignment may leave columns hanging off either end for free and may contain one internal gap (`open_penalty=3.0` for its first column, `ext_penalty=0.01` for each further column); both the candidate and its reverse complement are tried, and the reverse complement is reported when both score the same.

`percent_score` reproduces the web tool's "Percent Score" column: `100 * score / (2 * m)`, where `m` is the narrowest profile seen so far. `m` starts at the query's width and is lowered by each candidate in matrix-ID order, the order of the web tool's table, and is never raised again. The value therefore depends on the set of candidates, not only on the pair, and it can exceed 100. Use `sort_by="percent_score"` to rank by it.

The alignment is a semi-global variant of the Needleman-Wunsch algorithm that permits one internal gap, as implemented by the `matrix_aligner` program that the web tool runs (Sandelin et al., [Funct Integr Genomics 3:125-134, 2003](https://doi.org/10.1007/s10142-003-0086-6); source in the [`jaspar_tools`](https://bitbucket.org/CBGR/jaspar_tools) repository). The implementation here is independent; `gaps`, `offset` and `alignment_length` follow what that program reports.

### ProfileHit fields

The fields after `percent_score` are read from `result`, the comparison of the query with the candidate (the same object `align_motifs` returns).

| Field | Type | Description |
|-------|------|-------------|
| `matrix_id` | str | JASPAR matrix ID of the candidate |
| `name` | str | Name of the TF the candidate belongs to |
| `result` | `AlignmentResult` | The comparison of the query with this candidate |
| `score` | float | Alignment score (the web tool's "Score") |
| `percent_score` | float | The web tool's "Percent Score" (depends on the set of candidates, see above); None with `method="pearson"` |
| `is_reverse_complement` | bool | The candidate's reverse complement aligned better |
| `gaps` | int | Gap columns in the best alignment (0 if it has no gap); the gap is one run of that many columns |
| `width` | int | Number of columns of the candidate |
| `offset` | int | Position in the query minus position in the candidate at the first pair of aligned columns (from 1, in the orientation of the candidate that was aligned) |
| `alignment_length` | int | Columns of the alignment, gap columns included and the free overhangs excluded |

## Motif similarity

Compare two motifs using column-wise metrics.

### Pearson correlation

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import align_motifs, pearson_correlation

jdb = JasparDB()
m1 = jdb.fetch_motif_by_id("MA0001.1")
m2 = jdb.fetch_motif_by_id("MA0002.1")

# At a specific offset
score = pearson_correlation(m1, m2, offset=0)
print(f"Pearson r = {score:.4f}")

# Best alignment across all offsets and orientations
best = align_motifs(m1, m2, method="pearson", min_overlap=4)
print(f"Best: r={best.score:.4f}, offset={best.offset}, "
      f"reverse_complement={best.is_reverse_complement}")
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

## Motif alignment

`align_motifs()` compares two motifs with `method="matrix_align"` (default, the Matrix Align scoring described under [Profile search](#profile-search)) or `method="pearson"` (the best Pearson correlation over all offsets and both orientations), and returns the score and the aligned columns. `format_alignment()` renders the result as a gapped, side-by-side alignment of the two consensus sequences. Scores of the two methods are on different scales.

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import align_motifs, format_alignment

jdb = JasparDB()
m1 = jdb.fetch_motif_by_id("MA0001.1")
m2 = jdb.fetch_motif_by_id("MA0002.1")

result = align_motifs(m1, m2, method="pearson")
print(format_alignment(m1, m2, result))
# target            0 ------CCATAAATAG 10
#                   0 ------|.|||----- 16
# query             0 TAACCACAATA----- 11

print(f"score={result.score:.4f} offset={result.offset} "
      f"is_reverse_complement={result.is_reverse_complement}")
```

### Result fields

`align_motifs` returns an `AlignmentResult`.

| Field | Type | Description |
|-------|------|-------------|
| `motif1_id`, `motif2_id` | str | The two motifs' JASPAR matrix IDs |
| `score` | float | Matrix Align score, or the Pearson correlation |
| `method` | str | `"matrix_align"` or `"pearson"` |
| `alignment` | `AlignmentPath` | The aligned columns: zero-based column boundaries of both motifs (`coordinates`), without the columns that hang off either end |
| `offset` | int | Position in motif1 minus position in motif2 at the first aligned pair |
| `is_reverse_complement` | bool | Whether motif2 was reverse-complemented |
| `overlap`, `gaps`, `alignment_length` | int | Paired columns, gap columns (Pearson: always 0), and their sum |
| `parameters` | dict | The options used |

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
its closest *relatives*. The search runs on a JASPAR server (the
[JASPAR profile inference tool](https://github.com/wassermanlab/JASPAR-inference-tool)),
so it needs network access and takes several seconds.

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
