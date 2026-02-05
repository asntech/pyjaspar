# pyJASPAR

> A Pythonic interface to JASPAR transcription factor motifs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4509415.svg)](https://doi.org/10.5281/zenodo.4509415)
[![PyPI](https://img.shields.io/pypi/v/pyjaspar.svg)](https://pypi.python.org/pypi/pyjaspar)
[![Python](https://img.shields.io/pypi/pyversions/pyjaspar.svg)](https://www.python.org)
[![Bioconda](https://anaconda.org/bioconda/pyjaspar/badges/version.svg)](https://anaconda.org/bioconda/pyjaspar)
[![Issues](https://img.shields.io/github/issues/asntech/pyjaspar.svg)](https://github.com/asntech/pyjaspar/issues)

**pyJASPAR** uses BioPython and SQLite3 to provide a serverless interface to the [JASPAR database](http://jaspar.genereg.net) for querying TF motif profiles across multiple releases.

Supported releases: **JASPAR2026**, JASPAR2024, JASPAR2022, JASPAR2020, JASPAR2018, JASPAR2016, JASPAR2014.

## Installation

```bash
pip install pyjaspar
```

Or via conda/mamba:

```bash
#conda
conda install -c bioconda pyjaspar

#mamba
conda install -c bioconda pyjaspar

```

Or via Docker:

```bash
docker pull ghcr.io/asntech/pyjaspar
docker run ghcr.io/asntech/pyjaspar motif-by-id MA0095.2
```

Requires Python >= 3.10. Full documentation at [asntech.github.io/pyjaspar](https://asntech.github.io/pyjaspar).

## Quick Start — Python API

```python
from pyjaspar import JasparDB # jaspardb shoudl also work (<v5.0.0)

# Connect to the latest release (JASPAR2026)
jdb = JasparDB()

# Or specify a release
jdb = JasparDB(release="JASPAR2020")

# Context manager supported
with JasparDB() as jdb:
    motif = jdb.fetch_motif_by_id("MA0095.2")
    print(motif.name)  # YY1
```

### Fetch motifs

```python
from pyjaspar import JasparDB

jdb = JasparDB()

# By ID
motif = jdb.fetch_motif_by_id("MA0095.2")

# By name
motifs = jdb.fetch_motifs_by_name("CTCF")

# With filters
motifs = jdb.fetch_motifs(
    collection="CORE",
    tax_group=["vertebrates", "insects"],
    tf_class="Homeo domain factors",
    min_ic=12,
)

# All available releases
print(jdb.get_releases())
```

All methods return `Bio.motifs.jaspar.Motif` objects.

## Quick Start — CLI

```bash
# List releases
pyjaspar releases

# Fetch motif by ID
pyjaspar motif-by-id MA0095.2

# Fetch motifs by name
pyjaspar motifs-by-name CTCF

# Fetch CORE vertebrate motifs
pyjaspar motifs --collection CORE --tax-group Vertebrates

# Export metadata as JSON
pyjaspar metadata --tax-group Vertebrates --format json -o metadata.json

# Show help
pyjaspar --help

# Analysis via CLI
pyjaspar scan ACGTACGTACGT --motif-id MA0139.1 --threshold 0.5
pyjaspar similarity MA0001.1 MA0002.1 --best --format json
pyjaspar enrichment -fg foreground.fasta -bg background.fasta --motif-ids MA0001.1
```

## Analysis Tools

pyJASPAR provides few analysis tools, which are in development.

### Motif similarity

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import pearson_correlation, best_correlation

jdb = JasparDB()
m1 = jdb.fetch_motif_by_id("MA0001.1")
m2 = jdb.fetch_motif_by_id("MA0002.1")

# Column-wise Pearson correlation
score = pearson_correlation(m1, m2)

# Best alignment across all offsets and orientations
score, offset, is_rc = best_correlation(m1, m2)
```

### Sequence scanning

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import scan_sequence

jdb = JasparDB()
motif = jdb.fetch_motifs_by_name("CTCF")[0]

hits = scan_sequence("ACGTACGTACGTACGT" * 10, motif, threshold=0.7)
for hit in hits:
    print(f"Position {hit.position} ({hit.strand}): score={hit.score:.2f}")
```

### Enrichment analysis

```python
from pyjaspar import JasparDB
from pyjaspar.analysis import motif_enrichment

jdb = JasparDB()
motifs = jdb.fetch_motifs(collection="CORE", tax_group="vertebrates", min_ic=15)

results = motif_enrichment(
    foreground=["ACGT..." * 50],   # your sequences
    background=["TGCA..." * 50],   # control sequences
    motifs=motifs[:10],
    threshold=0.8,
)
for r in results:
    print(f"{r.motif_id} ({r.motif_name}): p={r.pvalue:.4f}, q={r.qvalue:.4f}")
```

## Deep Learning (DL) Collection

JASPAR 2026 introduces a Deep Learning collection derived from BPNet models trained on
ENCODE TF ChIP-seq data ([Rauluseviciute et al., NAR 2024](https://doi.org/10.1093/nar/gkad1059)).

Install with `pip install pyjaspar[dl]` (adds numpy).

The DL collection has three entity types:

| Entity | ID prefix | Description |
|--------|-----------|-------------|
| Profile | DL* | TF summary profile with primary and alternate motif patterns |
| Motif pattern | MO* | Cluster-level motif with PFM and CWM matrices |
| Model | BP* | BPNet model trained on a specific TF and cell line |

### Python API

```python
from pyjaspar import JasparDB

jdb = JasparDB()  # JASPAR2026

# Fetch a TF profile
profile = jdb.dl.fetch_profile("DL0001.1")
print(profile.tf_name, profile.tax_group)
print(profile.primary_motif.pfm.length)  # PFM matrix length

# Fetch a BPNet model with its motifs
model = jdb.dl.fetch_model("BP000001.1")
print(model.tf_name, model.cell_line, len(model.motifs))

# Fetch a single motif pattern
motif = jdb.dl.fetch_motif_pattern("MO000001.1")
print(motif.matrices.keys())  # dict_keys(['PFM', 'CWM'])

# Convert PFM to BioPython Motif object
bio_motif = motif.to_biopython()
print(bio_motif.pssm)

# Search profiles and models
profiles = jdb.dl.search_profiles(tf_name="REST")
models = jdb.dl.search_models(cell_line="K562")
```

The `.dl` property is only available on JASPAR2026+. Accessing it on older releases
raises `DLNotAvailableError`.

### CLI

```bash
# Fetch a profile
pyjaspar dl profile DL0001.1 --format json

# Search profiles by TF name
pyjaspar dl profiles --tf-name REST

# Fetch a model
pyjaspar dl model BP000001.1

# Search models by cell line
pyjaspar dl models --cell-line K562

# Fetch a motif pattern
pyjaspar dl motif MO000001.1 --format json
```

## Backward Compatibility

Existing code using `jaspardb` (lowercase) continues to work:

```python
from pyjaspar import jaspardb  # still works
jdb = jaspardb(release="JASPAR2026")
```

The new canonical name is `JasparDB`.

## Development

```bash
git clone https://github.com/asntech/pyjaspar.git
cd pyjaspar
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Lint
ruff check src/ tests/
ruff format src/ tests/
```

## Citation

Aziz Khan. pyJASPAR: a Pythonic interface to JASPAR transcription factor motifs. (2021). Zenodo, doi:10.5281/zenodo.4485856

