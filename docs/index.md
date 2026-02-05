# pyJASPAR

> A Pythonic interface to JASPAR transcription factor motifs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4509415.svg)](https://doi.org/10.5281/zenodo.4509415)
[![PyPI](https://img.shields.io/pypi/v/pyjaspar.svg)](https://pypi.python.org/pypi/pyjaspar)
[![Python](https://img.shields.io/pypi/pyversions/pyjaspar.svg)](https://www.python.org)
[![Bioconda](https://anaconda.org/bioconda/pyjaspar/badges/version.svg)](https://anaconda.org/bioconda/pyjaspar)
[![License](https://img.shields.io/github/license/asntech/pyjaspar)](https://github.com/asntech/pyjaspar/blob/main/LICENSE)

**pyJASPAR** uses BioPython and SQLite3 to provide a serverless interface to the [JASPAR database](http://jaspar.genereg.net) for querying TF motif profiles across multiple releases.

Supported releases: **JASPAR2026**, JASPAR2024, JASPAR2022, JASPAR2020, JASPAR2018, JASPAR2016, JASPAR2014.

## Quick Install

```bash
pip install pyjaspar[all]
```

## Features

<div class="grid cards" markdown>

-   **Core API**

    ---

    Query JASPAR motifs by ID, name, or advanced filters. Returns BioPython `Motif` objects.

    [:octicons-arrow-right-24: Core API](guide/core-api.md)

-   **Analysis Tools**

    ---

    Sequence scanning, motif similarity (Pearson, Euclidean, KL), and enrichment analysis.

    [:octicons-arrow-right-24: Analysis](guide/analysis.md)

-   **DL Collection**

    ---

    JASPAR 2026 Deep Learning collection: BPNet profiles, models, and motif patterns.

    [:octicons-arrow-right-24: DL Collection](guide/dl-collection.md)

-   **CLI**

    ---

    Command-line interface for motif retrieval, analysis, and DL queries.

    [:octicons-arrow-right-24: CLI Reference](guide/cli.md)

</div>

## Quick Example

```python
from pyjaspar import JasparDB

jdb = JasparDB()  # defaults to JASPAR2026

# Fetch a motif by ID
motif = jdb.fetch_motif_by_id("MA0095.2")
print(motif.name)  # YY1

# Fetch motifs with filters
motifs = jdb.fetch_motifs(
    collection="CORE",
    tax_group="vertebrates",
    min_ic=12,
)
print(f"Found {len(motifs)} motifs")
```

## Next Steps

- [Installation](getting-started/installation.md) -- all install options
- [Quick Start](getting-started/quickstart.md) -- get up and running
- [CLI Reference](guide/cli.md) -- command-line usage
