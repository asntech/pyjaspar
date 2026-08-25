# CLI Reference


## Motif retrieval

### motif-by-id

Fetch a single motif by JASPAR matrix ID.

```bash
pyjaspar motif-by-id MA0095.2

# Different output formats
pyjaspar motif-by-id MA0095.2 --motif-format pfm
pyjaspar motif-by-id MA0095.2 --motif-format transfac

# Metadata instead of matrix
pyjaspar motif-by-id MA0095.2 --metadata
pyjaspar motif-by-id MA0095.2 --metadata --format json

# Specific release
pyjaspar motif-by-id MA0095.2 --release 2020
```

### motifs-by-name

Fetch motifs by transcription factor name.

```bash
pyjaspar motifs-by-name CTCF
pyjaspar motifs-by-name CTCF --metadata --format json
```

### motifs

Fetch motifs with advanced filters.

```bash
# CORE vertebrate motifs
pyjaspar motifs --collection CORE --tax-group Vertebrates

# Save to file
pyjaspar motifs --collection CORE --tax-group Vertebrates -o motifs.jaspar

# PFM format
pyjaspar motifs --collection CORE --motif-format pfm

# Filter by minimum information content
pyjaspar motifs --min-ic 12 --min-length 8
```

### metadata

Export motif metadata.

```bash
pyjaspar metadata --tax-group Vertebrates
pyjaspar metadata --tax-group Vertebrates --format json
pyjaspar metadata --tax-group Vertebrates --format json -o metadata.json
```

## Analysis commands

Requires `pip install pyjaspar[analysis]`.

### scan

Scan a DNA sequence for motif occurrences.

```bash
# Literal DNA sequence
pyjaspar scan ACGTACGTACGTACGT --motif-id MA0139.1 --threshold 0.5

# FASTA file input
pyjaspar scan sequences.fasta --motif-id MA0139.1

# JSON output
pyjaspar scan ACGTACGT --motif-id MA0139.1 --format json

# Forward strand only
pyjaspar scan ACGTACGT --motif-id MA0139.1 --no-reverse
```

### similarity

Compare two motifs using similarity metrics.

```bash
# Pearson correlation (default)
pyjaspar similarity MA0001.1 MA0002.1

# Euclidean distance
pyjaspar similarity MA0001.1 MA0002.1 --metric euclidean

# KL divergence
pyjaspar similarity MA0001.1 MA0002.1 --metric kl

# Find best alignment (Pearson only)
pyjaspar similarity MA0001.1 MA0002.1 --best

# With specific offset
pyjaspar similarity MA0001.1 MA0002.1 --offset 3

# JSON output
pyjaspar similarity MA0001.1 MA0002.1 --best --format json
```

### enrichment

Test motif enrichment in foreground vs background sequences.

```bash
pyjaspar enrichment \
    --foreground fg_sequences.fasta \
    --background bg_sequences.fasta \
    --motif-ids MA0001.1,MA0002.1,MA0139.1 \
    --threshold 0.8

# JSON output
pyjaspar enrichment \
    -fg fg.fasta -bg bg.fasta \
    --motif-ids MA0001.1 \
    --format json
```

## Deep Learning (DL) commands

### dl profile

Fetch a DL profile by ID.

```bash
pyjaspar dl profile DL0001.1
pyjaspar dl profile DL0001.1 --format json
```

### dl model

Fetch a DL model by ID.

```bash
pyjaspar dl model BP000001.1
pyjaspar dl model BP000001.1 --format json
```

### dl motif

Fetch a DL motif pattern by MO* ID.

```bash
pyjaspar dl motif MO000001.1
pyjaspar dl motif MO000001.1 --format json
```

### dl profiles

Search DL profiles.

```bash
pyjaspar dl profiles --tf-name REST
pyjaspar dl profiles --tax-group vertebrates --format json
```

### dl models

Search DL models.

```bash
pyjaspar dl models --tf-name REST
pyjaspar dl models --cell-line K562
pyjaspar dl models --tf-name REST --cell-line K562 --format json
```

## General commands

### releases

List available JASPAR releases.

```bash
pyjaspar releases
pyjaspar releases --latest
```

### collections

List JASPAR collections available in a release. The set of collections
varies by release (e.g. `CORE` and `UNVALIDATED` only from 2022 onward).

```bash
pyjaspar collections
pyjaspar collections -r 2020
```

### cite

Show citation information.

```bash
pyjaspar cite
```

## Output formats

Most commands support `--format json` or `--format tsv` (default).

- **TSV**: Tab-separated values, suitable for piping to `cut`, `awk`, etc.
- **JSON**: Structured output, suitable for `jq` or programmatic use.

### Scripting examples

```bash
# Get motif IDs as a list
pyjaspar metadata --tax-group Vertebrates | cut -f1 | tail -n +2

# Count vertebrate motifs
pyjaspar metadata --format json --tax-group Vertebrates | jq '.count'

# Pipe to downstream tools
pyjaspar motif-by-id MA0095.2 --motif-format pfm > YY1.pfm
```
