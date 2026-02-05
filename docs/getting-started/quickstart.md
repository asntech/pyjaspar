# Quick Start

## Python API

```python
from pyjaspar import JasparDB

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

All methods return `Bio.motifs.jaspar.Motif` objects from BioPython.

### Working with motifs

```python
motif = jdb.fetch_motif_by_id("MA0095.2")

# Access counts matrix
print(motif.counts)

# Get Position Weight Matrix
pwm = motif.counts.normalize(pseudocounts=0.5)

# Get PSSM (log-odds)
pssm = pwm.log_odds()

# Consensus sequence
print(motif.consensus)

# Format output
print(motif.format("jaspar"))
print(motif.format("transfac"))
```

## CLI

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
```

## Next steps

- [Core API](../guide/core-api.md) -- detailed API guide
- [Analysis Tools](../guide/analysis.md) -- scanning, similarity, enrichment
- [DL Collection](../guide/dl-collection.md) -- JASPAR 2026 deep learning data
- [CLI Reference](../guide/cli.md) -- all CLI commands
