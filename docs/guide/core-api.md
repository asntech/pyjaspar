# Core API

The `JasparDB` class is the main entry point for querying JASPAR motifs.

## Connecting to a release

```python
from pyjaspar import JasparDB

# Latest release (JASPAR2026)
jdb = JasparDB()

# Specific release
jdb = JasparDB(release="JASPAR2024")

# Custom SQLite file
jdb = JasparDB(sqlite_db_path="/path/to/custom.sqlite")

# Context manager
with JasparDB() as jdb:
    motif = jdb.fetch_motif_by_id("MA0001.1")
```

## Fetch by ID

```python
# Full ID (base + version)
motif = jdb.fetch_motif_by_id("MA0095.2")

# Base ID only (returns latest version)
motif = jdb.fetch_motif_by_id("MA0095")

# Returns None if not found
motif = jdb.fetch_motif_by_id("MA9999.1")
assert motif is None
```

## Fetch by name

```python
motifs = jdb.fetch_motifs_by_name("CTCF")
for m in motifs:
    print(f"{m.matrix_id}: {m.name}")
```

## Advanced filtering

The `fetch_motifs()` method supports a wide range of filters:

```python
motifs = jdb.fetch_motifs(
    collection="CORE",        # JASPAR collection (default: "CORE")
    tf_name="CTCF",           # TF name(s)
    tf_class="Zinc finger",   # TF structural class(es)
    tf_family="C2H2",         # TF family(ies)
    tax_group="vertebrates",  # Taxonomic supergroup
    species=9606,             # Species taxonomy ID(s)
    data_type="ChIP-seq",     # Data type(s)
    min_ic=10,                # Minimum information content
    min_length=8,             # Minimum motif length
    min_sites=50,             # Minimum number of binding sites
    all_versions=True,        # Return all versions (not just latest)
)
```

All filter arguments accept either a single value or a list of values:

```python
# Multiple taxonomic groups
motifs = jdb.fetch_motifs(tax_group=["vertebrates", "insects"])

# Multiple species
motifs = jdb.fetch_motifs(species=[9606, 10090])  # human, mouse
```

### Filter parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `collection` | str | JASPAR collection (default: `"CORE"`, set to `None` for all) |
| `tf_name` | str or list | Transcription factor name(s) |
| `tf_class` | str or list | TF structural class(es) |
| `tf_family` | str or list | TF family(ies) |
| `matrix_id` | str or list | Specific JASPAR matrix ID(s) |
| `tax_group` | str or list | Taxonomic supergroup(s) |
| `species` | int or list | Species taxonomy ID(s) |
| `data_type` | str or list | Data type(s) used to build the matrix |
| `pazar_id` | str or list | PAZAR TF ID(s) |
| `medline` | str or list | PubMed ID(s) |
| `min_ic` | float | Minimum information content (specificity) |
| `min_length` | int | Minimum motif length |
| `min_sites` | int | Minimum number of binding sites |
| `all_versions` | bool | Return all versions of matching motifs |
| `all` | bool | Return every motif (ignores other filters) |

## Working with motif objects

All fetch methods return `Bio.motifs.jaspar.Motif` objects:

```python
motif = jdb.fetch_motif_by_id("MA0095.2")

# Identity
print(motif.matrix_id)   # MA0095.2
print(motif.base_id)     # MA0095
print(motif.name)        # YY1

# Annotations
print(motif.collection)  # CORE
print(motif.species)     # ['9606']
print(motif.species_name) # ['Homo sapiens']
print(motif.tf_class)    # ['C2H2 zinc finger factors']
print(motif.tf_family)   # ['More than 3 adjacent zinc fingers']
print(motif.tax_group)   # vertebrates
print(motif.data_type)   # ChIP-seq
print(motif.acc)         # UniProt accession(s)

# Counts matrix
print(motif.counts)
print(motif.length)      # number of positions

# Derived matrices
pwm = motif.counts.normalize(pseudocounts=0.5)
pssm = pwm.log_odds()
print(motif.consensus)

# Output formats
print(motif.format("jaspar"))
print(motif.format("pfm"))
print(motif.format("transfac"))
```

## Available releases

```python
releases = jdb.get_releases()
print(releases)
# ['JASPAR2026', 'JASPAR2024', 'JASPAR2022', 'JASPAR2020',
#  'JASPAR2018', 'JASPAR2016', 'JASPAR2014']
```

## Cross-release comparison

```python
jdb_2026 = JasparDB(release="JASPAR2026")
jdb_2020 = JasparDB(release="JASPAR2020")

motif_new = jdb_2026.fetch_motif_by_id("MA0001.1")
motif_old = jdb_2020.fetch_motif_by_id("MA0001.1")
```
