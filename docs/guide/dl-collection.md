# Deep Learning (DL) Collection

JASPAR 2026 introduces a Deep Learning collection derived from BPNet models trained on
ENCODE TF ChIP-seq data.

## Entity types

The DL collection has three entity types:

| Entity | ID prefix | Description |
|--------|-----------|-------------|
| Profile | DL* | TF summary profile with primary and alternate motif patterns |
| Motif pattern | MO* | Cluster-level motif with PFM and CWM matrices |
| Model | BP* | BPNet model trained on a specific TF and cell line |

## Python API

### Fetch a profile

```python
from pyjaspar import JasparDB

jdb = JasparDB()  # JASPAR2026

profile = jdb.dl.fetch_profile("DL0001.1")
print(profile.tf_name)         # ARNT2
print(profile.tax_group)       # vertebrates
print(profile.primary_motif)   # primary motif pattern
print(profile.alt_motifs)      # alternate motif patterns
print(profile.jaspar_matches)  # best classic JASPAR matches
print(profile.linked_model_ids)  # linked BPNet model IDs
```

### Fetch a model

```python
model = jdb.dl.fetch_model("BP000001.1")
print(model.tf_name)       # REST
print(model.cell_line)     # K562
print(model.data_type)     # ChIP-seq
print(len(model.motifs))   # number of motif patterns
```

### Fetch a motif pattern

```python
motif = jdb.dl.fetch_motif_pattern("MO000001.1")
print(motif.motif_id)       # MO000001.1
print(motif.matrices.keys())  # dict_keys(['PFM', 'CWM'])

# Access PFM and CWM matrices
pfm = motif.pfm   # Matrix object with .values (numpy array) and .length
cwm = motif.cwm   # Contribution Weight Matrix

# Convert PFM to BioPython Motif object
bio_motif = motif.to_biopython()
print(bio_motif.pssm)
```

### Search profiles and models

```python
# Search profiles by TF name
profiles = jdb.dl.search_profiles(tf_name="REST")
for p in profiles:
    print(f"{p.profile_id}: {p.tf_name} ({p.tax_group})")

# Search models by cell line
models = jdb.dl.search_models(cell_line="K562")
for m in models:
    print(f"{m.model_id}: {m.tf_name} ({m.cell_line})")

# Combined filters
models = jdb.dl.search_models(tf_name="REST", cell_line="K562")
```

### Search parameters

**`search_profiles()`**

| Parameter | Type | Description |
|-----------|------|-------------|
| `tf_name` | str | TF name filter |
| `tax_group` | str | Taxonomic group filter |

**`search_models()`**

| Parameter | Type | Description |
|-----------|------|-------------|
| `tf_name` | str | TF name filter |
| `cell_line` | str | Cell line filter |
| `tax_id` | int | Taxonomy ID filter |
| `data_type` | str | Data type filter |
| `model_name` | str | Model name filter |
| `source` | str | Source filter |

## Availability

The `.dl` property is only available on JASPAR2026+. Accessing it on older releases
raises `DLNotAvailableError`:

```python
from pyjaspar import JasparDB, DLNotAvailableError

try:
    jdb = JasparDB(release="JASPAR2020")
    jdb.dl  # raises DLNotAvailableError
except DLNotAvailableError as e:
    print(e)
```
