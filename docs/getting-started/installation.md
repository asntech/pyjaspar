# Installation

pyJASPAR requires Python >= 3.10.

## pip

```bash
pip install pyjaspar
```

## conda / mamba

```bash
conda install -c bioconda pyjaspar
```

## Docker

```bash
docker pull ghcr.io/asntech/pyjaspar

# Run CLI commands directly
docker run ghcr.io/asntech/pyjaspar motif-by-id MA0095.2
docker run ghcr.io/asntech/pyjaspar releases --latest
```

## Development install

```bash
git clone https://github.com/asntech/pyjaspar.git
cd pyjaspar
pip install -e ".[dev]"
```

This installs all extras (CLI, analysis, DL, docs) plus development tools (pytest, ruff, mypy, pre-commit).

## Verify installation

```python
from pyjaspar import __version__, JasparDB

print(__version__)  # 5.0.0
jdb = JasparDB()
print(jdb)  # JASPAR release:JASPAR2026:...
```
