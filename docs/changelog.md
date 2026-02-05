# Changelog

## v5.0.0

Major release with new architecture, features, and JASPAR 2026 DL support.

### New features

- **Deep Learning (DL) collection** -- access BPNet-derived motif profiles, models, and patterns via `jdb.dl.*`
- **Analysis tools** -- motif scanning, similarity metrics, and enrichment analysis
- **CLI** -- full command-line interface with rich-click for motif retrieval, analysis, and DL queries
- **Docker support** -- pre-built Docker image at `ghcr.io/asntech/pyjaspar`

### Architecture changes

- Migrated to `src/` layout with `pyproject.toml` (hatchling)
- Modular subpackages: `analysis/`, `cli/`, `dl/`
- Custom exception hierarchy (`PyJasparError`, `DatabaseError`, etc.)
- Optional extras: `[cli]`, `[analysis]`, `[dl]`, `[all]`
- Full type annotations (Python 3.10+)
- Comprehensive test suite (170+ tests)

### Breaking changes

- Minimum Python version is now 3.10
- Package import path unchanged (`from pyjaspar import JasparDB`)
- Backward compatibility alias `jaspardb` still works

## Previous releases

See the [GitHub releases](https://github.com/asntech/pyjaspar/releases) page for older versions.
