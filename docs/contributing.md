# Contributing

## Development setup

```bash
git clone https://github.com/asntech/pyjaspar.git
cd pyjaspar
pip install -e ".[dev]"
```

This installs all dependencies including testing, linting, and documentation tools.

## Running tests

```bash
# Full test suite
pytest tests/ -v

# With coverage
pytest --cov=pyjaspar --cov-report=html tests/

# Specific test file
pytest tests/test_api.py -v
```

## Linting and formatting

The project uses [ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
# Check for lint errors
ruff check src/ tests/

# Auto-fix lint errors
ruff check --fix src/ tests/

# Check formatting
ruff format --check src/ tests/

# Apply formatting
ruff format src/ tests/
```

## Pre-commit hooks

Pre-commit hooks are configured in `.pre-commit-config.yaml`:

```bash
pre-commit install
pre-commit run --all-files
```

## Type checking

```bash
mypy src/pyjaspar/
```

## Documentation

```bash
# Install docs dependencies
pip install -e ".[docs]"

# Serve locally (with live reload)
mkdocs serve

# Build static site
mkdocs build --strict
```

## Project structure

```
src/pyjaspar/
├── __init__.py          # Package exports
├── _constants.py        # Version, releases, defaults
├── _exceptions.py       # Custom exceptions
├── _queries.py          # SQL query builders
├── db.py                # JasparDB main class
├── utils.py             # Utility functions
├── analysis/            # Analysis tools (scanning, similarity, enrichment)
├── cli/                 # CLI commands (main, dl, analysis)
├── dl/                  # Deep Learning collection support
└── data/                # Bundled JASPAR SQLite databases
```

## Submitting changes

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes
4. Run tests and linting: `pytest && ruff check src/ tests/`
5. Commit and push
6. Open a pull request
