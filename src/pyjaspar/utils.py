"""Utility functions for pyjaspar."""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path


def data_dir() -> Path:
    """Return the path to the directory containing JASPAR SQLite databases."""
    return Path(__file__).parent / "data"


def get_jaspardb_path(filename: str) -> Path:
    """Return the full path to a JASPAR SQLite database file.

    Args:
        filename: The SQLite database filename (e.g. 'JASPAR2026.sqlite').

    Returns:
        Absolute path to the SQLite file.

    Raises:
        FileNotFoundError: If the database file does not exist.
    """
    path = data_dir() / filename
    if not path.exists():
        raise FileNotFoundError(f"JASPAR database not found: {path}")
    return path


def dict_list_to_tsv(
    dict_list: list[dict[str, object]],
    keys: list[str],
    separator: str = "\t",
) -> str:
    """Convert a list of dictionaries to a TSV-formatted string.

    Args:
        dict_list: List of dictionaries to convert.
        keys: Column names (dict keys) to include, in order.
        separator: Field separator (default: tab).

    Returns:
        A string containing the TSV data with headers.
    """

    def stringify_value(value: object) -> str:
        if isinstance(value, list):
            return ";".join(str(v) for v in value)
        return str(value) if value is not None else ""

    processed = [{k: stringify_value(v) for k, v in d.items()} for d in dict_list]

    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=keys, delimiter=separator)
    writer.writeheader()
    for row in processed:
        writer.writerow(row)
    return buffer.getvalue()
