"""Matrix data parsing for the DL collection."""

from __future__ import annotations

import json

import numpy as np
from numpy.typing import NDArray

_ALPHABET = ("A", "C", "G", "T")


def parse_matrix_json(matrix_data: str, matrix_type: str) -> NDArray[np.float64]:
    """Parse a DL_MATRIX.MATRIX_DATA JSON string into a numpy array.

    The JSON format is::

        {"A": "v1 v2 v3 ...", "C": "...", "G": "...", "T": "..."}

    where values are whitespace-separated numbers.

    For PFM matrices, values are non-negative integers (stored as strings).
    For CWM matrices, values can be negative floats including scientific notation.

    Args:
        matrix_data: Raw JSON string from DL_MATRIX.MATRIX_DATA.
        matrix_type: ``"PFM"`` or ``"CWM"`` (used in error messages only).

    Returns:
        numpy array of shape ``(4, L)`` where rows are A, C, G, T.

    Raises:
        ValueError: If the JSON is malformed, nucleotides are missing,
            values are non-numeric, or row lengths are inconsistent.
    """
    try:
        data = json.loads(matrix_data)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid {matrix_type} matrix JSON: {e}") from e

    if not isinstance(data, dict):
        raise ValueError(
            f"Expected a JSON object for {matrix_type} matrix, got {type(data).__name__}"
        )

    rows: list[list[float]] = []
    for nt in _ALPHABET:
        if nt not in data:
            raise ValueError(
                f"Missing nucleotide '{nt}' in {matrix_type} matrix data. "
                f"Expected keys: {', '.join(_ALPHABET)}"
            )
        values_str = str(data[nt]).strip()
        if not values_str:
            raise ValueError(f"Empty values for nucleotide '{nt}' in {matrix_type} matrix")
        try:
            values = [float(v) for v in values_str.split()]
        except ValueError as e:
            raise ValueError(
                f"Non-numeric value for nucleotide '{nt}' in {matrix_type} matrix: {e}"
            ) from e
        rows.append(values)

    lengths = {len(r) for r in rows}
    if len(lengths) != 1:
        raise ValueError(
            f"Inconsistent row lengths in {matrix_type} matrix: "
            f"A={len(rows[0])}, C={len(rows[1])}, G={len(rows[2])}, T={len(rows[3])}"
        )

    return np.array(rows, dtype=np.float64)
