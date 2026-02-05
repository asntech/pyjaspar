"""SQL query builder for JASPAR motif searches (private module)."""

from __future__ import annotations


def _ensure_list(value: str | int | list) -> list:
    """Normalize a single value or list to always be a list."""
    if isinstance(value, list):
        return value
    return [value]


def build_motif_query(
    *,
    collection: str | list[str] | None = None,
    tf_name: str | list[str] | None = None,
    tf_class: str | list[str] | None = None,
    tf_family: str | list[str] | None = None,
    species: str | int | list[str | int] | None = None,
    tax_group: str | list[str] | None = None,
    pazar_id: str | list[str] | None = None,
    medline: str | list[str] | None = None,
    data_type: str | list[str] | None = None,
) -> tuple[str, list[str]]:
    """Build a parameterized SQL query for fetching JASPAR motif internal IDs.

    All user-provided values use ``?`` placeholders to prevent SQL injection.

    Returns:
        A tuple of (sql_string, params) for use with ``cursor.execute()``.
    """
    tables: list[str] = ["MATRIX m"]
    where_clauses: list[str] = []
    params: list[str] = []

    # Select by MATRIX.COLLECTION
    if collection:
        values = [c.upper() for c in _ensure_list(collection)]
        placeholders = ",".join("?" * len(values))
        where_clauses.append(f"m.COLLECTION IN ({placeholders})")
        params.extend(values)

    # Select by MATRIX.NAME
    if tf_name:
        values = _ensure_list(tf_name)
        placeholders = ",".join("?" * len(values))
        where_clauses.append(f"m.NAME IN ({placeholders})")
        params.extend(values)

    # Select by MATRIX_SPECIES.TAX_ID
    if species:
        tables.append("MATRIX_SPECIES ms")
        where_clauses.append("m.ID = ms.ID")
        values = [str(s) for s in _ensure_list(species)]
        placeholders = ",".join("?" * len(values))
        where_clauses.append(f"ms.TAX_ID IN ({placeholders})")
        params.extend(values)

    # Annotation-based filters from MATRIX_ANNOTATION table.
    # Each filter follows the same pattern: join on annotation table,
    # match TAG, then match VAL(s).
    _annotation_filters: list[tuple[str | list[str] | None, str, str]] = [
        (tf_class, "ma1", "class"),
        (tf_family, "ma2", "family"),
        (pazar_id, "ma3", "pazar_tf_id"),
        (medline, "ma4", "medline"),
        (data_type, "ma5", "type"),
        (tax_group, "ma6", "tax_group"),
    ]

    for filter_value, alias, tag in _annotation_filters:
        if filter_value:
            tables.append(f"MATRIX_ANNOTATION {alias}")
            where_clauses.append(f"m.ID = {alias}.ID")
            where_clauses.append(f"{alias}.TAG = ?")
            params.append(tag)

            values = _ensure_list(filter_value)
            if tag == "tax_group":
                values = [str(v).lower() for v in values]
            placeholders = ",".join("?" * len(values))
            where_clauses.append(f"{alias}.VAL IN ({placeholders})")
            params.extend([str(v) for v in values])

    sql = f"SELECT DISTINCT(m.ID) FROM {', '.join(tables)}"
    if where_clauses:
        sql += f" WHERE {' AND '.join(where_clauses)}"
    sql += " COLLATE NOCASE"

    return sql, params
