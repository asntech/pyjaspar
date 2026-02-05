"""SQL queries for the JASPAR DL collection."""

from __future__ import annotations


def _ensure_list(value: str | int | list[str] | list[int]) -> list[str | int]:
    """Normalize a single value or list to always be a list."""
    if isinstance(value, list):
        return value
    return [value]


# --- Profile queries ---

FETCH_PROFILE_BY_ID = (
    "SELECT ID, BASE_ID, VERSION, PRIMARY_MOTIF_ID, TF_NAME, TAX_GROUP "
    "FROM DL_PROFILE_SUMMARY "
    "WHERE BASE_ID = ? AND VERSION = ?"
)

FETCH_PROFILE_LATEST_VERSION = (
    "SELECT VERSION FROM DL_PROFILE_SUMMARY WHERE BASE_ID = ? ORDER BY VERSION DESC LIMIT 1"
)

FETCH_PROFILE_CLUSTERS = (
    "SELECT CLUSTER_ID, CLUSTER_VERSION "
    "FROM DL_PROFILE_SUMMARY_CLUSTER "
    "WHERE BASE_ID = ? AND VERSION = ?"
)

# --- Cluster queries ---

FETCH_CLUSTER = (
    "SELECT ID, BASE_ID, VERSION, MOTIF_ID, TF_NAME "
    "FROM DL_CLUSTER "
    "WHERE BASE_ID = ? AND VERSION = ?"
)

FETCH_CLUSTER_LATEST_VERSION = (
    "SELECT VERSION FROM DL_CLUSTER WHERE BASE_ID = ? ORDER BY VERSION DESC LIMIT 1"
)

FETCH_CLUSTER_JASPAR_MATCH = (
    "SELECT JASPAR_MATCH FROM DL_CLUSTER_JASPAR_MATCH WHERE BASE_ID = ? AND VERSION = ?"
)

FETCH_CLUSTER_MODELS = (
    "SELECT MODEL_ID, MODEL_VERSION "
    "FROM DL_CLUSTER_MODEL "
    "WHERE CLUSTER_ID = ? AND CLUSTER_VERSION = ?"
)

# --- Model queries ---

FETCH_MODEL_BY_ID = (
    "SELECT ID, BASE_ID, VERSION, PRIMARY_MOTIF_ID, TF_NAME, CELL_LINE, "
    "TAX_ID, DATA_TYPE, MODEL_NAME, SOURCE, SOURCE_ID, SOURCE_URL "
    "FROM DL_MODEL "
    "WHERE BASE_ID = ? AND VERSION = ?"
)

FETCH_MODEL_LATEST_VERSION = (
    "SELECT VERSION FROM DL_MODEL WHERE BASE_ID = ? ORDER BY VERSION DESC LIMIT 1"
)

FETCH_MODEL_MOTIFS = (
    "SELECT ID, NAME, SOURCE, SOURCE_ID, SOURCE_VERSION, "
    "NUM_SEQLETS, NUM_HITS, MOTIF_TYPE "
    "FROM DL_MOTIF "
    "WHERE SOURCE_ID = ? AND SOURCE_VERSION = ? "
    "ORDER BY ID"
)

# --- Motif queries ---

FETCH_MOTIF_BY_ID = (
    "SELECT ID, NAME, SOURCE, SOURCE_ID, SOURCE_VERSION, "
    "NUM_SEQLETS, NUM_HITS, MOTIF_TYPE "
    "FROM DL_MOTIF WHERE ID = ?"
)

FETCH_MOTIF_MATRICES = "SELECT MATRIX_TYPE, MATRIX_DATA FROM DL_MATRIX WHERE MOTIF_ID = ?"


# --- Search query builders ---


def build_profile_search_query(
    *,
    tf_name: str | list[str] | None = None,
    tax_group: str | list[str] | None = None,
) -> tuple[str, list[str]]:
    """Build a parameterized search query for DL profiles.

    Returns:
        A tuple of ``(sql_string, params)``.
    """
    where_clauses: list[str] = []
    params: list[str] = []

    if tf_name is not None:
        values = _ensure_list(tf_name)
        placeholders = ",".join("?" * len(values))
        where_clauses.append(f"TF_NAME IN ({placeholders})")
        params.extend(str(v) for v in values)

    if tax_group is not None:
        values = [str(v).lower() for v in _ensure_list(tax_group)]
        placeholders = ",".join("?" * len(values))
        where_clauses.append(f"LOWER(TAX_GROUP) IN ({placeholders})")
        params.extend(values)

    sql = (
        "SELECT ID, BASE_ID, VERSION, PRIMARY_MOTIF_ID, TF_NAME, TAX_GROUP FROM DL_PROFILE_SUMMARY"
    )
    if where_clauses:
        sql += " WHERE " + " AND ".join(where_clauses)
    sql += " ORDER BY BASE_ID, VERSION"

    return sql, params


def build_model_search_query(
    *,
    tf_name: str | list[str] | None = None,
    cell_line: str | list[str] | None = None,
    tax_id: int | list[int] | None = None,
    data_type: str | list[str] | None = None,
    model_name: str | list[str] | None = None,
    source: str | list[str] | None = None,
) -> tuple[str, list[str]]:
    """Build a parameterized search query for DL models.

    Returns:
        A tuple of ``(sql_string, params)``.
    """
    where_clauses: list[str] = []
    params: list[str] = []

    filters: list[tuple[str | None, str]] = [
        (tf_name, "TF_NAME"),
        (cell_line, "CELL_LINE"),
        (data_type, "DATA_TYPE"),
        (model_name, "MODEL_NAME"),
        (source, "SOURCE"),
    ]

    for value, column in filters:
        if value is not None:
            values = _ensure_list(value)
            placeholders = ",".join("?" * len(values))
            where_clauses.append(f"{column} IN ({placeholders})")
            params.extend(str(v) for v in values)

    if tax_id is not None:
        values = _ensure_list(tax_id)
        placeholders = ",".join("?" * len(values))
        where_clauses.append(f"TAX_ID IN ({placeholders})")
        params.extend(str(v) for v in values)

    sql = (
        "SELECT ID, BASE_ID, VERSION, PRIMARY_MOTIF_ID, TF_NAME, CELL_LINE, "
        "TAX_ID, DATA_TYPE, MODEL_NAME, SOURCE, SOURCE_ID, SOURCE_URL "
        "FROM DL_MODEL"
    )
    if where_clauses:
        sql += " WHERE " + " AND ".join(where_clauses)
    sql += " ORDER BY BASE_ID, VERSION"

    return sql, params
