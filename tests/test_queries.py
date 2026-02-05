"""Tests for the SQL query builder."""

from __future__ import annotations

from pyjaspar._queries import build_motif_query


def test_no_filters():
    sql, params = build_motif_query()
    assert "SELECT DISTINCT(m.ID) FROM MATRIX m" in sql
    assert params == []


def test_collection_single():
    sql, params = build_motif_query(collection="CORE")
    assert "m.COLLECTION IN (?)" in sql
    assert params == ["CORE"]


def test_collection_list():
    sql, params = build_motif_query(collection=["CORE", "UNVALIDATED"])
    assert "m.COLLECTION IN (?,?)" in sql
    assert params == ["CORE", "UNVALIDATED"]


def test_tf_name_single():
    sql, params = build_motif_query(tf_name="CTCF")
    assert "m.NAME IN (?)" in sql
    assert "CTCF" in params


def test_species():
    sql, params = build_motif_query(species=9606)
    assert "MATRIX_SPECIES ms" in sql
    assert "ms.TAX_ID IN (?)" in sql
    assert "9606" in params


def test_species_list():
    sql, params = build_motif_query(species=[9606, 10090])
    assert "ms.TAX_ID IN (?,?)" in sql
    assert "9606" in params
    assert "10090" in params


def test_tf_class():
    sql, params = build_motif_query(tf_class="Homeo domain factors")
    assert "MATRIX_ANNOTATION ma1" in sql
    assert "ma1.TAG = ?" in sql
    assert "class" in params
    assert "Homeo domain factors" in params


def test_tax_group_lowercased():
    sql, params = build_motif_query(tax_group="Vertebrates")
    assert "MATRIX_ANNOTATION ma6" in sql
    assert "ma6.TAG = ?" in sql
    assert "tax_group" in params
    assert "vertebrates" in params  # lowercased


def test_combined_filters():
    sql, params = build_motif_query(
        collection="CORE",
        tf_name="CTCF",
        tax_group="vertebrates",
    )
    assert "m.COLLECTION IN (?)" in sql
    assert "m.NAME IN (?)" in sql
    assert "MATRIX_ANNOTATION ma6" in sql
    assert "CORE" in params
    assert "CTCF" in params
    assert "vertebrates" in params


def test_no_string_interpolation():
    """Verify no SQL injection: params should use ? placeholders, not string formatting."""
    sql, params = build_motif_query(
        collection="'; DROP TABLE MATRIX; --",
        tf_name="'; DROP TABLE MATRIX; --",
    )
    # The SQL should contain ? placeholders, not the injected values
    assert "'; DROP TABLE MATRIX; --" not in sql
    assert "?" in sql
    # The dangerous values should be in params (safely passed to cursor.execute)
    assert "'; DROP TABLE MATRIX; --".upper() in params
