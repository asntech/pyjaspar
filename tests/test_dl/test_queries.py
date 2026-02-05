"""Tests for DL query builders."""

from __future__ import annotations

from pyjaspar.dl._queries import build_model_search_query, build_profile_search_query


def test_profile_search_no_filters():
    sql, params = build_profile_search_query()
    assert "WHERE" not in sql
    assert params == []
    assert "ORDER BY" in sql


def test_profile_search_tf_name():
    sql, params = build_profile_search_query(tf_name="REST")
    assert "TF_NAME IN (?)" in sql
    assert params == ["REST"]


def test_profile_search_tax_group():
    sql, params = build_profile_search_query(tax_group="Vertebrates")
    assert "LOWER(TAX_GROUP) IN (?)" in sql
    assert params == ["vertebrates"]


def test_profile_search_combined():
    sql, params = build_profile_search_query(tf_name="REST", tax_group="vertebrates")
    assert "TF_NAME IN (?)" in sql
    assert "LOWER(TAX_GROUP) IN (?)" in sql
    assert len(params) == 2


def test_profile_search_list_values():
    sql, params = build_profile_search_query(tf_name=["REST", "CTCF"])
    assert "TF_NAME IN (?,?)" in sql
    assert params == ["REST", "CTCF"]


def test_model_search_no_filters():
    sql, params = build_model_search_query()
    assert "WHERE" not in sql
    assert params == []


def test_model_search_tf_name():
    sql, params = build_model_search_query(tf_name="REST")
    assert "TF_NAME IN (?)" in sql
    assert params == ["REST"]


def test_model_search_cell_line():
    sql, params = build_model_search_query(cell_line="K562")
    assert "CELL_LINE IN (?)" in sql
    assert params == ["K562"]


def test_model_search_tax_id():
    sql, params = build_model_search_query(tax_id=9606)
    assert "TAX_ID IN (?)" in sql
    assert params == ["9606"]


def test_model_search_combined():
    sql, params = build_model_search_query(tf_name="REST", cell_line="K562")
    assert "TF_NAME IN (?)" in sql
    assert "CELL_LINE IN (?)" in sql
    assert len(params) == 2


def test_no_string_interpolation():
    """Ensure no raw string interpolation in queries."""
    sql, params = build_profile_search_query(tf_name="'; DROP TABLE --")
    assert "'; DROP TABLE --" not in sql
    assert "'; DROP TABLE --" in params
