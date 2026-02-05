"""Tests for DL matrix parsing."""

from __future__ import annotations

import numpy as np
import pytest

from pyjaspar.dl._parser import parse_matrix_json


def test_parse_pfm_matrix():
    """Parse a valid PFM JSON string."""
    json_str = (
        '{"A": "769 1511 681", "C": "2530 1715 648", "G": "1424 901 590", "T": "689 1285 3493"}'
    )
    result = parse_matrix_json(json_str, "PFM")
    assert result.shape == (4, 3)
    assert result[0, 0] == 769.0
    assert result[3, 2] == 3493.0


def test_parse_cwm_matrix_negative():
    """CWM matrices can have negative values."""
    json_str = '{"A": "-0.001 0.002", "C": "0.003 -0.004", "G": "0.005 0.006", "T": "-0.007 0.008"}'
    result = parse_matrix_json(json_str, "CWM")
    assert result.shape == (4, 2)
    assert result[0, 0] < 0
    assert result[0, 1] > 0


def test_parse_cwm_scientific_notation():
    """CWM values can use scientific notation."""
    json_str = '{"A": "3.84e-05", "C": "6.36e-05", "G": "1.65e-01", "T": "9.33e-07"}'
    result = parse_matrix_json(json_str, "CWM")
    assert result.shape == (4, 1)
    np.testing.assert_allclose(result[0, 0], 3.84e-05)


def test_parse_missing_nucleotide():
    """Missing nucleotide key should raise ValueError."""
    json_str = '{"A": "1 2", "C": "3 4", "G": "5 6"}'
    with pytest.raises(ValueError, match="Missing nucleotide 'T'"):
        parse_matrix_json(json_str, "PFM")


def test_parse_empty_values():
    """Empty value string should raise ValueError."""
    json_str = '{"A": "", "C": "1 2", "G": "3 4", "T": "5 6"}'
    with pytest.raises(ValueError, match="Empty values"):
        parse_matrix_json(json_str, "PFM")


def test_parse_inconsistent_lengths():
    """Inconsistent row lengths should raise ValueError."""
    json_str = '{"A": "1 2 3", "C": "4 5", "G": "6 7 8", "T": "9 10 11"}'
    with pytest.raises(ValueError, match="Inconsistent row lengths"):
        parse_matrix_json(json_str, "PFM")


def test_parse_malformed_json():
    """Malformed JSON should raise ValueError."""
    with pytest.raises(ValueError, match="Invalid"):
        parse_matrix_json("not json", "PFM")


def test_parse_non_numeric():
    """Non-numeric values should raise ValueError."""
    json_str = '{"A": "abc", "C": "1", "G": "2", "T": "3"}'
    with pytest.raises(ValueError, match="Non-numeric"):
        parse_matrix_json(json_str, "PFM")


def test_parse_result_dtype():
    """Result should be float64."""
    json_str = '{"A": "1", "C": "2", "G": "3", "T": "4"}'
    result = parse_matrix_json(json_str, "PFM")
    assert result.dtype == np.float64
