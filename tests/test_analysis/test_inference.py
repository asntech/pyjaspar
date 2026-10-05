"""Tests for profile inference (network calls are mocked)."""

from __future__ import annotations

import io
import json
import urllib.error
import urllib.request

import pytest

from pyjaspar.analysis.inference import INFER_HOSTS, InferenceHit, infer_profiles

RESPONSE = {
    "count": 2,
    "results": [
        {
            "matrix_id": "MA0162.2",
            "name": "EGR1",
            "evalue": 0.0,
            "dbd": 1.0,
            "url": "https://jaspar2024.elixir.no/api/v1/matrix/MA0162.2/",
            "sequence_logo": "https://jaspar2024.elixir.no/static/logos/svg/MA0162.2.svg",
        },
        {
            "matrix_id": "MA0732.1",
            "name": "EGR3",
            "evalue": 2.43e-90,
            "dbd": 0.884,
            "url": "https://jaspar2024.elixir.no/api/v1/matrix/MA0732.1/",
            "sequence_logo": "https://jaspar2024.elixir.no/static/logos/svg/MA0732.1.svg",
        },
    ],
}
EMPTY = {"count": 0, "results": None}


@pytest.fixture
def fake_api(monkeypatch):
    """Replace urlopen; ``fake_api.payload`` is returned, ``fake_api.urls`` records requests."""

    class Api:
        payload = RESPONSE
        urls: list[str] = []

    def urlopen(url, timeout=None):
        Api.urls.append(url)
        return io.BytesIO(json.dumps(Api.payload).encode())

    Api.urls = []
    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    return Api


def test_infer_profiles_parses_hits(fake_api):
    hits = infer_profiles("MKLAA")
    assert hits == [
        InferenceHit(
            "MA0162.2",
            "EGR1",
            0.0,
            1.0,
            "https://jaspar2024.elixir.no/static/logos/svg/MA0162.2.svg",
            "JASPAR2024",
        ),
        InferenceHit(
            "MA0732.1",
            "EGR3",
            2.43e-90,
            0.884,
            "https://jaspar2024.elixir.no/static/logos/svg/MA0732.1.svg",
            "JASPAR2024",
        ),
    ]


def test_infer_profiles_empty_result(fake_api):
    fake_api.payload = EMPTY
    assert infer_profiles("MKLAA") == []


def test_infer_profiles_uses_release_host_and_cleans_sequence(fake_api):
    infer_profiles(" mkl\naa* ", release="JASPAR2022")
    assert fake_api.urls == [f"{INFER_HOSTS['JASPAR2022']}/api/v1/infer/MKLAA/?format=json"]


def test_infer_profiles_sets_release_on_hits(fake_api):
    hits = infer_profiles("MKLAA", release="JASPAR2020")
    assert {h.release for h in hits} == {"JASPAR2020"}


def test_infer_profiles_release_without_service(fake_api):
    with pytest.raises(ValueError, match="not available"):
        infer_profiles("MKLAA", release="JASPAR2026")
    assert fake_api.urls == []


@pytest.mark.parametrize("sequence", ["", "MK1AA", "MK-AA"])
def test_infer_profiles_invalid_sequence(fake_api, sequence):
    with pytest.raises(ValueError, match="amino acid"):
        infer_profiles(sequence)
    assert fake_api.urls == []
