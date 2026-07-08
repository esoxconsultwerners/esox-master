"""Contract tests for AP_A7b (pipeline/ingest_relab_returned.py).

Two documented deviations from the handover's literal test spec, both forced
by the genuine contents of this RELAB PDS4 bundle (verified by parsing all
23,606 cached labels and raw-grepping every label XML on 2026-07-08):

  1. Bodies: the handover expected "at least Itokawa and Ryugu". This bundle
     contains Ryugu only (31 spectra, sample C0002 / Hayabusa2); Itokawa and
     Bennu are absent (0 mentions anywhere). Per the project rule "report
     values as archived, never invent", we assert Ryugu present and report
     Itokawa/Bennu absence rather than fail - the same "report if absent, do
     not fail" stance the handover already granted Bennu.

  2. Wavelengths: the handover suggested (0.3, 3.0), but the archive spectra
     are kept pure (no truncation - a standing hard rule) and the FTIR-
     combined Ryugu products extend to ~167 um. We therefore assert a >= 0.3
     lower bound and require VIS-NIR coverage, but do not cap the IR tail.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
PROV = ROOT / "data" / "provenance"

SAMPLES = INTERIM / "relab_returned_samples.parquet"
SPECTRA = INTERIM / "relab_returned_spectra.parquet"


@pytest.fixture(scope="session")
def samples():
    if not SAMPLES.exists():
        pytest.fail(f"{SAMPLES} missing - run pipeline/ingest_relab_returned.py first")
    return pd.read_parquet(SAMPLES)


@pytest.fixture(scope="session")
def spectra():
    if not SPECTRA.exists():
        pytest.fail(f"{SPECTRA} missing - run pipeline/ingest_relab_returned.py first")
    return pd.read_parquet(SPECTRA)


def test_spectrum_ids_unique(samples):
    assert samples["relab_spectrum_id"].is_unique


def test_samples_spectra_consistent(samples, spectra):
    a = set(samples["relab_spectrum_id"])
    b = set(spectra["relab_spectrum_id"].unique())
    assert a == b, f"{len(a - b)} samples without points, {len(b - a)} orphan spectra"


def test_minimum_precious_count(samples):
    # Soft floor - a small precious set, deliberately not over-constrained.
    assert len(samples) >= 10, f"only {len(samples)} returned-sample spectra"


def test_ryugu_present_bodies_reported(samples):
    bodies = set(samples["parent_body"])
    assert "Ryugu" in bodies, f"Ryugu missing; bodies present: {bodies}"
    # Itokawa / Bennu are absent from this bundle - a legitimate finding,
    # not a failure. Guard only that classification did not silently fail.
    assert "unknown" not in bodies, f"unclassified parent_body present: {bodies}"


def test_reflectance_range(spectra):
    r = spectra["reflectance"]
    assert ((r > 0) & (r < 1.5)).all(), (
        f"reflectance outside (0, 1.5): {r[(r <= 0) | (r >= 1.5)].head().tolist()}"
    )


def test_wavelength_coverage(spectra, samples):
    wl = spectra["wavelength_um"]
    assert (wl >= 0.3).all(), f"wavelengths below 0.3 um: {wl[wl < 0.3].head().tolist()}"
    assert (wl < 200).all(), f"wavelengths >= 200 um: {wl[wl >= 200].head().tolist()}"
    # VIS-NIR (the C2 space-weathering window) must be covered by most spectra.
    visnir = samples[(samples["wl_min_um"] <= 0.8) & (samples["wl_max_um"] >= 2.4)]
    assert len(visnir) >= 10, f"only {len(visnir)} spectra span the VIS-NIR window"


def test_provenance_complete():
    path = PROV / "relab_returned.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for k in ("source", "urn", "version", "date", "row_count", "selection",
              "bodies_present", "citation", "license_terms"):
        assert doc.get(k) not in (None, "", [], {}), f"relab_returned.json missing {k}"
    assert "Ryugu" in doc["bodies_present"]
