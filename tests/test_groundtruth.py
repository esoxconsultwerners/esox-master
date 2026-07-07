"""Contract tests for AP_A5 (pipeline/ingest_groundtruth_spectra.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "data" / "interim" / "groundtruth_summary.parquet"
SPECTRA = ROOT / "data" / "interim" / "groundtruth_spectra.parquet"


@pytest.fixture(scope="session")
def summary():
    if not SUMMARY.exists():
        pytest.fail(f"{SUMMARY} missing - run pipeline/ingest_groundtruth_spectra.py")
    return pd.read_parquet(SUMMARY)


@pytest.fixture(scope="session")
def spectra():
    if not SPECTRA.exists():
        pytest.fail(f"{SPECTRA} missing - run pipeline/ingest_groundtruth_spectra.py")
    return pd.read_parquet(SPECTRA)


def test_summary_keys_unique(summary):
    key = summary["number_mp"].astype("string").fillna(summary["designation"])
    dups = summary.assign(key=key).duplicated(subset=["key", "gt_source"]).sum()
    assert dups == 0, f"{dups} duplicate (object, gt_source) summary rows"


def test_ecas_row_count(summary):
    n = int((summary["gt_source"] == "ecas").sum())
    assert 550 <= n <= 600, f"ECAS row count {n}"


def test_smass2_object_count(summary):
    n = int((summary["gt_source"] == "smass2").sum())
    assert 1200 <= n <= 1600, f"SMASSII object count {n}"


def test_mithneos_object_count(summary):
    n = int((summary["gt_source"] == "mithneos").sum())
    assert n > 500, f"MITHNEOS object count {n}"


def test_wavelength_range(spectra):
    wl = spectra["wavelength_um"].dropna()
    assert ((wl > 0.3) & (wl < 2.6)).all(), (
        f"wavelengths outside (0.3, 2.6) um: "
        f"{wl[(wl <= 0.3) | (wl >= 2.6)].head().tolist()}"
    )


def test_reflectance_range(spectra):
    # Archive values are kept as published. A handful of MITHNEOS points
    # (visible anchor at 0.66 um) are non-physical noise, but every one of
    # them carries an archived reflectance_err > 1 - the archive itself
    # flags them. Strict (0, 5) applies to all points the archive trusts.
    ok = spectra.dropna(subset=["reflectance"])
    trusted = ok[ok["reflectance_err"].isna() | (ok["reflectance_err"] < 1)]
    r = trusted["reflectance"]
    assert ((r > 0) & (r < 5)).all(), (
        f"trusted reflectance outside (0, 5): {r[(r <= 0) | (r >= 5)].head().tolist()}"
    )
    outliers = ok[(ok["reflectance"] <= 0) | (ok["reflectance"] >= 5)]
    assert len(outliers) <= 20, f"{len(outliers)} out-of-range points (expected <= 20)"
    assert (outliers["reflectance_err"] > 1).all(), (
        "out-of-range reflectance without a large archived error"
    )


def test_min_points_per_spectrum(spectra):
    pts = spectra.groupby(["gt_source", "gt_spectrum_id"]).size()
    small = pts[pts < 5]
    assert small.empty, f"spectra with < 5 points: {small.head().to_dict()}"


def test_apophis_in_mithneos(summary, spectra):
    apo = summary[(summary["gt_source"] == "mithneos") & (summary["number_mp"] == 99942)]
    assert len(apo) == 1, "Apophis missing from MITHNEOS summary"
    pts = spectra[(spectra["gt_source"] == "mithneos") & (spectra["number_mp"] == 99942)]
    assert len(pts) >= 5, f"Apophis MITHNEOS spectrum has only {len(pts)} points"


@pytest.mark.parametrize("name", ["smass2", "mithneos", "ecas"])
def test_provenance(name):
    path = ROOT / "data" / "provenance" / f"{name}.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for key in ("url", "urn", "version", "doi", "date", "sha256", "citation"):
        assert doc.get(key), f"{name} provenance missing {key}"
    assert doc.get("row_count", 0) > 0
