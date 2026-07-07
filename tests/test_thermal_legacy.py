"""Contract tests for AP_A4 (pipeline/ingest_thermal_legacy.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
AKARI = ROOT / "data" / "interim" / "akari.parquet"
IRAS = ROOT / "data" / "interim" / "iras_simps.parquet"


@pytest.fixture(scope="session")
def akari():
    if not AKARI.exists():
        pytest.fail(f"{AKARI} missing - run pipeline/ingest_thermal_legacy.py first")
    return pd.read_parquet(AKARI)


@pytest.fixture(scope="session")
def iras():
    if not IRAS.exists():
        pytest.fail(f"{IRAS} missing - run pipeline/ingest_thermal_legacy.py first")
    return pd.read_parquet(IRAS)


def test_akari_number_unique(akari):
    numbered = akari.loc[akari["number_mp"].notna(), "number_mp"]
    assert numbered.is_unique, "duplicate number_mp in AKARI"


def test_iras_number_unique(iras):
    numbered = iras.loc[iras["number_mp"].notna(), "number_mp"]
    assert numbered.is_unique, "duplicate number_mp in IRAS/SIMPS"


def test_akari_row_count(akari):
    assert 4500 <= len(akari) <= 6000, f"AKARI row count {len(akari):,}"


def test_iras_row_count(iras):
    assert 2000 <= len(iras) <= 2600, f"IRAS row count {len(iras):,}"


@pytest.mark.parametrize("src,alb,diam", [
    ("akari", "akari_albedo", "akari_diameter_km"),
    ("iras", "iras_albedo", "iras_diameter_km"),
])
def test_value_ranges(akari, iras, src, alb, diam):
    df = akari if src == "akari" else iras
    a = df[alb].dropna()
    d = df[diam].dropna()
    assert ((a > 0) & (a <= 1.1)).all(), (
        f"{src} albedo out of (0, 1.1]: {a[(a <= 0) | (a > 1.1)].head().tolist()}"
    )
    assert ((d > 0.1) & (d < 1100)).all(), (
        f"{src} diameter out of (0.1, 1100): {d[(d <= 0.1) | (d >= 1100)].head().tolist()}"
    )


def test_ceres_present_in_both(akari, iras):
    assert (akari["number_mp"] == 1).any(), "Ceres missing from AKARI"
    assert (iras["number_mp"] == 1).any(), "Ceres missing from IRAS/SIMPS"


def test_apophis_absent_in_both(akari, iras):
    assert not (akari["number_mp"] == 99942).any(), "Apophis unexpectedly in AKARI"
    assert not (iras["number_mp"] == 99942).any(), "Apophis unexpectedly in IRAS/SIMPS"


def test_err_columns_non_negative(akari, iras):
    for df, cols in ((akari, ["akari_diameter_err", "akari_albedo_err"]),
                     (iras, ["iras_diameter_err", "iras_albedo_err"])):
        for c in cols:
            v = df[c].dropna()
            assert (v >= 0).all(), f"negative values in {c}: {v[v < 0].head().tolist()}"


@pytest.mark.parametrize("name,extra_keys", [
    ("akari", ["version", "bibcode", "source_used"]),
    ("iras_simps", ["version", "doi", "dataset_id"]),
])
def test_provenance(name, extra_keys):
    path = ROOT / "data" / "provenance" / f"{name}.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for key in ["url", "date", "sha256", "row_count", "citation"] + extra_keys:
        assert doc.get(key), f"{name} provenance missing {key}"
    assert doc["row_count"] > 0
