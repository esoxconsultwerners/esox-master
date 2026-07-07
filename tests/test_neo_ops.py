"""Contract tests for AP_A6 (pipeline/ingest_neo_ops.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
NEO_OPS = ROOT / "data" / "interim" / "neo_ops.parquet"


@pytest.fixture(scope="session")
def df():
    if not NEO_OPS.exists():
        pytest.fail(f"{NEO_OPS} missing - run pipeline/ingest_neo_ops.py first")
    return pd.read_parquet(NEO_OPS)


def key(df):
    return df["number_mp"].astype("string").fillna(df["designation"])


def test_keys_unique(df):
    dups = key(df).duplicated().sum()
    assert dups == 0, f"{dups} duplicate keys after per-source dedup"


def test_nhats_row_count(df):
    # The NHATS list grows over time: 6,983 objects at first retrieval
    # (2026-07); the handover's [1500, 6000] window was already stale,
    # widened with headroom.
    n = int(df["nhats_min_dv_kms"].notna().sum())
    assert 1500 <= n <= 12000, f"NHATS row count {n:,}"


def test_neocc_row_count(df):
    n = int(df["neocc_risk_ps_max"].notna().sum())
    assert 800 <= n <= 3000, f"NEOCC risk row count {n:,}"


def test_nhats_dv_range(df):
    dv = df["nhats_min_dv_kms"].dropna()
    assert ((dv > 3.0) & (dv < 15.0)).all(), (
        f"min dv outside (3, 15) km/s: {dv[(dv <= 3) | (dv >= 15)].head().tolist()}"
    )


def test_impact_probabilities(df):
    for col in ("neocc_ip_max", "neocc_ip_cum", "sentry_ip"):
        ip = df[col].dropna()
        assert ((ip > 0) & (ip < 1)).all(), (
            f"{col} outside (0, 1): {ip[(ip <= 0) | (ip >= 1)].head().tolist()}"
        )


def test_palermo_scale(df):
    for col in ("neocc_risk_ps_max", "neocc_risk_ps_cum",
                "sentry_ps_max", "sentry_ps_cum"):
        ps = df[col].dropna()
        assert (ps < 2).all(), f"{col} >= 2: {ps[ps >= 2].head().tolist()}"


def test_apophis_absent_from_risk_layers(df):
    apo = df[(df["number_mp"] == 99942) | (df["designation"] == "2004 MN4")]
    if len(apo):
        row = apo.iloc[0]
        assert pd.isna(row["neocc_risk_ps_max"]), "Apophis on NEOCC risk list"
        assert pd.isna(row["sentry_ps_max"]), "Apophis on Sentry list"


@pytest.mark.parametrize("name", ["nhats", "neocc_risk", "sentry"])
def test_provenance(name):
    path = ROOT / "data" / "provenance" / f"{name}.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for k in ("url", "date", "retrieved_utc", "sha256", "row_count", "citation"):
        assert doc.get(k), f"{name} provenance missing {k}"
    assert "T" in doc["retrieved_utc"], "retrieved_utc is not a datetime"
    assert doc["row_count"] > 0
