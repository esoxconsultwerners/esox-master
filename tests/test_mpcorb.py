"""Contract tests for AP_A1 (pipeline/ingest_mpcorb.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BACKBONE = ROOT / "data" / "interim" / "orbits_backbone.parquet"
GASP_KEYS = ROOT / "data" / "interim" / "gasp_core_keys.parquet"


@pytest.fixture(scope="session")
def df():
    if not BACKBONE.exists():
        pytest.fail(f"{BACKBONE} missing - run pipeline/ingest_mpcorb.py first")
    return pd.read_parquet(BACKBONE)


def test_number_mp_unique_among_non_null(df):
    numbered = df.loc[df["number_mp"].notna(), "number_mp"]
    assert numbered.is_unique, "duplicate number_mp values found"


def test_row_count(df):
    assert len(df) > 1_300_000, f"only {len(df):,} rows"


def test_gasp_core_match_rate(df):
    gasp = pd.read_parquet(GASP_KEYS)
    matched = gasp["number_mp"].isin(df.loc[df["number_mp"].notna(), "number_mp"]).sum()
    rate = matched / len(gasp)
    assert rate > 0.995, f"GASP core match rate {rate:.4f} ({matched}/{len(gasp)})"


def test_orbital_element_ranges(df):
    a = df["mpc_a"].dropna()
    e = df["mpc_e"].dropna()
    i = df["mpc_i"].dropna()
    assert (a > 0).all(), f"non-positive mpc_a: {a[a <= 0].head().tolist()}"
    # MPCORB contains a few real extreme scattered-disk / inner-Oort objects
    # beyond 2000 au (2014 FE72, 2015 FF504, 2015 FW539); allow that tail
    # but keep the window strict for the rest of the catalog.
    n_extreme = int((a >= 2000).sum())
    assert n_extreme <= 10, (
        f"{n_extreme} objects with mpc_a >= 2000 au: "
        f"{a[a >= 2000].head().tolist()}"
    )
    assert ((e >= 0) & (e < 1.2)).all(), (
        f"mpc_e out of [0, 1.2): {e[(e < 0) | (e >= 1.2)].head().tolist()}"
    )
    assert ((i >= 0) & (i <= 180)).all(), (
        f"mpc_i out of [0, 180]: {i[(i < 0) | (i > 180)].head().tolist()}"
    )


def test_neo_count_sanity(df):
    n_neo = int(df["mpc_neo_flag"].sum())
    assert 30_000 < n_neo < 60_000, f"NEO count {n_neo:,} outside 2026 sanity window"


@pytest.mark.parametrize("name", ["mpcorb", "sbdb"])
def test_provenance(name):
    path = ROOT / "data" / "provenance" / f"{name}.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    assert doc.get("url"), "provenance missing url"
    assert doc.get("date"), "provenance missing date"
    assert doc.get("row_count", 0) > 0, "provenance missing row count"
