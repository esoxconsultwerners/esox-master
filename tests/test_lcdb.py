"""Contract tests for AP_A2 (pipeline/ingest_lcdb.py)."""

import json
import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
LCDB = ROOT / "data" / "interim" / "lcdb.parquet"
GASP_KEYS = ROOT / "data" / "interim" / "gasp_core_keys.parquet"

U_CODE = re.compile(r"^[123][+-]?$")


@pytest.fixture(scope="session")
def df():
    if not LCDB.exists():
        pytest.fail(f"{LCDB} missing - run pipeline/ingest_lcdb.py first")
    return pd.read_parquet(LCDB)


def test_number_mp_unique_among_non_null(df):
    numbered = df.loc[df["number_mp"].notna(), "number_mp"]
    assert numbered.is_unique, (
        f"{numbered.duplicated().sum()} duplicate number_mp values "
        "(dedup rule: keep highest lcdb_quality_numeric)"
    )


def test_row_count(df):
    assert len(df) > 15_000, f"only {len(df):,} rows"


def test_period_ranges(df):
    p = df["lcdb_period_h"].dropna()
    # LCDB contains real ultra-fast small-NEA rotators below 0.01 h
    # (14 objects, e.g. 2014 RC at 0.004389 h = 15.8 s; fastest entry
    # 0.0033 h = 11.9 s), so the sanity floor is 0.002 h (7.2 s).
    assert (p > 0.002).all(), f"periods <= 0.002 h: {p[p <= 0.002].head().tolist()}"
    assert (p < 20_000).all(), f"periods >= 20000 h: {p[p >= 20_000].head().tolist()}"


def test_quality_codes(df):
    u = df["lcdb_quality_u"].dropna()
    bad = u[~u.str.match(U_CODE)]
    assert bad.empty, f"invalid U codes: {bad.unique().tolist()[:10]}"


def test_gasp_core_match_rate(df):
    gasp = pd.read_parquet(GASP_KEYS)
    matched = gasp["number_mp"].isin(df.loc[df["number_mp"].notna(), "number_mp"]).sum()
    rate = matched / len(gasp)
    assert rate > 0.20, f"GASP core match rate {rate:.4f} ({matched}/{len(gasp)})"


def test_apophis_present_with_period(df):
    apo = df[df["number_mp"] == 99942]
    assert len(apo) == 1, "Apophis (99942) missing from LCDB table"
    assert pd.notna(apo.iloc[0]["lcdb_period_h"]), "Apophis has no period value"


def test_provenance():
    path = ROOT / "data" / "provenance" / "lcdb.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    assert doc.get("url"), "provenance missing url"
    assert doc.get("version"), "provenance missing PDS version"
    assert doc.get("doi"), "provenance missing DOI"
    assert doc.get("date"), "provenance missing date"
    assert doc.get("row_count", 0) > 0, "provenance missing row count"
