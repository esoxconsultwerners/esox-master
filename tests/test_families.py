"""Contract tests for AP_A3 (pipeline/ingest_families.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ROOT / "data" / "interim" / "families.parquet"
GASP_KEYS = ROOT / "data" / "interim" / "gasp_core_keys.parquet"


@pytest.fixture(scope="session")
def df():
    if not FAMILIES.exists():
        pytest.fail(f"{FAMILIES} missing - run pipeline/ingest_families.py first")
    return pd.read_parquet(FAMILIES)


def test_membership_key_unique(df):
    numbered = df[df["number_mp"].notna()]
    dups = numbered.duplicated(subset=["number_mp", "fam_id"]).sum()
    assert dups == 0, (
        f"{dups} duplicate (number_mp, fam_id) pairs "
        "(ingest drops exact duplicates; investigate source overlap)"
    )


def test_total_members(df):
    assert len(df) > 100_000, f"only {len(df):,} members"


def test_distinct_families(df):
    n = df["fam_id"].nunique()
    assert n > 100, f"only {n} distinct families"


def test_fam_name_non_null(df):
    with_id = df[df["fam_id"].notna()]
    n_null = int(with_id["fam_name"].isna().sum())
    assert n_null == 0, f"{n_null} rows with fam_id but null fam_name"


def test_gasp_core_match_rate(df):
    gasp = pd.read_parquet(GASP_KEYS)
    members = df.loc[df["number_mp"].notna(), "number_mp"].unique()
    matched = gasp["number_mp"].isin(members).sum()
    rate = matched / len(gasp)
    assert rate > 0.30, f"GASP core match rate {rate:.4f} ({matched}/{len(gasp)})"


def test_vesta_family(df):
    vesta = df[df["fam_name"].str.contains("Vesta", case=False, na=False)]
    assert len(vesta) > 10_000, f"Vesta family has only {len(vesta):,} members"


def test_apophis_no_family(df):
    apo = df[df["number_mp"] == 99942]
    assert len(apo) == 0, (
        f"Apophis unexpectedly assigned to family {apo['fam_name'].tolist()} - "
        "NEOs should not appear in main-belt family lists"
    )


def test_provenance():
    path = ROOT / "data" / "provenance" / "families.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    assert doc.get("urn"), "provenance missing urn"
    assert doc.get("version"), "provenance missing version"
    assert doc.get("doi"), "provenance missing DOI"
    assert doc.get("date"), "provenance missing date"
    assert doc.get("row_count", 0) > 0, "provenance missing row count"
