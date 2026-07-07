"""Contract tests for AP_A5b (pipeline/ingest_taxonomy.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
LONG = ROOT / "data" / "interim" / "taxonomy_long.parquet"
WIDE = ROOT / "data" / "interim" / "taxonomy_wide.parquet"


@pytest.fixture(scope="session")
def long_df():
    if not LONG.exists():
        pytest.fail(f"{LONG} missing - run pipeline/ingest_taxonomy.py first")
    return pd.read_parquet(LONG)


@pytest.fixture(scope="session")
def wide_df():
    if not WIDE.exists():
        pytest.fail(f"{WIDE} missing - run pipeline/ingest_taxonomy.py first")
    return pd.read_parquet(WIDE)


def test_long_rows_unique(long_df):
    key = long_df["number_mp"].astype("string").fillna(long_df["designation"])
    dups = long_df.assign(key=key).duplicated(
        subset=["key", "taxon_system", "taxon_class"]).sum()
    assert dups == 0, (
        f"{dups} duplicate (object, system, class) rows "
        "(source is one row per object; ingest keeps first occurrence)"
    )


def test_systems_present(long_df):
    n = long_df["taxon_system"].nunique()
    assert n >= 4, f"only {n} taxonomy systems present"


def test_distinct_objects(long_df):
    key = long_df["number_mp"].astype("string").fillna(long_df["designation"])
    n = key.nunique()
    assert n > 2_000, f"only {n} distinct objects"


def test_vesta_and_ceres(wide_df):
    vesta = wide_df[wide_df["number_mp"] == 4].iloc[0]
    ceres = wide_df[wide_df["number_mp"] == 1].iloc[0]
    assert vesta["taxon_tholen"] == "V", f"Vesta tholen {vesta['taxon_tholen']!r}"
    assert vesta["taxon_bus"] == "V", f"Vesta bus {vesta['taxon_bus']!r}"
    assert ceres["taxon_tholen"] == "G", f"Ceres tholen {ceres['taxon_tholen']!r}"
    assert ceres["taxon_bus"] == "C", f"Ceres bus {ceres['taxon_bus']!r}"


def test_class_strings_clean(long_df):
    cls = long_df["taxon_class"]
    assert cls.notna().all(), "null taxon_class in long table"
    assert (cls.str.len() > 0).all(), "empty taxon_class"
    assert (cls == cls.str.strip()).all(), "whitespace padding in taxon_class"


def test_provenance():
    path = ROOT / "data" / "provenance" / "taxonomy.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for k in ("url", "urn", "version", "doi", "date", "sha256",
              "row_count", "per_system_counts", "citation"):
        assert doc.get(k), f"provenance missing {k}"
    assert doc["row_count"] > 0
    assert len(doc["per_system_counts"]) >= 4
