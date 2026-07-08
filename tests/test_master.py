"""Contract tests for Phase B (pipeline/build_master.py)."""

from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "data" / "final"
I = ROOT / "data" / "interim"
PROV = ROOT / "data" / "provenance"

FULL = FINAL / "esox_master_full.parquet"
CORE = FINAL / "esox_master_core.parquet"

BEST_SOURCE = [("diameter_best", "diameter_source"),
               ("albedo_best", "albedo_source"),
               ("period_best", "period_source"),
               ("h_best", "h_source")]


@pytest.fixture(scope="session")
def full():
    if not FULL.exists():
        pytest.fail(f"{FULL} missing - run pipeline/build_master.py first")
    return pd.read_parquet(FULL)


@pytest.fixture(scope="session")
def core():
    if not CORE.exists():
        pytest.fail(f"{CORE} missing - run pipeline/build_master.py first")
    return pd.read_parquet(CORE)


def test_object_key_unique(full, core):
    assert full["object_key"].is_unique, "object_key not unique in full"
    assert core["object_key"].is_unique, "object_key not unique in core"


def test_core_row_count(core):
    assert len(core) == 19190, f"core has {len(core)} rows, expected 19,190"


@pytest.mark.parametrize("best,src", BEST_SOURCE)
def test_best_source_consistency(full, best, src):
    bad = full[full[best].notna() & full[src].isna()]
    assert bad.empty, f"{len(bad)} rows with {best} non-null but {src} null"


def test_full_within_1pct_of_backbone(full):
    backbone = len(pd.read_parquet(I / "orbits_backbone.parquet", columns=["number_mp"]))
    assert abs(len(full) - backbone) / backbone <= 0.01, (
        f"full {len(full)} vs backbone {backbone} deviates > 1%"
    )


@pytest.mark.parametrize("num", [1, 4, 433, 101955, 162173, 99942])
def test_spot_objects_present(full, num):
    assert (full["object_key"] == str(num)).sum() == 1, f"object {num} not uniquely present"


def test_complete_profile_positive_and_deterministic(core):
    cp = int(core["complete_profile"].sum())
    assert cp > 0, "no complete physical profiles"
    # determinism/consistency: the stored flag must equal recomputation from
    # its component flags (a fixed function of the inputs, no randomness).
    recomputed = (core["has_gaia_spectrum"] & core["has_period_best"]
                  & core["has_family"] & core["has_diameter_best"]
                  & core["has_taxon_any"])
    assert (recomputed == core["complete_profile"]).all(), "complete_profile not deterministic"


def test_precedence_rules_updated():
    text = (PROV / "precedence_rules.md").read_text()
    assert "10.1093/pasj/psu037" in text, "Usui DOI missing"
    assert "10.1051/0004-6361/200912693" in text, "Durech DOI missing"
    assert "_pending_" not in text, "realized-outcomes placeholder still present"
    assert "| `diameter_best` | AKARI |" in text, "realized table not filled"
