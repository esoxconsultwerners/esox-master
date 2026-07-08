"""Contract tests for C4 (analysis/c4_classifier/c4_classifier.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "data" / "final" / "esox_master_core.parquet"
PROV = ROOT / "data" / "provenance" / "taxon_esox.json"
COMPLEXES = {"S", "C", "X", "V", "D", "A", "K/L"}


@pytest.fixture(scope="session")
def core():
    return pd.read_parquet(CORE)


@pytest.fixture(scope="session")
def prov():
    return json.loads(PROV.read_text())


def test_taxon_esox_present_and_valid(core):
    assert "taxon_esox" in core.columns
    assert core["taxon_esox"].notna().all() and len(core) == 19190
    assert set(core["taxon_esox"].unique()) <= COMPLEXES | {"unclassified"}


def test_conf_range_and_unclassified_matches_threshold(core, prov):
    c = core["taxon_esox_conf"]
    assert ((c >= 0) & (c <= 1)).all(), "taxon_esox_conf outside [0,1]"
    thr = prov["confidence_threshold"]
    is_unc = core["taxon_esox"] == "unclassified"
    assert (core.loc[is_unc, "taxon_esox_conf"] < thr).all(), "unclassified above threshold"
    assert (core.loc[~is_unc, "taxon_esox_conf"] >= thr).all(), "classified below threshold"


def test_heldout_agreement(prov):
    assert prov["heldout_pds_agreement"] >= 0.80, prov["heldout_pds_agreement"]


def test_confident_coverage(core):
    cov = (core["taxon_esox"] != "unclassified").mean()
    assert 0.45 <= cov <= 0.65, f"confident coverage {cov:.3f} outside [0.45, 0.65]"


def test_completeness_payoff_and_deterministic(core):
    new_cp = int(core["complete_profile_esox"].sum())
    assert new_cp > 204, f"complete profiles {new_cp} not greater than 204"
    recomputed = (core["has_gaia_spectrum"] & core["has_period_best"] & core["has_family"]
                  & core["has_diameter_best"] & (core["taxon_esox"] != "unclassified"))
    assert int(recomputed.sum()) == new_cp, "complete_profile_esox not deterministic"


def test_model_derived_and_separate_from_literature(core, prov):
    assert "MODEL-DERIVED" in prov["kind"]
    # literature taxonomy columns preserved and distinct from the model column
    assert "taxon_mahlke" in core.columns and "taxon_esox" in core.columns
    both = core["taxon_mahlke"].notna()
    assert (core.loc[both, "taxon_mahlke"] != core.loc[both, "taxon_esox"]).any(), \
        "taxon_esox must not be a copy of the literature label"
