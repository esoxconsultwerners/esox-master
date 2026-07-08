"""Contract tests for C2 (analysis/c2_composition/c2_composition.py)."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "data" / "final" / "esox_master_core.parquet"
PROV = ROOT / "data" / "provenance" / "analog_composition.json"
PROPS = ROOT / "data" / "interim" / "meteorite_group_properties.csv"
DEGENERATE = {"CI", "CK", "CO", "R", "mesosiderite", "pallasite"}
# always emitted for every Gaia-spectrum object; oc_subgroup/properties are
# legitimately sparse (only for ordinary_chondrite / confident non-degenerate).
ALWAYS_COLS = ["analog_group_top", "analog_group_conf", "analog_distribution",
               "analog_status"]
ANALOG_COLS = ALWAYS_COLS + ["oc_subgroup"]


@pytest.fixture(scope="session")
def core():
    return pd.read_parquet(CORE)


@pytest.fixture(scope="session")
def prov():
    return json.loads(PROV.read_text())


def test_columns_exist_for_gaia_objects(core):
    spec = core[core["has_gaia_spectrum"]]
    for c in ANALOG_COLS:
        assert c in core.columns, f"{c} column missing"
    for c in ALWAYS_COLS:
        assert spec[c].notna().all(), f"{c} null for some spectrum object"


def test_top_values_valid(core, prov):
    valid = set(prov["separable_subset"]) | {"degenerate"}
    bad = set(core["analog_group_top"].dropna().unique()) - valid
    assert not bad, f"invalid analog_group_top values: {bad}"


def test_no_confident_oc_subgroup(core):
    assert set(core["oc_subgroup"].dropna().unique()) <= {"unresolved"}, \
        "oc_subgroup must only ever be 'unresolved' (never H/L/LL)"
    # ordinary_chondrite tops carry oc_subgroup=unresolved
    oc = core["analog_group_top"] == "ordinary_chondrite"
    assert (core.loc[oc, "oc_subgroup"] == "unresolved").all()


def test_conf_and_distribution(core):
    c = core["analog_group_conf"].dropna()
    assert ((c >= 0) & (c <= 1)).all(), "analog_group_conf outside [0,1]"
    for js in core["analog_distribution"].dropna().head(300):
        d = json.loads(js)
        if d:
            assert abs(sum(d.values()) - 1.0) < 0.02, f"distribution sums to {sum(d.values())}"


def test_degenerate_never_ok(core):
    assert (core.loc[core["analog_group_top"] == "degenerate", "analog_status"] != "ok").all()
    ok = core[core["analog_status"] == "ok"]
    assert not set(ok["analog_group_top"]) & DEGENERATE, "degenerate group marked ok"


def test_external_reproduces_c14_apophis():
    sys.path.insert(0, str(ROOT / "analysis" / "c2_composition"))
    import c2_composition as c2
    gts = pd.read_parquet(ROOT / "data" / "interim" / "groundtruth_spectra.parquet")
    apo = gts[gts["number_mp"] == 99942]
    r1 = c2.match_external_spectrum((apo["wavelength_um"] * 1000).to_numpy(),
                                    apo["reflectance"].to_numpy())
    assert r1["top"] == "ordinary_chondrite", f"Apophis top {r1['top']} (expected OC)"
    # determinism: identical result on a second call
    r2 = c2.match_external_spectrum((apo["wavelength_um"] * 1000).to_numpy(),
                                    apo["reflectance"].to_numpy())
    assert r1["distribution"] == r2["distribution"], "external match not deterministic"


def test_properties_citations_present():
    p = pd.read_csv(PROPS)
    assert len(p) > 5 and "citation" in p.columns
    assert p["citation"].astype(str).str.strip().ne("").all(), "blank citation"
    # no fabricated DOI: every row is either a real citation or an explicit TODO
    assert p["citation"].astype(str).str.contains("TODO").all(), \
        "citations must be explicit TODO placeholders (no fabricated refs)"


def test_model_derived_provenance(prov):
    assert "MODEL-DERIVED" in prov["kind"]
    assert "ordinary_chondrite" in prov["separable_subset"]
    assert set(prov["degenerate_groups"]) == DEGENERATE
