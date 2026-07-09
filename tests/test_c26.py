"""Contract tests for C2.6 - external NIR path migrated onto the C2.5 matcher.

The tests encode what Step 1 MEASURED (nir/c26_numbers.json), not a hoped-for
result: the NIR ground-truth path is domain-valid (mostly in-manifold), Apophis
is in-manifold, and confident status still requires clearing the manifold gate.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis" / "c1_degeneration"))
sys.path.insert(0, str(ROOT / "analysis" / "c2_composition"))
NUM = ROOT / "analysis" / "c2_composition" / "nir" / "c26_numbers.json"
GTS = ROOT / "data" / "interim" / "groundtruth_spectra.parquet"


@pytest.fixture(scope="session")
def nums():
    return json.loads(NUM.read_text())


@pytest.fixture(scope="session")
def apo():
    import c2_composition as c2
    g = pd.read_parquet(GTS)
    a = g[g["number_mp"] == 99942]
    return c2.match_external_spectrum((a["wavelength_um"] * 1000).to_numpy(),
                                      a["reflectance"].to_numpy())


def test_nir_manifold_recorded(nums):
    assert 0 <= nums["nir_inside_frac_pct"] <= 100
    assert nums["apophis_inside"] == (nums["apophis_d2"] <= nums["nir_selfq95"])
    # NIR ground-truth is domain-valid, unlike the ~5%-inside Gaia core
    assert nums["nir_inside_frac_pct"] > 50


def test_external_uses_c25_matcher(apo):
    assert apo["matcher"] == "c2.5-per-group", "RF must be retired for composition"


def test_confident_requires_manifold(apo):
    # confident status can only be emitted inside the manifold - same gate as the core
    if apo["analog_status"] == "ok":
        assert apo["within_manifold"]
    if apo["confident_analog"] not in ("unresolved", "degenerate"):
        assert apo["analog_status"] == "ok" and apo["within_manifold"]


def test_apophis_matches_measured_verdict(apo, nums):
    assert abs(apo["manifold_d2"] - nums["apophis_d2"]) < 0.05
    assert apo["within_manifold"] == nums["apophis_inside"]
    assert apo["top"] == "ordinary_chondrite"
    # measured: Apophis is in-manifold but NOT confidently single-group -> unresolved
    assert apo["analog_status"] == "unresolved"


def test_no_confident_oc_subgroup(apo):
    assert apo["oc_subgroup"] in (None, "unresolved")
    if apo["confident_analog"] == "ordinary_chondrite":
        assert apo["oc_subgroup"] == "unresolved"


def test_deterministic(apo):
    import c2_composition as c2
    g = pd.read_parquet(GTS)
    a = g[g["number_mp"] == 99942]
    r2 = c2.match_external_spectrum((a["wavelength_um"] * 1000).to_numpy(),
                                    a["reflectance"].to_numpy())
    assert r2["distribution"] == apo["distribution"]
    assert r2["manifold_d2"] == apo["manifold_d2"]
