"""Contract tests for CM v0 (analysis/coverage_matrix.py)."""

from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "data" / "final"
CSV = FINAL / "coverage_matrix_v0.csv"
CORE = FINAL / "esox_master_core.parquet"

EXPECTED_PROPS = {"orbit", "h_best", "phase_curve", "gaia_spectrum", "diameter_best",
                  "albedo_best", "family", "period_best", "taxon_any",
                  "groundtruth_spectrum", "complete_physical_profile"}
FIGURES = ["coverage_matrix_v0.png", "coverage_matrix_v0.svg",
           "completeness_funnel_v0.png", "completeness_funnel_v0.svg"]


@pytest.fixture(scope="session")
def table():
    if not CSV.exists():
        pytest.fail(f"{CSV} missing - run analysis/coverage_matrix.py first")
    return pd.read_csv(CSV)


def test_all_properties_present(table):
    assert set(table["property"]) == EXPECTED_PROPS, (
        f"missing/extra: {EXPECTED_PROPS ^ set(table['property'])}"
    )


def test_pct_in_range(table):
    assert ((table["pct_core"] >= 0) & (table["pct_core"] <= 100)).all()


def test_key_counts(table):
    g = table.set_index("property")["n_core"]
    assert g["gaia_spectrum"] == 19190, f"gaia_spectrum {g['gaia_spectrum']} != 19190"
    assert g["complete_physical_profile"] == 204, (
        f"complete profile {g['complete_physical_profile']} != 204"
    )


def test_funnel_monotonic_ends_at_204():
    core = pd.read_parquet(CORE)
    steps = ["has_gaia_spectrum", "has_diameter_best", "has_family",
             "has_period_best", "has_taxon_any"]
    mask = pd.Series(True, index=core.index)
    survivors = []
    for c in steps:
        mask = mask & core[c].astype(bool)
        survivors.append(int(mask.sum()))
    assert all(a >= b for a, b in zip(survivors, survivors[1:])), (
        f"funnel not monotonic non-increasing: {survivors}"
    )
    assert survivors[0] == 19190 and survivors[-1] == 204, survivors


def test_figures_exist_nonempty():
    for f in FIGURES:
        p = FINAL / f
        assert p.exists(), f"{p} missing"
        assert p.stat().st_size > 0, f"{p} empty"


def test_deterministic_table():
    # The table is a pure function of the core parquet: recomputing must match.
    core = pd.read_parquet(CORE)
    n = len(core)
    tbl = pd.read_csv(CSV).set_index("property")
    checks = {"gaia_spectrum": "has_gaia_spectrum", "diameter_best": "has_diameter_best",
              "family": "has_family", "period_best": "has_period_best",
              "taxon_any": "has_taxon_any"}
    for prop, col in checks.items():
        expect = int(core[col].astype(bool).sum())
        assert tbl.loc[prop, "n_core"] == expect, f"{prop}: {tbl.loc[prop,'n_core']} != {expect}"
        assert abs(tbl.loc[prop, "pct_core"] - round(100 * expect / n, 1)) < 0.05
