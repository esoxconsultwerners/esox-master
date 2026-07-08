"""Contract tests for AP_A7 (pipeline/ingest_relab.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ROOT / "data" / "interim" / "relab_samples.parquet"
SPECTRA = ROOT / "data" / "interim" / "relab_spectra.parquet"
GROUP_MAP = ROOT / "data" / "interim" / "relab_group_map.csv"


@pytest.fixture(scope="session")
def samples():
    if not SAMPLES.exists():
        pytest.fail(f"{SAMPLES} missing - run pipeline/ingest_relab.py first")
    return pd.read_parquet(SAMPLES)


@pytest.fixture(scope="session")
def spectra():
    if not SPECTRA.exists():
        pytest.fail(f"{SPECTRA} missing - run pipeline/ingest_relab.py first")
    return pd.read_parquet(SPECTRA)


def test_spectrum_id_unique(samples):
    assert samples["relab_spectrum_id"].is_unique


def test_samples_and_spectra_consistent(samples, spectra):
    sample_ids = set(samples["relab_spectrum_id"])
    spectra_ids = set(spectra["relab_spectrum_id"].unique())
    assert sample_ids == spectra_ids, (
        f"{len(sample_ids - spectra_ids)} samples without points, "
        f"{len(spectra_ids - sample_ids)} orphan spectra"
    )


def test_kept_counts(samples):
    assert len(samples) > 800, f"only {len(samples)} kept spectra"
    n_met = samples["meteorite_name"].str.lower().nunique()
    assert n_met > 300, f"only {n_met} distinct meteorites"


def test_viable_groups(samples):
    counts = samples[samples["relab_group"] != "unmapped"]["relab_group"].value_counts()
    viable = counts[counts >= 5]
    assert len(viable) >= 12, f"only {len(viable)} groups with >= 5 spectra"


@pytest.mark.parametrize("group", ["LL", "CM", "CV", "H", "L", "eucrite"])
def test_key_groups_covered(samples, group):
    # After the intermediate-class soft-label split, "LL-like coverage" means
    # the union {LL, L-LL} (the L/LL intermediates moved to the L-LL soft
    # group). The pure-LL count dropped from 91 to 89 but stays well above 5.
    groups = {"LL": ["LL", "L-LL"]}.get(group, [group])
    n = int(samples["relab_group"].isin(groups).sum())
    assert n >= 5, f"group(s) {groups} have only {n} spectra"


def test_value_ranges(spectra):
    r = spectra["reflectance"]
    wl = spectra["wavelength_um"]
    assert ((r > 0) & (r < 1.5)).all(), (
        f"reflectance outside (0, 1.5): {r[(r <= 0) | (r >= 1.5)].head().tolist()}"
    )
    # Kept spectra retain their full archived range (no truncation):
    # BDR data starts at 0.275 um and 10 FTIR-combined products extend
    # to ~100 um. Sanity bounds are therefore wide; the VIS/NIR window
    # (0.3-3.0 um) must still dominate, and the Gaia span 0.40-1.00 um
    # is asserted per spectrum in test_gaia_range_span.
    assert ((wl > 0.2) & (wl < 120)).all(), (
        f"wavelengths outside (0.2, 120) um: "
        f"{wl[(wl <= 0.2) | (wl >= 120)].head().tolist()}"
    )
    frac_visnir = ((wl > 0.3) & (wl < 3.0)).mean()
    assert frac_visnir > 0.95, f"only {frac_visnir:.1%} of points in 0.3-3.0 um"


def test_gaia_range_span(samples):
    bad = samples[(samples["wl_min_um"] > 0.40) | (samples["wl_max_um"] < 1.00)]
    assert bad.empty, (
        f"{len(bad)} kept spectra do not span 0.40-1.00 um: "
        f"{bad['relab_spectrum_id'].head().tolist()}"
    )


def test_grain_sizes(samples):
    both = samples.dropna(subset=["grain_size_min_um", "grain_size_max_um"])
    assert (both["grain_size_min_um"] <= both["grain_size_max_um"]).all(), \
        "grain_size_min > grain_size_max"
    for col in ("grain_size_min_um", "grain_size_max_um"):
        v = samples[col].dropna()
        assert ((v >= 0) & (v < 10_000)).all(), f"{col} outside [0, 10000)"


def test_unmapped_fraction(samples):
    frac = (samples["relab_group"] == "unmapped").mean()
    assert frac < 0.20, f"unmapped fraction {frac:.1%} >= 20%"


def test_group_map_exists():
    assert GROUP_MAP.exists()
    gm = pd.read_csv(GROUP_MAP)
    assert {"raw", "canonical", "group_kind"} <= set(gm.columns)
    assert len(gm) > 10
    assert set(gm["group_kind"]) <= {"primary", "soft", "coarse"}, \
        f"unexpected group_kind values: {set(gm['group_kind'])}"
    # the intermediate soft-label groups must be present and flagged soft
    soft = gm[gm["group_kind"] == "soft"]["canonical"].unique()
    assert "L-LL" in soft and "H-L" in soft, f"soft groups missing: {soft}"


def test_provenance():
    path = ROOT / "data" / "provenance" / "relab.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for k in ("url", "urn", "version", "date", "row_count", "funnel",
              "citation", "license_terms"):
        assert doc.get(k), f"provenance missing {k}"
    assert doc["row_count"] > 800
