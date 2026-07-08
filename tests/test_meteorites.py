"""Contract tests for AP_A8 (pipeline/ingest_meteorites.py).

MetBull source actually used: NASA Open Data "Meteorite Landings" (the frozen
~2013 fallback), because the live Meteoritical Bulletin search is behind a
Cloudflare JS challenge. Entry-count bound below therefore uses the > 40,000
fallback threshold, not the > 60,000 live threshold.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
PROV = ROOT / "data" / "provenance"

METBULL = INTERIM / "metbull.parquet"
ORBITS = INTERIM / "meteorite_orbits.parquet"
CNEOS = INTERIM / "cneos_fireballs.parquet"

CANONICAL_GROUPS = {
    "H", "L", "LL", "OC-ung", "CI", "CM", "CO", "CV", "CK", "CR", "CH", "CB",
    "C-ungrouped", "E", "EH", "EL", "R", "K", "eucrite", "howardite",
    "diogenite", "aubrite", "ureilite", "angrite", "brachinite", "winonaite",
    "lodranite-acapulcoite", "iron", "mesosiderite", "pallasite",
    "lunar-meteorite", "martian",
    # soft-label groups for genuine ordinary-chondrite intermediates
    # (e.g. the L/LL5 falls Golden and Dingle Dell); see group-map fix.
    "L-LL", "H-L", "H-LL",
}


@pytest.fixture(scope="session")
def metbull():
    if not METBULL.exists():
        pytest.fail(f"{METBULL} missing - run pipeline/ingest_meteorites.py first")
    return pd.read_parquet(METBULL)


@pytest.fixture(scope="session")
def orbits():
    if not ORBITS.exists():
        pytest.fail(f"{ORBITS} missing - run pipeline/ingest_meteorites.py first")
    return pd.read_parquet(ORBITS)


@pytest.fixture(scope="session")
def cneos():
    if not CNEOS.exists():
        pytest.fail(f"{CNEOS} missing - run pipeline/ingest_meteorites.py first")
    return pd.read_parquet(CNEOS)


def test_metbull_name_unique(metbull):
    assert metbull["metbull_name"].is_unique, "metbull_name not unique after dedup"


def test_metbull_entry_count(metbull):
    # NASA Open Data fallback (frozen ~2013); live MetBull would be > 60,000.
    assert len(metbull) > 40_000, (
        f"only {len(metbull)} MetBull entries - expected > 40,000 with the "
        f"NASA fallback"
    )


def test_metbull_falls_count(metbull):
    falls = int((metbull["metbull_fall_find"] == "fall").sum())
    assert 1000 <= falls <= 2000, f"falls count {falls} outside [1000, 2000]"


def test_orbits_min_events(orbits):
    assert len(orbits) >= 40, f"only {len(orbits)} orbit events"


def test_orbits_every_row_has_ref(orbits):
    missing = orbits[orbits["morb_orbit_ref"].isna()
                     | (orbits["morb_orbit_ref"].astype(str).str.strip() == "")]
    assert missing.empty, (
        f"{len(missing)} orbit rows without morb_orbit_ref: "
        f"{missing['morb_meteorite_name'].tolist()}"
    )


def test_orbits_element_ranges(orbits):
    a = orbits["morb_a_au"].dropna()
    e = orbits["morb_e"].dropna()
    i = orbits["morb_i_deg"].dropna()
    assert ((a > 0.5) & (a < 6)).all(), f"a outside (0.5, 6): {a[(a <= 0.5) | (a >= 6)].tolist()}"
    assert ((e >= 0) & (e < 1)).all(), f"e outside [0, 1): {e[(e < 0) | (e >= 1)].tolist()}"
    assert ((i >= 0) & (i < 90)).all(), f"i outside [0, 90): {i[(i < 0) | (i >= 90)].tolist()}"


def test_orbits_class_maps_to_group(orbits):
    bad = orbits[~orbits["morb_group"].isin(CANONICAL_GROUPS)]
    assert bad.empty, (
        f"{len(bad)} orbit events do not map to a canonical group: "
        f"{list(zip(bad['morb_meteorite_name'], bad['morb_class'], bad['morb_group']))}"
    )


def test_cneos_min_rows(cneos):
    assert len(cneos) > 500, f"only {len(cneos)} CNEOS fireballs"


def test_provenance_complete():
    for name in ("metbull", "meteorite_orbits", "cneos_fireballs"):
        path = PROV / f"{name}.json"
        assert path.exists(), f"{path} missing"
        doc = json.loads(path.read_text())
        for k in ("source", "url", "retrieval_date", "row_count", "citation",
                  "license_terms"):
            assert doc.get(k) not in (None, "", []), f"{name}.json missing {k}"


def test_metbull_source_documented():
    doc = json.loads((PROV / "metbull.json").read_text())
    assert "staleness" in doc and doc["staleness"], "metbull staleness not documented"
    assert "NASA Open Data" in doc["source"], "metbull source not the NASA fallback"
