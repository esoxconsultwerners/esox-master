"""Contract tests for AP_A9 (pipeline/ingest_damit.py)."""

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
PROV = ROOT / "data" / "provenance"

MODELS = INTERIM / "damit_models.parquet"
BEST = INTERIM / "damit_best.parquet"


@pytest.fixture(scope="session")
def models():
    if not MODELS.exists():
        pytest.fail(f"{MODELS} missing - run pipeline/ingest_damit.py first")
    return pd.read_parquet(MODELS)


@pytest.fixture(scope="session")
def best():
    if not BEST.exists():
        pytest.fail(f"{BEST} missing - run pipeline/ingest_damit.py first")
    return pd.read_parquet(BEST)


def test_model_id_unique_and_referenced(models):
    assert models["damit_model_id"].is_unique, "damit_model_id not unique"
    ref = models["damit_ref"].astype("string").str.strip()
    missing = models[ref.isna() | (ref == "")]
    assert missing.empty, f"{len(missing)} model rows without damit_ref"


def test_counts(models):
    assert len(models) > 4000, f"only {len(models)} models"
    n_ast = models.assign(
        _k=models["number_mp"].astype("string").fillna(models["designation"])
    )["_k"].nunique()
    assert n_ast > 3000, f"only {n_ast} distinct asteroids"


def test_value_ranges(models):
    lam = models["damit_lambda_deg"].dropna()
    beta = models["damit_beta_deg"].dropna()
    per = models["damit_period_h"].dropna()
    assert ((lam >= 0) & (lam < 360)).all(), (
        f"lambda outside [0, 360): {lam[(lam < 0) | (lam >= 360)].tolist()[:5]}"
    )
    assert ((beta >= -90) & (beta <= 90)).all(), (
        f"beta outside [-90, 90]: {beta[(beta < -90) | (beta > 90)].tolist()[:5]}"
    )
    assert ((per > 0.5) & (per < 5000)).all(), (
        f"period outside (0.5, 5000): {per[(per <= 0.5) | (per >= 5000)].tolist()[:5]}"
    )


def test_lcdb_cross_check(best):
    lcdb = pd.read_parquet(INTERIM / "lcdb.parquet")[["number_mp", "lcdb_period_h"]].dropna()
    lcdb = lcdb.drop_duplicates("number_mp")
    xc = best[["number_mp", "damit_period_h"]].dropna().merge(lcdb, on="number_mp")
    assert len(xc) > 100, f"only {len(xc)} asteroids shared with LCDB"
    rel = ((xc["damit_period_h"] - xc["lcdb_period_h"]).abs() / xc["lcdb_period_h"])
    med = float(rel.median())
    assert med < 0.005, (
        f"median relative period difference {med:.5f} >= 0.5% - unit/parsing error?"
    )


def test_best_one_row_per_asteroid(best):
    key = best["number_mp"].astype("string").fillna(best["designation"])
    assert key.is_unique, "damit_best has more than one row for some asteroid"
    assert len(best) > 3000, f"damit_best only {len(best)} rows"


def test_provenance_complete():
    path = PROV / "damit.json"
    assert path.exists(), f"{path} missing"
    doc = json.loads(path.read_text())
    for k in ("source", "url", "retrieval_date", "row_counts", "citation",
              "license_terms", "method_derivation"):
        assert doc.get(k) not in (None, "", [], {}), f"damit.json missing {k}"
    assert "CC-BY" in doc["license_terms"], "license not CC-BY"
