"""Guard tests for C1 degeneration analysis artifacts (no stage re-runs)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

C1 = Path(__file__).resolve().parents[1] / "analysis" / "c1_degeneration"

ARTIFACTS = [
    "c1_0_census.md", "c1_1_composition_ceiling.md", "pairwise_auc.csv",
    "c1_1_confusion_noisefree.png", "c1_1_confusion_noisy.png",
    "c1_2_weathering.md", "c1_2_auc_delta.csv",
    "c1_3_taxonomy.md", "c1_3_confusion_complex.png", "c1_3_confusion_class.png",
    "c1_4_apophis_demo.md", "C1_FINDINGS.md",
]


@pytest.mark.parametrize("name", ARTIFACTS)
def test_artifact_exists(name):
    p = C1 / name
    assert p.exists() and p.stat().st_size > 0, f"{name} missing or empty"


def test_pairwise_auc_valid():
    m = pd.read_csv(C1 / "pairwise_auc.csv", index_col=0)
    assert m.shape[0] == m.shape[1] and m.shape[0] >= 6, "pairwise matrix too small"
    v = m.to_numpy()
    off = v[~np.eye(len(v), dtype=bool)]
    off = off[np.isfinite(off)]
    assert ((off >= 0.5) & (off <= 1.0)).all(), "AUC outside [0.5, 1]"
    assert np.allclose(v, v.T, equal_nan=True), "pairwise AUC not symmetric"


def test_findings_states_verdicts():
    t = (C1 / "C1_FINDINGS.md").read_text().lower()
    assert "composition" in t and "taxonomy" in t
    assert "complex" in t and "catalog" in t and "methods" in t


def test_weathering_citation_is_placeholder():
    t = (C1 / "c1_2_weathering.md").read_text()
    assert "TODO" in t, "weathering citation must remain a TODO placeholder (no fabrication)"
