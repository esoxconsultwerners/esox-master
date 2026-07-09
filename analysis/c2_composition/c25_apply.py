"""C2.5 apply - honest composition layer (Werner decision: complex-level +
honest-unresolved core).

Regenerates the Gaia-core analog_* columns with the count-unbiased per-group
matcher (c25_matcher) plus an absolute-fit MANIFOLD gate. Per the C2.5 finding
(diag/C25_MATCHER_FINDING.md) the Gaia core sits ~8x beyond the RELAB manifold,
so confident Gaia-only meteorite-group analogs are not honestly supportable:
the core is reported overwhelmingly analog_status='unresolved', the soft
analog_distribution is kept as an INDICATOR ONLY, and confident 'ok' is emitted
only where an object genuinely falls inside the RELAB manifold. The primary
composition product is the complex-level taxon_esox (C4). Append-only: only the
analog_* / properties columns are rewritten; all other bytes preserved.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "analysis" / "c1_degeneration"))
sys.path.insert(0, str(HERE))
import c1_common as cc
import c2_composition as c2
import c25_matcher as m5

FINAL = cc.ROOT / "data" / "final"
PROV = cc.ROOT / "data" / "provenance"
CORE = FINAL / "esox_master_core.parquet"
CONF_THR = 0.55
MANIFOLD_Q = 0.95


def build():
    _, meta, _ = cc.load_relab_resampled()
    V16 = c2.relab_features(meta, cc.BANDS, cc.REF_NM)
    ok = np.isfinite(V16).all(1)
    V16, meta = V16[ok], meta[ok].reset_index(drop=True)
    y = meta["relab_group"].map(c2.merge_label).to_numpy()
    groups = meta["meteorite_name"].to_numpy()
    F = m5.features_from_16(V16)
    mm = m5.GroupGaussianMatcher().fit(F, y)
    ba, ece, T, conf, corr = mm.fit_calibration(F, y, groups)
    cidx = {c: i for i, c in enumerate(mm.classes_)}
    d2_self = np.array([(F[i] - mm.mu_[cidx[y[i]]]) @ mm.prec_ @ (F[i] - mm.mu_[cidx[y[i]]])
                        for i in range(len(y))])
    q95 = float(np.quantile(d2_self, MANIFOLD_Q))
    return mm, ba, ece, T, q95


def main():
    mm, ba, ece, T, q95 = build()
    classes = mm.classes_
    sep, deg, sep_merged = c2.separable_sets()
    deg_merged = {c2.merge_label(x) for x in deg}

    core = pd.read_parquet(CORE)
    before = core["analog_status"].value_counts().to_dict() if "analog_status" in core else {}
    before_cm = int(((core.get("analog_group_top") == "CM") &
                     (core.get("analog_status") == "ok")).sum()) if "analog_status" in core else 0

    g = pd.read_parquet(cc.GASP_CAT)
    BAND_COLS = [f"refl_{b}" for b in cc.BANDS]
    X = g.set_index("number_mp").loc[core["number_mp"], BAND_COLS].to_numpy()
    F = m5.features_from_16(X)
    raw = mm.predict_proba(F)
    cx = core["taxon_esox"].astype(object).to_numpy()

    tops, softs, confs, dists, ocsub, status, d2s, evid = [], [], [], [], [], [], [], []
    for i in range(len(core)):
        post, _ = c2.apply_prior(raw[i], classes, cx[i])
        ti = int(post.argmax())
        soft = classes[ti]
        d2 = float((F[i] - mm.mu_[ti]) @ mm.prec_ @ (F[i] - mm.mu_[ti]))
        conf = float(mm.calibrate([float(post.max())])[0])
        within = d2 <= q95
        if soft in deg_merged:
            st, top = "degenerate", "degenerate"
            ev = "degenerate_sink"
        elif within and conf >= CONF_THR:
            st, top = "ok", soft
            ev = "gaia_within_manifold"
        else:
            st, top = "unresolved", "unresolved"
            ev = "offmanifold" if not within else "low_confidence"
        tops.append(top); softs.append(soft); confs.append(round(conf, 6))
        dists.append(c2.distribution_json(post, classes, sep_merged))
        ocsub.append("unresolved" if top == "ordinary_chondrite" else None)
        status.append(st); d2s.append(round(d2, 3)); evid.append(ev)

    core["analog_group_top"] = tops
    core["analog_top_soft"] = softs
    core["analog_group_conf"] = confs
    core["analog_distribution"] = dists
    core["oc_subgroup"] = ocsub
    core["analog_status"] = status
    core["analog_manifold_d2"] = d2s
    core["analog_evidence"] = evid

    c2.properties_bridge(core, sep_merged, deg_merged)
    core.to_parquet(CORE, index=False)

    after = pd.Series(status).value_counts().to_dict()
    after_cm = int(((np.array(tops) == "CM")).sum())
    n = len(core)
    n_ok = after.get("ok", 0)

    write_provenance(classes, sep_merged, deg_merged, ba, ece, T, q95, n_ok, n)
    write_rerun(before, before_cm, after, after_cm, n, ba, ece, T, q95, status, tops, cx)

    print("=== C2.5 apply complete ===")
    print(f"matcher: shared-cov Gaussian, edge-trimmed 13-band, bal_acc {ba:.3f}, ECE {ece:.3f}, T {T:.1f}")
    print(f"manifold q95 d2 = {q95:.1f}")
    print(f"BEFORE (RF C2): {before}  CM-confident {before_cm}")
    print(f"AFTER (honest): {after}  CM-top {after_cm}")
    print(f"confident 'ok' = {n_ok} ({100*n_ok/n:.2f}% of core); unresolved = {after.get('unresolved',0)} ({100*after.get('unresolved',0)/n:.1f}%)")
    return 0


def write_provenance(classes, sep_merged, deg_merged, ba, ece, T, q95, n_ok, n):
    (PROV / "analog_composition.json").write_text(json.dumps({
        "columns": ["analog_group_top", "analog_top_soft", "analog_group_conf",
                    "analog_distribution", "oc_subgroup", "analog_status",
                    "analog_manifold_d2", "analog_evidence", "weathering_strength",
                    "weathering_flag", "inferred_grain_density_range",
                    "inferred_porosity_note", "properties_source"],
        "kind": "MODEL-DERIVED from RELAB meteorite templates (resampled to the 16 "
                "Gaia bands under the documented 44 nm Gaussian bandpass "
                "approximation; edge bands 374/1034 nm dropped for matching). "
                "Distinct from measured quantities.",
        "matcher": "C2.5 shared-covariance Gaussian (per-group centroid + one pooled "
                   "Ledoit-Wolf covariance -> Mahalanobis; uniform class priors so "
                   "template count cannot win; temperature-scaled, isotonic-calibrated). "
                   "Replaces the C2 RandomForest, which was a template-density catch-all "
                   f"(diag/C25_MATCHER_FINDING.md). self-CV bal_acc {round(ba,3)}, ECE {round(ece,3)}, T {round(T,1)}.",
        "honesty_policy": "The Gaia DR3 asteroid-reflectance core sits ~8x beyond the "
                          "RELAB manifold (median Mahalanobis d2 ~189 vs RELAB q95 "
                          f"{round(q95,1)}); a rigid de-bias does not close it. Confident "
                          "Gaia-only meteorite-group analogs are therefore NOT honestly "
                          "supportable. analog_status is overwhelmingly 'unresolved'; the "
                          "soft analog_distribution is an INDICATOR ONLY; 'ok' is emitted "
                          "only inside the RELAB manifold (d2<=q95) and above conf threshold. "
                          "Primary composition product is complex-level taxon_esox (C4). "
                          f"confident 'ok' fraction {round(100*n_ok/n,2)}%.",
        "manifold_gate_q95_d2": round(q95, 3),
        "confidence_threshold": CONF_THR,
        "separable_subset": sorted(sep_merged),
        "degenerate_groups": sorted(deg_merged),
        "oc_policy": "H/L/LL merged into ordinary_chondrite; oc_subgroup always 'unresolved'",
        "external_nir_path": "confident analogs for objects with ground-truth NIR spectra "
                             "(Apophis etc.) come via c2_composition.match_external_spectrum; "
                             "migrating that path to the C2.5 matcher is a documented follow-up.",
    }, indent=2))


def write_rerun(before, before_cm, after, after_cm, n, ba, ece, T, q95, status, tops, cx):
    tops = np.array(tops, dtype=object); status = np.array(status)
    okm = status == "ok"
    okvc = pd.Series(tops[okm]).value_counts().to_dict()
    S = cx == "S"
    s_cm = int(((tops == "CM") & okm & S).sum())
    def tbl(d):
        return "".join(f"| {k} | {v:,} |\n" for k, v in sorted(d.items(), key=lambda x: -x[1]))
    (HERE / "c25_3_rerun.md").write_text(
        "# C2.5 Step 3 - honest re-run (complex-level + unresolved core)\n\n"
        f"New matcher: shared-covariance Gaussian, edge-trimmed 13-band, uniform priors, "
        f"manifold-gated. self-CV bal_acc {ba:.3f}, ECE {ece:.3f}, T {T:.1f}, manifold q95 d2 {q95:.1f}.\n\n"
        f"## Before (RF C2 density catch-all)\n\n| status | n |\n|---|---:|\n" + tbl(before) +
        f"\nConfident CM: **{before_cm:,}** (97% of confident).\n\n"
        f"## After (honest)\n\n| status | n |\n|---|---:|\n" + tbl(after) +
        f"\nConfident 'ok': **{after.get('ok',0):,}** ({100*after.get('ok',0)/n:.2f}% of core). "
        f"Unresolved: **{after.get('unresolved',0):,}** ({100*after.get('unresolved',0)/n:.1f}%). "
        f"CM as a confident top: **{int((tops[okm]=='CM').sum())}** (was {before_cm:,}).\n\n"
        f"## Confident 'ok' top-analog distribution\n\n| analog | n |\n|---|---:|\n" + tbl(okvc) +
        f"\nS-complex getting confident CM: **{s_cm}** (target 0). The core is now "
        f"overwhelmingly unresolved - the honest state given the Gaia<->RELAB domain "
        f"divergence. The composition product users should rely on is the complex-level "
        f"taxon_esox (C4, 57% coverage, 81% held-out agreement); the meteorite-analog "
        f"distribution is a soft indicator, and confident group-level analogs are reserved "
        f"for in-manifold objects and the external NIR ground-truth path (Apophis).\n")


if __name__ == "__main__":
    sys.exit(main())
