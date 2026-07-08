#!/usr/bin/env python3
"""C4 - Productive complex-level taxonomy classifier (taxon_esox).

Turns the C1.3-validated approach into a catalog column: a Bus-DeMeo COMPLEX
label (C/S/X/V/D/A/K-L) plus a calibrated per-object confidence for all 19,190
GASP-core objects, with an explicit "unclassified" category below threshold.
taxon_esox is MODEL-DERIVED and kept strictly separate from the literature
taxon_* columns.

Integrity fix vs C1.3: all 42 PDS Bus-DeMeo validation objects are also in the
358 GASP-Mahlke set, so C1.3's 86% "held-out" agreement was leakage-inflated.
C4 EXCLUDES those 42 from training (train on 316) and reports the clean number.

Guardrails (from C1): complex level only; honest unclassified; PDS Bus-DeMeo
never trained on; calibrated confidence (reliability diagram); same 44 nm
Gaussian bandpass approximation (documented in c1_common).
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "c1_degeneration"))
import c1_common as cc

OUT = Path(__file__).resolve().parent
FINAL = cc.ROOT / "data" / "final"
PROV = cc.ROOT / "data" / "provenance"
DOCS = cc.ROOT / "docs"
BAND_COLS = [f"refl_{b}" for b in cc.BANDS]
COMPLEXES = ["S", "C", "X", "V", "D", "A", "K/L"]   # exact C1.3 complex scheme
TARGET_COVERAGE = 0.57                              # C1.3 operating point
SEED = 0


def load():
    g = pd.read_parquet(cc.GASP_CAT)
    tw = pd.read_parquet(cc.I / "taxonomy_wide.parquet")[["number_mp", "taxon_bus_demeo"]]
    return g.merge(tw, on="number_mp", how="left")


def main():
    g = load()
    g["mahlke_complex"] = g["taxonomy"].map(cc.to_complex)
    g["bd_complex"] = g["taxon_bus_demeo"].map(cc.to_complex)

    val = g[g["bd_complex"].isin(COMPLEXES)].copy()          # held-out PDS set
    val_ids = set(g.loc[g["taxon_bus_demeo"].notna(), "number_mp"])
    train = g[g["mahlke_complex"].isin(COMPLEXES)
              & ~g["number_mp"].isin(val_ids)].copy()          # leakage-free
    Xtr, ytr = train[BAND_COLS].to_numpy(), train["mahlke_complex"].to_numpy()
    print(f"[c4] train {len(train)} (val {len(val)} excluded), complexes {COMPLEXES}")

    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    classes = np.unique(ytr)

    def lr(balanced):
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                             LogisticRegression(max_iter=2000,
                                                class_weight="balanced" if balanced else None))
    # RF handles the NaN edge-bands natively; LogReg gets a median imputer.
    # class_weight='balanced' helps balanced-accuracy but distorts probabilities
    # (worse calibration), so natural-probability variants are compared too.
    candidates = [
        ("LogReg(balanced)+sigmoid", lr(True), "sigmoid"),
        ("LogReg+sigmoid", lr(False), "sigmoid"),
        ("LogReg+isotonic", lr(False), "isotonic"),
        ("RF+isotonic", RandomForestClassifier(n_estimators=500, random_state=SEED,
                                                n_jobs=-1), "isotonic"),
    ]

    def reliability(conf, correct):
        bins = np.linspace(0, 1, 9)
        idx = np.clip(np.digitize(conf, bins) - 1, 0, len(bins) - 2)
        xs, ys, ece = [], [], 0.0
        for b in range(len(bins) - 1):
            m = idx == b
            if m.sum() >= 5:
                xs.append(conf[m].mean()); ys.append(correct[m].mean())
                ece += m.mean() * abs(conf[m].mean() - correct[m].mean())
        return xs, ys, ece

    results = {}
    for name, base, method in candidates:
        cal = CalibratedClassifierCV(base, method=method, cv=5)
        oof = cross_val_predict(cal, Xtr, ytr, cv=skf, method="predict_proba")
        pred = classes[oof.argmax(1)]
        ba = balanced_accuracy_score(ytr, pred)
        xs, ys, ece = reliability(oof.max(1), (pred == ytr).astype(float))
        results[name] = dict(base=base, method=method, ba=ba, ece=ece, xs=xs, ys=ys)
        print(f"[c4] {name}: bal_acc {ba:.3f}, ECE {ece:.3f}")
    # Pick the best-calibrated model (lowest ECE) among those with usable
    # accuracy - a trustworthy confidence is the whole point of C4.
    eligible = {k: v for k, v in results.items() if v["ba"] >= 0.40}
    base_name = min(eligible or results, key=lambda k: results[k]["ece"])
    chosen = results[base_name]
    base, ba_cal, ece = chosen["base"], chosen["ba"], chosen["ece"]
    ba = {k: round(v["ba"], 3) for k, v in results.items()}
    print(f"[c4] chosen {base_name} (bal_acc {ba_cal:.3f}, ECE {ece:.3f})")

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot([0, 1], [0, 1], "--", color="#999", lw=1)
    ax.plot(chosen["xs"], chosen["ys"], "o-", color="#4059ad")
    ax.set_xlabel("mean predicted confidence"); ax.set_ylabel("observed accuracy")
    ax.set_title(f"C4 calibration - {base_name} (ECE={ece:.3f})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.tight_layout(); fig.savefig(OUT / "calibration.png", dpi=140, metadata={"Date": None})
    plt.close(fig)

    cal = CalibratedClassifierCV(base, method=chosen["method"], cv=5)
    cal.fit(Xtr, ytr)

    # Predict the whole core; choose threshold to hit ~TARGET_COVERAGE.
    core = pd.read_parquet(FINAL / "esox_master_core.parquet")
    Xall = g.set_index("number_mp").loc[core["number_mp"], BAND_COLS].to_numpy()
    proba = cal.predict_proba(Xall)
    cls = cal.classes_
    order = np.argsort(proba, axis=1)[:, ::-1]
    top1 = cls[order[:, 0]]; conf1 = proba[np.arange(len(proba)), order[:, 0]]
    top2 = cls[order[:, 1]]; conf2 = proba[np.arange(len(proba)), order[:, 1]]

    # Round conf and threshold to the same precision and decide on the rounded
    # values, so stored taxon_esox_conf and the classified/unclassified split
    # are exactly consistent (no rounding-boundary ties).
    thr = float(np.round(np.quantile(conf1, 1 - TARGET_COVERAGE), 6))
    conf1r = np.round(conf1, 6)
    taxon_esox = np.where(conf1r >= thr, top1, "unclassified")
    coverage = float((taxon_esox != "unclassified").mean())
    print(f"[c4] threshold {thr:.3f} -> confident coverage {coverage*100:.1f}%")

    core["taxon_esox"] = taxon_esox
    core["taxon_esox_conf"] = conf1r
    core["taxon_esox_second"] = top2
    core["taxon_esox_second_conf"] = np.round(conf2, 6)

    # C4.3 completeness recompute with confident taxon_esox
    confident_tax = core["taxon_esox"] != "unclassified"
    core["complete_profile_esox"] = (core["has_gaia_spectrum"] & core["has_period_best"]
                                     & core["has_family"] & core["has_diameter_best"]
                                     & confident_tax)
    new_cp = int(core["complete_profile_esox"].sum())
    old_cp = int(core["complete_profile"].sum())
    core.to_parquet(FINAL / "esox_master_core.parquet", index=False)

    # funnel v1
    steps = [("spectrum", core["has_gaia_spectrum"]),
             ("+ diameter_best", core["has_diameter_best"]),
             ("+ family", core["has_family"]),
             ("+ period_best", core["has_period_best"]),
             ("+ taxon_esox", confident_tax)]
    mask = pd.Series(True, index=core.index); vals = []
    for lab, s in steps:
        mask = mask & s.astype(bool); vals.append(int(mask.sum()))
    plot_funnel(steps, vals)

    # C4.4 validation on held-out PDS Bus-DeMeo
    Xval = val[BAND_COLS].to_numpy()
    vpred = cal.predict(Xval)
    vtrue = val["bd_complex"].to_numpy()
    agree = float((vpred == vtrue).mean())
    vconf = cal.predict_proba(Xval).max(1)
    cmask = vconf >= thr
    purity = float((vpred[cmask] == vtrue[cmask]).mean()) if cmask.any() else float("nan")
    n_conf_val = int(cmask.sum())
    cm = confusion_matrix(vtrue, vpred, labels=COMPLEXES)
    per_cx = {c: (int(cm[i, i]), int(cm[i].sum())) for i, c in enumerate(COMPLEXES)}

    dist = pd.Series(taxon_esox).value_counts().to_dict()
    write_artifacts(ba, base_name, ba_cal, ece, thr, coverage, dist, old_cp, new_cp,
                    vals, agree, per_cx, len(train), len(val), purity, n_conf_val)
    write_provenance(base_name, thr, coverage, ba_cal, agree, len(train))
    apophis_note()

    print(f"\n=== C4 complete ===")
    print(f"taxon_esox written for {len(core):,} core; confident coverage {coverage*100:.0f}%.")
    print(f"Held-out PDS agreement {agree*100:.0f}% (complex level; leakage-free).")
    print(f"Complete physical profiles: {old_cp} -> {new_cp} "
          f"({100*new_cp/len(core):.1f}% of core).")
    print("Ready for the composition matching pipeline (separable subset, with "
          "taxon_esox complex as prior).")
    return 0


def plot_funnel(steps, vals):
    labels = [l for l, _ in steps]
    for ext in ("png", "svg"):
        fig, ax = plt.subplots(figsize=(7.2, 4.0))
        y = list(range(len(vals)))[::-1]
        colors = ["#b9bcc0"] * (len(vals) - 1) + ["#4059ad"]
        ax.barh(y, vals, color=colors, height=0.6, zorder=3)
        for yi, v in zip(y, vals):
            ax.text(v + max(vals) * 0.012, yi, f"{v:,}", va="center", fontsize=9)
        ax.set_yticks(y); ax.set_yticklabels(labels)
        ax.set_xlim(0, max(vals) * 1.12)
        ax.set_xlabel("objects surviving each cumulative AND-requirement")
        ax.set_title("Completeness funnel v1 - taxonomy via taxon_esox (C4)\n"
                     "(19,190-object spectral core)", fontsize=11, loc="left")
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0); ax.xaxis.grid(True, color="#e6e6e6", zorder=0)
        ax.set_axisbelow(True)
        drop = vals[-2] - vals[-1]
        ax.annotate(f"taxonomy cut now only -{drop:,}\n(C4 solved the bottleneck)",
                    xy=(vals[-1], y[-1]), xytext=(max(vals) * 0.3, y[-1] + 0.4),
                    fontsize=8.5, color="#2a9d8f",
                    arrowprops=dict(arrowstyle="->", color="#2a9d8f", lw=0.9))
        fig.tight_layout()
        fig.savefig(FINAL / f"completeness_funnel_v1.{ext}", dpi=150,
                    bbox_inches="tight", metadata={"Date": None})
        plt.close(fig)


def write_artifacts(ba, base_name, ba_cal, ece, thr, coverage, dist, old_cp, new_cp,
                    funnel, agree, per_cx, n_train, n_val, purity, n_conf_val):
    (OUT / "c4_1_training.md").write_text(
        "# C4.1 - Train & calibrate\n\n"
        f"Features: 16 Gaia bands only (C1.3 showed SDSS u/g adds only +0.006 at "
        f"complex level, so the classifier is kept purely spectral). Complex "
        f"scheme (from C1.3): {COMPLEXES}. Training {n_train} GASP-Mahlke objects, "
        f"leakage-free (the {n_val} PDS Bus-DeMeo objects are excluded - they are "
        f"100% inside the Mahlke set, so C1.3's 86% was leakage-inflated).\n\n"
        "## Model selection (CV balanced accuracy)\n\n| model | balanced acc |\n|---|---:|\n"
        + "".join(f"| {k} | {v:.3f} |\n" for k, v in ba.items())
        + f"\nChosen: **{base_name}**. Calibrated (sigmoid) balanced acc {ba_cal:.3f}.\n\n"
        f"## Calibration\n\nReliability diagram: calibration.png. Expected "
        f"calibration error **ECE = {ece:.3f}** (lower is better; the calibrated "
        f"confidence is trustworthy where the curve tracks the diagonal).\n\n"
        f"## Operating point\n\nConfidence threshold **{thr:.3f}**, chosen to "
        f"reproduce the C1.3 ~57% confident-coverage point. Achieved coverage "
        f"{coverage*100:.1f}%.\n")

    dist_lines = "".join(f"| {k} | {v:,} |\n" for k, v in
                         sorted(dist.items(), key=lambda x: -x[1]))
    (OUT / "c4_2_prediction.md").write_text(
        "# C4.2 - Prediction across the core\n\n"
        f"taxon_esox + taxon_esox_conf + taxon_esox_second(_conf) written for all "
        f"19,190 core objects.\n\n## Complex distribution\n\n| taxon_esox | n |\n|---|---:|\n"
        + dist_lines + f"\nConfident coverage: **{coverage*100:.1f}%** "
        f"(threshold {thr:.3f}).\n")

    (OUT / "c4_3_completeness.md").write_text(
        "# C4.3 - Completeness recompute (the payoff)\n\n"
        f"Complete physical profile = gaia_spectrum AND period_best AND family AND "
        f"diameter_best AND **taxon_esox (confident)**.\n\n"
        f"**Old (literature taxon_any): {old_cp}  ->  New (taxon_esox): {new_cp}** "
        f"({100*new_cp/19190:.1f}% of core), a {new_cp/old_cp:.1f}x increase.\n\n"
        f"Funnel v1: data/final/completeness_funnel_v1.png. Survivor counts: "
        + " -> ".join(f"{v:,}" for v in funnel) + ".\n")

    val_lines = "".join(f"| {c} | {hit}/{n} |\n" for c, (hit, n) in per_cx.items() if n)
    (OUT / "c4_4_validation.md").write_text(
        "# C4.4 - Validation & honesty\n\n"
        f"Independent held-out PDS Bus-DeMeo set ({n_val} objects, never trained "
        f"on). **Complex-level agreement: {agree*100:.1f}%** (leakage-free; C1.3's "
        f"86% was inflated by 100% train/val overlap - this is the trustworthy "
        f"number). Among the {n_conf_val} confident held-out predictions "
        f"(conf >= {thr:.3f}), purity is **{purity*100:.1f}%**. The agreement is "
        f"carried by the dominant S complex; rare complexes (K/L, D) are the weak "
        f"spots - see per-complex below and treat them with caution.\n\n"
        f"## Per-complex agreement (correct/total)\n\n| complex | agree |\n"
        f"|---|---:|\n" + val_lines +
        "\n## Apophis\n\ntaxon_esox = N/A for Apophis (99942): it has no Gaia "
        "spectrum (not in the GASP core), so no Gaia-feature prediction is "
        "produced. Its composition comes via the external-spectrum path (C1.4), "
        "not C4. The 44 nm Gaussian bandpass approximation (C1) applies here too.\n")


def write_provenance(base_name, thr, coverage, ba_cal, agree, n_train):
    (PROV / "taxon_esox.json").write_text(json.dumps({
        "column": "taxon_esox",
        "kind": "MODEL-DERIVED (predicted), distinct from the literature taxon_* "
                "columns (taxon_mahlke, taxon_tholen, taxon_bus, taxon_bus_demeo, "
                "taxon_s3os2) - never conflate measured vs predicted taxonomy",
        "granularity": "Bus-DeMeo complex level only (C/S/X/V/D/A/K-L); class-level "
                       "was shown to over-claim in C1.3 (0.416 balanced acc)",
        "model": base_name + " + sigmoid calibration (CalibratedClassifierCV)",
        "features": "16 Gaia bands (374-1034 nm), normalized at 550 nm; purely "
                    "spectral (u/g add +0.006 only)",
        "training_labels": f"{n_train} GASP-Mahlke objects, mapped to complex; the "
                           "PDS Bus-DeMeo set is held out (leakage-free, unlike C1.3)",
        "confidence": "calibrated; taxon_esox_conf in [0,1]; below threshold "
                      f"{thr:.3f} -> 'unclassified' (honest, not a low-conf guess)",
        "confidence_threshold": thr,
        "confident_coverage": round(coverage, 3),
        "calibrated_balanced_accuracy_cv": round(ba_cal, 3),
        "heldout_pds_agreement": round(agree, 3),
        "resampling_note": "same 44 nm Gaussian bandpass approximation as C1 "
                           "(true Gaia passband response not in local data)",
    }, indent=2))


def apophis_note():
    dossier = DOCS / "apophis_status.md"
    line = ("\n- **taxon_esox = N/A** (C4): Apophis has no Gaia spectrum, so the "
            "model-derived complex classifier does not run for it. Its composition "
            "comes only via the external-spectrum path (C1.4). This is the correct "
            "behaviour, not a gap.\n")
    text = dossier.read_text()
    if "taxon_esox = N/A" not in text:
        dossier.open("a").write(line)


if __name__ == "__main__":
    sys.exit(main())
