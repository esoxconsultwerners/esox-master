#!/usr/bin/env python3
"""C1.1 - Composition information ceiling (synthetic RELAB, kill criterion).

Measures how separable the viable RELAB primary meteorite groups are at Gaia
16-band resolution, noise-free (information ceiling) and with realistic Gaia
noise. Produces the pairwise separability table (the real scientific output)
and the composition kill-criterion verdict. Does NOT build a matching pipeline.
"""

import itertools
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c1_common as cc

OUT = Path(__file__).resolve().parent
N_SPLITS = 5
N_NOISE = 3            # noise realizations averaged for stability
AUC_THR = 0.75
rng_master = np.random.default_rng(20260708)


def resample_extra(meta, centers):
    """Resample RELAB at extra centers (e.g. SDSS u=354, g=477), normalized by
    each spectrum's own 550 nm value - identical convention to the 16 bands."""
    SP = pd.read_parquet(cc.I / "relab_spectra.parquet")
    SP = SP[SP["relab_spectrum_id"].isin(set(meta["relab_spectrum_id"]))].copy()
    SP["wl_nm"] = SP["wavelength_um"] * 1000.0
    vals = {}
    for sid, grp in SP.groupby("relab_spectrum_id"):
        v = cc.gaussian_resample(grp["wl_nm"].to_numpy(), grp["reflectance"].to_numpy(),
                                 centers + [cc.REF_NM])
        ref = v[-1]
        vals[sid] = (v[:-1] / ref) if (np.isfinite(ref) and ref) else np.full(len(centers), np.nan)
    return np.vstack([vals[s] for s in meta["relab_spectrum_id"]])


def cv_predict_proba(clf, X, y, groups):
    gkf = GroupKFold(n_splits=N_SPLITS)
    return cross_val_predict(clf, X, y, groups=groups, cv=gkf, method="predict_proba")


def multiclass_metrics(X, y, groups, labels):
    lda = cv_predict_proba(LinearDiscriminantAnalysis(), X, y, groups)
    rf = cv_predict_proba(RandomForestClassifier(n_estimators=300, random_state=0,
                                                 n_jobs=-1), X, y, groups)
    yhat_lda = np.array(labels)[lda.argmax(1)]
    yhat_rf = np.array(labels)[rf.argmax(1)]
    yb = pd.get_dummies(pd.Categorical(y, categories=labels)).to_numpy()
    per_class_auc = {}
    for k, lab in enumerate(labels):
        try:
            per_class_auc[lab] = roc_auc_score(yb[:, k], lda[:, k])
        except ValueError:
            per_class_auc[lab] = np.nan
    return {
        "bal_acc_lda": balanced_accuracy_score(y, yhat_lda),
        "bal_acc_rf": balanced_accuracy_score(y, yhat_rf),
        "per_class_auc": per_class_auc,
        "cm": confusion_matrix(y, yhat_lda, labels=labels),
    }


def pairwise_auc(X, y, groups, labels):
    mat = pd.DataFrame(np.nan, index=labels, columns=labels)
    for a, b in itertools.combinations(labels, 2):
        m = np.isin(y, [a, b])
        Xs, ys, gs = X[m], y[m], groups[m]
        yb = (ys == b).astype(int)
        ns = min(N_SPLITS, len(np.unique(gs)))
        try:
            proba = cross_val_predict(LinearDiscriminantAnalysis(), Xs, yb, groups=gs,
                                      cv=GroupKFold(n_splits=ns), method="predict_proba")[:, 1]
            auc = roc_auc_score(yb, proba)
            auc = max(auc, 1 - auc)
        except Exception:
            auc = np.nan
        mat.loc[a, b] = mat.loc[b, a] = auc
    return mat


def avg_noisy(fn, X, y, groups, labels, band_errors):
    accum = []
    for s in range(N_NOISE):
        rng = np.random.default_rng(1000 + s)
        Xn = cc.inject_gaia_noise(X, band_errors, rng)
        accum.append(fn(Xn, y, groups, labels))
    return accum


def plot_cm(cm, labels, title, path):
    cmn = cm / cm.sum(1, keepdims=True).clip(min=1)
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(cmn, cmap="magma_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("predicted"); ax.set_ylabel("true")
    ax.set_title(title, fontsize=11)
    fig.colorbar(im, fraction=0.046, label="row-normalized (recall)")
    fig.tight_layout()
    fig.savefig(path, dpi=140, metadata={"Date": None})
    plt.close(fig)


def main():
    X, meta, labels = cc.load_relab_resampled()
    y = meta["relab_group"].to_numpy()
    groups = meta["meteorite_name"].to_numpy()
    band_errors = cc.gasp_band_errors()
    print(f"[c1.1] {len(X)} resampled spectra, {len(labels)} primary groups, "
          f"{len(np.unique(groups))} distinct meteorites")

    mc_clean = multiclass_metrics(X, y, groups, labels)
    mc_noisy = avg_noisy(lambda Xn, y, g, l: multiclass_metrics(Xn, y, g, l),
                         X, y, groups, labels, band_errors)
    bal_lda_noisy = np.mean([m["bal_acc_lda"] for m in mc_noisy])
    bal_rf_noisy = np.mean([m["bal_acc_rf"] for m in mc_noisy])
    pca_noisy = {lab: np.nanmean([m["per_class_auc"][lab] for m in mc_noisy]) for lab in labels}

    plot_cm(mc_clean["cm"], labels, "C1.1 composition confusion - noise-free (ceiling)",
            OUT / "c1_1_confusion_noisefree.png")
    plot_cm(mc_noisy[0]["cm"], labels, "C1.1 composition confusion - with Gaia noise",
            OUT / "c1_1_confusion_noisy.png")
    pd.DataFrame(mc_clean["cm"], index=labels, columns=labels).to_csv(
        OUT / "c1_1_confusion_noisefree.csv")
    pd.DataFrame(mc_noisy[0]["cm"], index=labels, columns=labels).to_csv(
        OUT / "c1_1_confusion_noisy.csv")

    pw_clean = pairwise_auc(X, y, groups, labels)
    pw_noisy_list = avg_noisy(lambda Xn, y, g, l: pairwise_auc(Xn, y, g, l),
                              X, y, groups, labels, band_errors)
    pw_noisy = sum(pw_noisy_list) / len(pw_noisy_list)
    pw_noisy.to_csv(OUT / "pairwise_auc.csv")

    sep_frac = {}
    for lab in labels:
        row = pw_noisy.loc[lab].drop(lab)
        sep_frac[lab] = float((row >= AUC_THR).mean())
    separable = {l: f for l, f in sep_frac.items() if f > 0.5}
    frac_groups_separable = len(separable) / len(labels)
    catalog_supported = frac_groups_separable >= 0.30

    def offdiag(m):
        v = m.to_numpy(); return np.nanmean(v[~np.eye(len(m), dtype=bool)])
    mean_auc_clean, mean_auc_noisy = offdiag(pw_clean), offdiag(pw_noisy)

    Xug = np.hstack([X, resample_extra(meta, [354, 477])])
    ok = np.isfinite(Xug).all(1)
    pw_ug = sum(avg_noisy(lambda Xn, y, g, l: pairwise_auc(Xn, y, g, l),
                          Xug[ok], y[ok], groups[ok], labels,
                          {b: band_errors[b] for b in cc.BANDS})) / N_NOISE
    # noise only added to the 16 Gaia bands; u/g kept (blue advantage is real signal)
    lift = offdiag(pw_ug) - offdiag(pw_noisy)

    lines = ["# C1.1 - Composition information ceiling", "",
             f"Target: {len(labels)} viable RELAB **primary** meteorite groups "
             f"({len(X)} spectra, {len(np.unique(groups))} distinct meteorites). "
             "GroupKFold by meteorite (guardrail 1); 550 nm normalization "
             "(guardrail 3); Gaussian-bandpass resampling to the 16 Gaia bands "
             "(guardrail 6, see method note in c1_common.py).", "",
             "## Overall separability", "",
             "| metric | noise-free (ceiling) | with Gaia noise |", "|---|---:|---:|",
             f"| balanced accuracy (LDA) | {mc_clean['bal_acc_lda']:.3f} | {bal_lda_noisy:.3f} |",
             f"| balanced accuracy (RandomForest) | {mc_clean['bal_acc_rf']:.3f} | {bal_rf_noisy:.3f} |",
             f"| mean pairwise ROC-AUC | {mean_auc_clean:.3f} | {mean_auc_noisy:.3f} |", "",
             "## Per-class ROC-AUC (with Gaia noise, one-vs-rest)", "",
             "| group | AUC (noisy) | pairwise-separable fraction (AUC>=0.75) |",
             "|---|---:|---:|"]
    for lab in sorted(labels, key=lambda l: -pca_noisy[l]):
        lines.append(f"| {lab} | {pca_noisy[lab]:.3f} | {sep_frac[lab]*100:.0f}% |")

    lines += ["", "## Pairwise separability (the scientific core)", "",
              f"- mean off-diagonal pairwise AUC: {mean_auc_noisy:.3f} (noisy), "
              f"{mean_auc_clean:.3f} (ceiling).",
              f"- groups separable (AUC>=0.75 vs a majority of others): "
              f"**{len(separable)}/{len(labels)}** ({frac_groups_separable*100:.0f}%).",
              f"- least separable groups (collapse at Gaia resolution): "
              f"{', '.join(sorted(sep_frac, key=sep_frac.get)[:5])}.",
              f"- most separable: {', '.join(sorted(sep_frac, key=sep_frac.get, reverse=True)[:5])}.",
              "", "Full pair matrix in pairwise_auc.csv.", "",
              "## SDSS u/g blue-extension lift", "",
              f"Adding synthetic SDSS u (354 nm) + g (477 nm) changes the mean "
              f"pairwise AUC by **{lift:+.3f}** ({int(ok.sum())} spectra cover the "
              f"blue). This quantifies the Esox-specific advantage over Gaia-only work.",
              "", "## KILL CRITERION (composition)", "",
              f"Threshold: catalog framing supported if >=30% of primary groups are "
              f"pairwise-separable (AUC>=0.75 vs most others) WITH Gaia noise. "
              f"Observed: {frac_groups_separable*100:.0f}%.", "",
              f"**VERDICT: the data support a "
              f"{'CATALOG paper (confident meteorite analogs for a well-defined subset)' if catalog_supported else 'METHODS / limits paper'}.** "
              + ("A large fraction of primary groups stay distinct at 16-band Gaia "
                 "resolution with noise; a matching pipeline is justified as a "
                 "separate later package (not built here)." if catalog_supported else
                 "Most groups collapse under Gaia noise; the honest framing is a "
                 "methods/limits paper. No matching pipeline is justified."), "",
              "No matching pipeline is built in C1 regardless of verdict (by design).", ""]
    (OUT / "c1_1_composition_ceiling.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:40]))
    print(f"\n[c1.1] catalog_supported={catalog_supported} "
          f"({frac_groups_separable*100:.0f}% separable), u/g lift {lift:+.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
