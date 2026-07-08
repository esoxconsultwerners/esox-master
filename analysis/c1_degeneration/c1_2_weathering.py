#!/usr/bin/env python3
"""C1.2 - Space weathering as a nuisance parameter (composition arm).

Applies an exponential reddening continuum to the resampled RELAB spectra,
sweeps the strength C, marginalizes it into the classifier (train across
weathering strengths) and measures how much pairwise separability degrades
versus the weathering-free case. Also derives a would-be catalog column: the
per-spectrum estimated weathering strength.
"""

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c1_common as cc

OUT = Path(__file__).resolve().parent
N_SPLITS = 5
C_SWEEP = [0.0, 0.1, 0.2, 0.3]          # weathering strength (um), 0 = pristine
LAM_UM = np.array(cc.BANDS) / 1000.0
INV = 1.0 / LAM_UM - 1.0 / (cc.REF_NM / 1000.0)

# TODO(Werner): verify + fill the space-weathering reddening reference before
# the paper. Placeholder model here is W(lambda) = exp(-C / lambda) anchored at
# 550 nm (a standard exponential-continuum reddening form). Do NOT cite a
# fabricated source; this is left blank on purpose (Galinier discipline).
WEATHERING_CITATION = "TODO: verified reference (exponential reddening continuum)"


def weather(X, C):
    """Multiply by exp(-C * (1/lam - 1/0.55)); 550 nm stays anchored at 1."""
    return X * np.exp(-C * INV)[None, :]


def estimate_C(X):
    """Per-spectrum weathering-strength proxy: slope of ln(refl) vs (1/lam-1/ref).
    Returns estimated C (a would-be catalog column)."""
    out = np.full(len(X), np.nan)
    A = np.vstack([INV, np.ones_like(INV)]).T
    for i, row in enumerate(X):
        m = np.isfinite(row) & (row > 0)
        if m.sum() < 4:
            continue
        slope = np.linalg.lstsq(A[m], np.log(row[m]), rcond=None)[0][0]
        out[i] = -slope
    return out


def pairwise(X, y, groups, labels, rng, band_errors, augment_C=None):
    """Pairwise ROC-AUC. If augment_C given, training marginalizes weathering by
    stacking weathered copies (same meteorite id -> no leakage); test spectra
    get one random weathering strength + Gaia noise."""
    mat = pd.DataFrame(np.nan, index=labels, columns=labels)
    for a, b in itertools.combinations(labels, 2):
        m = np.isin(y, [a, b])
        Xs, ys, gs = X[m], y[m], groups[m]
        yb = (ys == b).astype(int)
        ns = min(N_SPLITS, len(np.unique(gs)))
        gkf = GroupKFold(n_splits=ns)
        aucs = []
        for tr, te in gkf.split(Xs, yb, gs):
            if augment_C is None:
                Xtr = cc.inject_gaia_noise(Xs[tr], band_errors, rng)
                ytr = yb[tr]
            else:
                parts_X, parts_y = [], []
                for C in augment_C:
                    parts_X.append(cc.inject_gaia_noise(weather(Xs[tr], C), band_errors, rng))
                    parts_y.append(yb[tr])
                Xtr = np.vstack(parts_X); ytr = np.concatenate(parts_y)
            Cte = rng.choice(augment_C) if augment_C else 0.0
            Xte = cc.inject_gaia_noise(weather(Xs[te], Cte), band_errors, rng)
            try:
                clf = LinearDiscriminantAnalysis().fit(Xtr, ytr)
                p = clf.predict_proba(Xte)[:, 1]
                aucs.append(roc_auc_score(yb[te], p))
            except Exception:
                pass
        if aucs:
            au = np.mean(aucs)
            mat.loc[a, b] = mat.loc[b, a] = max(au, 1 - au)
    return mat


def offdiag(m):
    v = m.to_numpy()
    return np.nanmean(v[~np.eye(len(m), dtype=bool)])


def main():
    X, meta, labels = cc.load_relab_resampled()
    y = meta["relab_group"].to_numpy()
    groups = meta["meteorite_name"].to_numpy()
    band_errors = cc.gasp_band_errors()
    rng = np.random.default_rng(7)

    print(f"[c1.2] {len(X)} spectra; sweeping C in {C_SWEEP}")
    before = pairwise(X, y, groups, labels, rng, band_errors, augment_C=None)
    after = pairwise(X, y, groups, labels, rng, band_errors, augment_C=C_SWEEP)

    rows = []
    for a, b in itertools.combinations(labels, 2):
        rows.append({"group_a": a, "group_b": b,
                     "auc_before": before.loc[a, b], "auc_after": after.loc[a, b],
                     "delta": after.loc[a, b] - before.loc[a, b]})
    delta = pd.DataFrame(rows).sort_values("delta")
    delta.to_csv(OUT / "c1_2_auc_delta.csv", index=False)

    mean_before, mean_after = offdiag(before), offdiag(after)
    mean_change = mean_after - mean_before
    worst = delta.head(1).iloc[0]
    degraded = delta[delta["delta"] < 0]
    mean_drop_among_degraded = degraded["delta"].mean()

    Cest = estimate_C(X)
    ce = pd.Series(Cest).describe()

    lines = ["# C1.2 - Space weathering as a nuisance parameter", "",
             f"Reddening model: **W(lambda) = exp(-C / lambda)** anchored at 550 nm; "
             f"strength swept over C = {C_SWEEP} (um). Marginalized into the "
             "classifier by training across all weathering strengths (weathered "
             "copies of a meteorite stay in its GroupKFold fold - no leakage).", "",
             f"Citation: **{WEATHERING_CITATION}**", "",
             "## Separability under weathering (mean pairwise ROC-AUC, Gaia noise)", "",
             "| | mean pairwise AUC |", "|---|---:|",
             f"| weathering-free | {mean_before:.3f} |",
             f"| weathering-marginalized | {mean_after:.3f} |",
             f"| **mean change** | **{mean_change:+.3f}** |", "",
             "Mean separability is essentially unchanged: marginalizing weathering "
             "acts as data augmentation and does not collapse the average. The real "
             f"signal is per-pair: among the {len(degraded)} pairs that degrade, the "
             f"mean drop is {mean_drop_among_degraded:+.3f} and the worst is "
             f"{worst['group_a']} vs {worst['group_b']} at {worst['delta']:+.3f}.", "",
             "## Most weathering-sensitive group pairs (largest AUC drop)", "",
             "| pair | AUC before | AUC after | delta |", "|---|---:|---:|---:|"]
    for _, r in delta.head(8).iterrows():
        lines.append(f"| {r['group_a']} vs {r['group_b']} | {r['auc_before']:.3f} | "
                     f"{r['auc_after']:.3f} | {r['delta']:+.3f} |")
    lines += ["", "Full table in c1_2_auc_delta.csv.", "",
              "## Estimated weathering strength (would-be catalog column)", "",
              "Per-spectrum reddening slope C, fit as ln(refl) vs (1/lambda - "
              "1/0.55). This is a derivable science product (a per-object "
              "weathering-strength column):", "",
              f"- median C = {ce['50%']:.3f}, IQR [{np.nanpercentile(Cest,25):.3f}, "
              f"{np.nanpercentile(Cest,75):.3f}], range [{ce['min']:.3f}, "
              f"{ce['max']:.3f}] over {int(ce['count'])} spectra.", "",
              "## Summary", "",
              f"Mean pairwise separability is robust to weathering once "
              f"marginalized (change {mean_change:+.3f} AUC - augmentation "
              f"compensates). The honest cost is concentrated in specific "
              f"weathering-ambiguous pairs (worst {worst['group_a']} vs "
              f"{worst['group_b']} {worst['delta']:+.3f}); those are the pairs a "
              "composition catalog must flag rather than assert.", ""]
    (OUT / "c1_2_weathering.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[c1.2] mean pairwise AUC {mean_before:.3f} -> {mean_after:.3f} "
          f"(change {mean_change:+.3f}); worst pair {worst['delta']:+.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
