#!/usr/bin/env python3
"""C1.4 - External-spectrum demonstration (the Apophis path).

Apophis has no Gaia spectrum (not in the GASP core), only MITHNEOS NIR spectra
(770-2485 nm). It can therefore be served only via the external-spectrum path.
This stage proves the composition machinery works on external ground truth:
restrict to the 7 Gaia bands >= 770 nm that Apophis covers (normalized at
770 nm, since 550 nm is unavailable), train the RELAB matcher on those bands,
and classify Apophis. Both a validation (literature says Sq -> LL/ordinary
chondrite) and the showcase mechanism.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c1_common as cc

OUT = Path(__file__).resolve().parent
DOCS = cc.ROOT / "docs"
NIR_BANDS = [770, 814, 858, 902, 946, 990, 1034]
NIR_REF = 770
OC = {"H", "L", "LL"}


def nir_resample(wl_nm, refl):
    v = cc.gaussian_resample(wl_nm, refl, NIR_BANDS)
    ri = NIR_BANDS.index(NIR_REF)
    return v / v[ri] if (np.isfinite(v[ri]) and v[ri]) else np.full(len(NIR_BANDS), np.nan)


def relab_nir(meta):
    SP = pd.read_parquet(cc.I / "relab_spectra.parquet")
    SP = SP[SP["relab_spectrum_id"].isin(set(meta["relab_spectrum_id"]))].copy()
    SP["wl_nm"] = SP["wavelength_um"] * 1000.0
    feats = {}
    for sid, grp in SP.groupby("relab_spectrum_id"):
        feats[sid] = nir_resample(grp["wl_nm"].to_numpy(), grp["reflectance"].to_numpy())
    X = np.vstack([feats[s] for s in meta["relab_spectrum_id"]])
    ok = np.isfinite(X).all(1)
    return X[ok], meta[ok].reset_index(drop=True)


def main():
    _, meta, labels = cc.load_relab_resampled()
    X, meta = relab_nir(meta)
    y = meta["relab_group"].to_numpy()
    print(f"[c1.4] RELAB NIR training: {len(X)} spectra x {len(NIR_BANDS)} bands")

    gts = pd.read_parquet(cc.I / "groundtruth_spectra.parquet")
    apo = gts[gts["number_mp"] == 99942]
    a_vec = nir_resample((apo["wavelength_um"] * 1000).to_numpy(),
                         apo["reflectance"].to_numpy())
    if not np.isfinite(a_vec).all():
        print("[c1.4] Apophis NIR resample incomplete"); return 2

    clf = RandomForestClassifier(n_estimators=600, random_state=0, n_jobs=-1,
                                 class_weight="balanced").fit(X, y)
    proba = clf.predict_proba(a_vec.reshape(1, -1))[0]
    order = np.argsort(proba)[::-1]
    classes = clf.classes_
    top = [(classes[i], proba[i]) for i in order[:6]]
    oc_mass = float(sum(proba[i] for i, c in enumerate(classes) if c in OC))
    ll_p = float(proba[list(classes).index("LL")]) if "LL" in classes else 0.0
    top_analog = top[0][0]

    lines = ["# C1.4 - Apophis external-spectrum composition demo", "",
             f"Apophis (99942) has no Gaia spectrum; served via the external path "
             f"from {len(apo)} MITHNEOS NIR points ("
             f"{apo['wavelength_um'].min()*1000:.0f}-{apo['wavelength_um'].max()*1000:.0f} nm). "
             f"Matcher restricted to the {len(NIR_BANDS)} Gaia bands >= 770 nm "
             "(normalized at 770 nm), trained on the RELAB primary groups.", "",
             "## Composition analog distribution (top 6)", "",
             "| rank | analog group | probability |", "|---|---|---:|"]
    for k, (g, p) in enumerate(top, 1):
        lines.append(f"| {k} | {g} | {p:.3f} |")
    lines += ["", f"- Ordinary-chondrite (H+L+LL) probability mass: **{oc_mass:.3f}**.",
              f"- LL probability: {ll_p:.3f}. Top analog: **{top_analog}**.", "",
              "## Interpretation", "",
              f"The literature classifies Apophis as **Sq**, the spectral bridge to "
              f"**LL ordinary chondrites**. The external-spectrum matcher's top "
              f"analog is **{top_analog}** with ordinary-chondrite mass {oc_mass:.2f}. "
              + ("This is consistent with the Sq->LL/ordinary-chondrite expectation "
                 "and validates the external-spectrum mechanism."
                 if oc_mass >= 0.4 or top_analog in OC else
                 "This is only partly consistent - the NIR-only, 7-band, no-2um-band "
                 "restriction limits ordinary-chondrite discrimination.") +
              " Caveat: only the 770-1034 nm Gaia window is used (no visible slope, "
              "no 2 um pyroxene band), so this is a lower bound on what a full "
              "external spectrum would resolve.", ""]
    (OUT / "c1_4_apophis_demo.md").write_text("\n".join(lines) + "\n")

    block = ["\n## Composition analog (external-spectrum path)", "",
             f"Via C1.4: Apophis has no Gaia spectrum, so its composition analog is "
             f"recovered from MITHNEOS NIR spectra resampled to the 7 Gaia bands "
             f">= 770 nm. Top analog **{top_analog}** (p={top[0][1]:.2f}); "
             f"ordinary-chondrite mass {oc_mass:.2f}, LL p={ll_p:.2f}. "
             + ("Consistent with the literature Sq -> LL story; the external-spectrum "
                "showcase mechanism works." if oc_mass >= 0.4 or top_analog in OC else
                "Partly consistent given the NIR-only 7-band restriction.")
             + " (Lower bound: no visible slope, no 2 um band at Gaia resolution.)", ""]
    dossier = DOCS / "apophis_status.md"
    heading = "## Composition analog (external-spectrum path)"
    text = dossier.read_text()
    if heading in text:
        text = text[:text.index(heading)].rstrip() + "\n"
        dossier.write_text(text)
    dossier.open("a").write("\n".join(block) + "\n")

    print("\n".join(lines))
    print(f"\n[c1.4] top analog {top_analog}, OC mass {oc_mass:.3f}, LL p {ll_p:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
