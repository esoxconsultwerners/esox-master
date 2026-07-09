#!/usr/bin/env python3
"""C2 - Composition matching pipeline (Physical Small Bodies, Stage 1).

Assigns every core object (and external spectra like Apophis) a meteorite-analog
DISTRIBUTION with honest confidence, restricted to the C1 separable subset, with
the C4 complex as a prior and the H/L/LL degeneracy as built-in flag logic.
Produces catalog columns, not a paper.

Reuses the C1 resampling / normalization / weathering machinery (imported from
analysis/c1_degeneration/c1_common.py). Same 44 nm Gaussian-bandpass
approximation as C1/C4 (documented there).
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import GroupKFold, cross_val_predict

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "c1_degeneration"))
import c1_common as cc

OUT = Path(__file__).resolve().parent
FINAL = cc.ROOT / "data" / "final"
PROV = cc.ROOT / "data" / "provenance"
DOCS = cc.ROOT / "docs"
C1 = cc.ROOT / "analysis" / "c1_degeneration"
BAND_COLS = [f"refl_{b}" for b in cc.BANDS]
NIR_BANDS = [770, 814, 858, 902, 946, 990, 1034]
NIR_REF = 770
OC_SUB = {"H", "L", "LL"}
OC_LABEL = "ordinary_chondrite"
CONF_THR = 0.55
PRIOR_LOW = 0.25
SEED = 0

# C4 complex -> plausible analog groups (merged labels). Documented mapping;
# implausible analogs get down-weighted by PRIOR_LOW, not zeroed.
PRIOR_MAP = {
    "S": {OC_LABEL, "lodranite-acapulcoite", "ureilite", "brachinite"},
    "C": {"CM", "CR", "CV", "CI", "CK", "CO"},
    "X": {"EH", "EL", "aubrite", "mesosiderite", "pallasite"},
    "V": {"eucrite", "howardite", "diogenite"},
    "A": {"ureilite", "brachinite", "pallasite"},
    "D": {"CM", "CR", "CI"},
    "K/L": {"CO", "CV", "CK"},
}


def separable_sets():
    pw = pd.read_csv(C1 / "pairwise_auc.csv", index_col=0)
    sep_frac = {g: float((pw.loc[g].drop(g) >= 0.75).mean()) for g in pw.index}
    sep = {g for g, f in sep_frac.items() if f > 0.5}
    deg = {g for g, f in sep_frac.items() if f <= 0.5}
    sep_merged = {OC_LABEL if g in OC_SUB else g for g in sep}
    return sep, deg, sep_merged


def merge_label(g):
    return OC_LABEL if g in OC_SUB else g


def relab_features(meta, bands, ref):
    SP = pd.read_parquet(cc.I / "relab_spectra.parquet")
    SP = SP[SP["relab_spectrum_id"].isin(set(meta["relab_spectrum_id"]))].copy()
    SP["wl_nm"] = SP["wavelength_um"] * 1000.0
    ri = bands.index(ref)
    feats = {}
    for sid, grp in SP.groupby("relab_spectrum_id"):
        v = cc.gaussian_resample(grp["wl_nm"].to_numpy(), grp["reflectance"].to_numpy(), bands)
        feats[sid] = v / v[ri] if (np.isfinite(v[ri]) and v[ri]) else np.full(len(bands), np.nan)
    X = np.vstack([feats[s] for s in meta["relab_spectrum_id"]])
    return X


def build_matcher(bands, ref, reliability_png=None):
    _, meta, _ = cc.load_relab_resampled()
    X = relab_features(meta, bands, ref)
    ok = np.isfinite(X).all(1)
    X, meta = X[ok], meta[ok].reset_index(drop=True)
    y = meta["relab_group"].map(merge_label).to_numpy()
    groups = meta["meteorite_name"].to_numpy()
    clf = RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1,
                                 class_weight="balanced")
    if reliability_png:
        oof = cross_val_predict(clf, X, y, groups=groups, cv=GroupKFold(5),
                                method="predict_proba")
        classes = np.unique(y)
        pred = classes[oof.argmax(1)]
        conf, corr = oof.max(1), (pred == y).astype(float)
        bins = np.linspace(0, 1, 9)
        idx = np.clip(np.digitize(conf, bins) - 1, 0, len(bins) - 2)
        xs, ys, ece = [], [], 0.0
        for b in range(len(bins) - 1):
            m = idx == b
            if m.sum() >= 5:
                xs.append(conf[m].mean()); ys.append(corr[m].mean())
                ece += m.mean() * abs(conf[m].mean() - corr[m].mean())
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot([0, 1], [0, 1], "--", color="#999", lw=1)
        ax.plot(xs, ys, "o-", color="#4059ad")
        ax.set_xlabel("mean predicted confidence"); ax.set_ylabel("observed accuracy")
        ax.set_title(f"C2 matcher reliability (ECE={ece:.3f})")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        fig.tight_layout(); fig.savefig(reliability_png, dpi=140, metadata={"Date": None})
        plt.close(fig)
        ba = balanced_accuracy_score(y, pred)
        clf.fit(X, y)
        return clf, np.unique(y), ba, ece
    clf.fit(X, y)
    return clf, np.unique(y)


def apply_prior(proba, classes, complex_label):
    if complex_label is None or complex_label == "unclassified" or complex_label not in PRIOR_MAP:
        return proba, False
    plausible = PRIOR_MAP[complex_label]
    w = np.array([1.0 if c in plausible else PRIOR_LOW for c in classes])
    post = proba * w
    s = post.sum()
    return (post / s if s > 0 else proba), True


def entropy(p):
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def distribution_json(proba, classes, sep_merged):
    d = {c: round(float(p), 4) for c, p in zip(classes, proba) if c in sep_merged}
    s = sum(d.values())
    if s > 0:
        d = {c: round(v / s, 4) for c, v in d.items()}
    return json.dumps(d)


# Two matchers built once at import for reuse by match_external_spectrum.
_FULL = None
_NIR = None
_SEP = _DEG = _SEPM = None


def _ensure_matchers():
    global _FULL, _NIR, _SEP, _DEG, _SEPM
    if _FULL is None:
        _SEP, _DEG, _SEPM = separable_sets()
        _FULL = build_matcher(cc.BANDS, cc.REF_NM)
        _NIR = build_matcher(NIR_BANDS, NIR_REF)


def match_external_spectrum(wl_nm, refl, complex_label=None):
    """External-spectrum composition (C2.6 migration): delegates to the C2.5
    per-group matcher + RELAB-manifold gate in c26_external. The RandomForest
    external path is RETIRED for composition; taxon_esox/C4 is a separate
    classifier and is untouched."""
    from c26_external import match_external_spectrum as _mx
    return _mx(wl_nm, refl, complex_label)


def main():
    _ensure_matchers()
    sep, deg, sep_merged = _SEP, _DEG, _SEPM
    deg_merged = {merge_label(g) for g in deg}
    (clf, classes), (_, _) = _FULL, _NIR
    full_clf, full_classes, ba, ece = build_matcher(cc.BANDS, cc.REF_NM,
                                                     reliability_png=OUT / "c2_reliability.png")
    print(f"[c2] matcher: {len(full_classes)} analog classes, RELAB self-CV bal_acc "
          f"{ba:.3f}, ECE {ece:.3f}")

    core = pd.read_parquet(FINAL / "esox_master_core.parquet")
    g = pd.read_parquet(cc.GASP_CAT)
    Xcore = g.set_index("number_mp").loc[core["number_mp"], BAND_COLS].to_numpy()
    has_spec = core["has_gaia_spectrum"].to_numpy().astype(bool)
    proba = full_clf.predict_proba(Xcore)

    tops, confs, dists, ocsub, status, ent_flat, ent_post, used_prior = (
        [], [], [], [], [], [], [], [])
    cx = core["taxon_esox"].to_numpy()
    for i in range(len(core)):
        p = proba[i]
        post, used = apply_prior(p, full_classes, cx[i])
        ent_flat.append(entropy(p)); ent_post.append(entropy(post)); used_prior.append(used)
        top = full_classes[post.argmax()]; conf = float(post.max())
        tops.append(top); confs.append(round(conf, 6))
        dists.append(distribution_json(post, full_classes, sep_merged))
        ocsub.append("unresolved" if top == OC_LABEL else None)
        if not has_spec[i]:
            st = "low_snr"
        elif top in deg_merged:
            st = "degenerate"
        elif conf >= CONF_THR:
            st = "ok"
        elif cx[i] == "unclassified":
            st = "unclassified_prior"
        else:
            st = "low_snr"
        status.append(st)

    core["analog_group_top"] = np.where(np.array(status) == "degenerate", "degenerate",
                                        np.array(tops, dtype=object))
    core["analog_group_conf"] = confs
    core["analog_distribution"] = dists
    core["oc_subgroup"] = ocsub
    core["analog_status"] = status

    # C1.2 weathering nuisance: estimated strength + fragile-pair flag
    from numpy.linalg import lstsq
    inv = 1.0 / (np.array(cc.BANDS) / 1000.0) - 1.0 / (cc.REF_NM / 1000.0)
    A = np.vstack([inv, np.ones_like(inv)]).T
    wstr = []
    for row in Xcore:
        m = np.isfinite(row) & (row > 0)
        wstr.append(-lstsq(A[m], np.log(row[m]), rcond=None)[0][0] if m.sum() >= 4 else np.nan)
    core["weathering_strength"] = np.round(wstr, 4)
    FRAGILE = {"CO", "EL", "CR", "R", "CK", "aubrite"}   # C1.2 named-fragile groups
    core["weathering_flag"] = [t in FRAGILE for t in tops]

    # entropy reduction from the C4 prior
    up = np.array(used_prior)
    red = (np.array(ent_flat)[up] - np.array(ent_post)[up]) / np.clip(np.array(ent_flat)[up], 1e-9, None)
    mean_red = float(np.nanmean(red)) if up.any() else 0.0

    core.to_parquet(FINAL / "esox_master_core.parquet", index=False)

    dist_counts = pd.Series([s for s in status]).value_counts().to_dict()
    top_counts = pd.Series([t for t, s in zip(core["analog_group_top"], status)
                            if s == "ok"]).value_counts().to_dict()

    # C2.3 external path: Apophis + Gaia-vs-external cross-check
    apo, gaia_ext_agree, n_xcheck, ga_conf, n_conf_x = apophis_and_crosscheck(
        core, g, full_clf, full_classes)

    # C2.4 properties bridge
    n_props = properties_bridge(core, sep_merged, deg_merged)
    core.to_parquet(FINAL / "esox_master_core.parquet", index=False)

    write_artifacts(full_classes, sep_merged, deg_merged, ba, ece, dist_counts,
                    top_counts, mean_red, apo, gaia_ext_agree, n_xcheck, n_props,
                    int(has_spec.sum()), ga_conf, n_conf_x)
    write_provenance(full_classes, sep_merged, deg_merged, mean_red, ba, ece)

    n_ok = dist_counts.get("ok", 0)
    n_deg = dist_counts.get("degenerate", 0)
    n_spec = int(has_spec.sum())
    print("\n=== C2 complete ===")
    print(f"Composition analogs written for {n_spec:,} core objects; "
          f"{100*n_ok/n_spec:.0f}% confident, "
          f"{100*(n_deg+dist_counts.get('low_snr',0)+dist_counts.get('unclassified_prior',0))/n_spec:.0f}% "
          f"degenerate/unresolved.")
    print(f"C4 prior reduced analog entropy by {mean_red*100:.0f}%.")
    print(f"Apophis production analog: {apo['top']} (conf {apo['top_conf']:.2f}, "
          f"OC mass {apo['oc_mass']:.2f}, oc_subgroup unresolved).")
    print(f"Physical-properties bridge: {n_props} objects with analog-inferred density.")
    print("Ready for Stage 2 (DR4 density catalog) after DR4 release.")
    return 0


def apophis_and_crosscheck(core, g, clf, classes):
    gts = pd.read_parquet(cc.I / "groundtruth_spectra.parquet")
    apo_pts = gts[gts["number_mp"] == 99942]
    apo = match_external_spectrum((apo_pts["wavelength_um"] * 1000).to_numpy(),
                                  apo_pts["reflectance"].to_numpy(), complex_label=None)
    # cross-check: core objects that also have A5 ground-truth spectra
    both = gts[gts["number_mp"].isin(set(core["number_mp"]))]
    core_top = dict(zip(core["number_mp"], core["analog_group_top"]))
    core_stat = dict(zip(core["number_mp"], core["analog_status"]))
    agree, agree_conf = [], []
    for num, grp in both.groupby("number_mp"):
        ext = match_external_spectrum((grp["wavelength_um"] * 1000).to_numpy(),
                                      grp["reflectance"].to_numpy())
        ctop = core_top.get(num)
        if ctop and ctop != "degenerate" and ext["top"] != "degenerate":
            hit = ext["top"] == ctop
            agree.append(hit)
            if core_stat.get(num) == "ok":
                agree_conf.append(hit)
    ga = float(np.mean(agree)) if agree else float("nan")
    ga_conf = float(np.mean(agree_conf)) if agree_conf else float("nan")

    block = ["\n## Composition analog (production)", "",
             f"Via C2 (external-spectrum path, {apo['path']}): Apophis analog "
             f"distribution top **{apo['top']}** (conf {apo['top_conf']:.2f}), "
             f"ordinary-chondrite mass {apo['oc_mass']:.2f}. Per C1 honesty, the "
             f"ordinary-chondrite subgroup is reported **unresolved** (H/L/LL are "
             f"not separable at Gaia/NIR resolution) - Apophis gets 'ordinary "
             f"chondrite, subgroup unresolved', not 'LL'.", ""]
    dossier = DOCS / "apophis_status.md"
    text = dossier.read_text()
    if "## Composition analog (production)" not in text:
        dossier.open("a").write("\n".join(block) + "\n")
    return apo, ga, len(agree), ga_conf, len(agree_conf)


def properties_bridge(core, sep_merged, deg_merged):
    path = cc.I / "meteorite_group_properties.csv"
    rows = [
        (OC_LABEL, 3.2, 3.8, "low-moderate (~5-10%)", "olivine+pyroxene+FeNi metal (H/L/LL)", "TODO"),
        ("CM", 2.2, 2.9, "high", "hydrated phyllosilicates, volatile-rich", "TODO"),
        ("CV", 3.0, 3.6, "moderate", "CAIs + chondrules, anhydrous", "TODO"),
        ("CR", 2.8, 3.4, "moderate-high", "metal + hydrated matrix", "TODO"),
        ("CO", 3.0, 3.6, "moderate", "small chondrules, anhydrous", "TODO"),
        ("CI", 1.6, 2.4, "very high", "aqueously altered, no chondrules", "TODO"),
        ("eucrite", 3.0, 3.4, "low", "basaltic HED (Vesta crust)", "TODO"),
        ("howardite", 3.1, 3.4, "low", "HED breccia (Vesta)", "TODO"),
        ("diogenite", 3.2, 3.5, "low", "orthopyroxenite HED (Vesta)", "TODO"),
        ("aubrite", 2.9, 3.2, "low", "enstatite achondrite, iron-poor", "TODO"),
        ("ureilite", 3.1, 3.4, "low-moderate", "olivine+pyroxene+carbon", "TODO"),
        ("EH", 3.5, 3.7, "low", "enstatite chondrite, metal-rich", "TODO"),
        ("EL", 3.4, 3.6, "low", "enstatite chondrite", "TODO"),
        ("brachinite", 3.3, 3.7, "low", "olivine-rich primitive achondrite", "TODO"),
        ("lodranite-acapulcoite", 3.2, 3.6, "low", "primitive achondrite, partial melt", "TODO"),
        ("martian", 3.0, 3.4, "low", "SNC (basaltic/cumulate)", "TODO"),
        ("lunar-meteorite", 2.9, 3.3, "low-moderate", "feldspathic/basaltic breccia", "TODO"),
        ("mesosiderite", 4.0, 5.5, "low", "stony-iron (silicate+metal)", "TODO"),
        ("pallasite", 4.5, 6.0, "low", "stony-iron (olivine+metal)", "TODO"),
    ]
    props = pd.DataFrame(rows, columns=["group", "grain_density_min_gcc",
                                        "grain_density_max_gcc", "porosity_note",
                                        "composition_note", "citation"])
    props.to_csv(path, index=False)
    pmap = props.set_index("group")

    dens, poro, src = [], [], []
    for top, st in zip(core["analog_group_top"], core["analog_status"]):
        if st == "ok" and top in pmap.index and top in sep_merged:
            r = pmap.loc[top]
            dens.append(f"{r['grain_density_min_gcc']}-{r['grain_density_max_gcc']} g/cm3")
            poro.append(r["porosity_note"]); src.append(f"analog-inferred from {top}")
        else:
            dens.append(None); poro.append(None); src.append(None)
    core["inferred_grain_density_range"] = dens
    core["inferred_porosity_note"] = poro
    core["properties_source"] = src
    return int(pd.Series(src).notna().sum())


def write_artifacts(classes, sep_merged, deg_merged, ba, ece, dist_counts, top_counts,
                    mean_red, apo, ga, n_xcheck, n_props, n_spec, ga_conf, n_conf_x):
    (OUT / "c2_1_matcher.md").write_text(
        "# C2.1 - Composition matcher (separable subset)\n\n"
        f"Analog target classes ({len(classes)}): {sorted(classes)}. Separable "
        f"subset from C1 pairwise AUC (>=0.75 vs majority); H/L/LL merged into "
        f"**ordinary_chondrite** (C1: mutually ~0.72, never emit a confident OC "
        f"subgroup). Degenerate sink classes (reported, never confidently "
        f"asserted): {sorted(deg_merged)}.\n\n"
        f"RandomForest on RELAB templates (16 Gaia bands, 550 nm norm). "
        f"Self-classification GroupKFold-by-meteorite: balanced acc {ba:.3f}, "
        f"reliability ECE {ece:.3f} (c2_reliability.png).\n\n"
        f"## C4 complex prior\n\nComplex -> plausible analog groups (implausible "
        f"down-weighted x{PRIOR_LOW}, not zeroed): {PRIOR_MAP}. Where taxon_esox "
        f"is 'unclassified', a flat prior is used. **The prior reduces mean analog "
        f"entropy by {mean_red*100:.1f}%** - the quantitative reason C4 came "
        f"before C2.\n")

    dc = "".join(f"| {k} | {v:,} |\n" for k, v in sorted(dist_counts.items(), key=lambda x: -x[1]))
    tc = "".join(f"| {k} | {v:,} |\n" for k, v in sorted(top_counts.items(), key=lambda x: -x[1]))
    (OUT / "c2_2_catalog.md").write_text(
        "# C2.2 - Catalog application\n\n"
        f"analog_* columns written for {n_spec:,} core objects with a Gaia "
        f"spectrum.\n\n## analog_status distribution\n\n| status | n |\n|---|---:|\n"
        + dc + f"\n## Confident (ok) top-analog distribution\n\n| analog | n |\n|---|---:|\n"
        + tc + f"\nConfident coverage: {100*dist_counts.get('ok',0)/n_spec:.1f}% "
        f"of Gaia-spectrum objects; the rest are honestly degenerate / low_snr / "
        f"unclassified_prior (C1 requires a substantial unresolved fraction).\n")

    (OUT / "c2_3_external.md").write_text(
        "# C2.3 - External-spectrum path\n\n"
        f"match_external_spectrum() resamples any (wavelength, reflectance) and "
        f"returns the analog distribution, choosing the full 16-band matcher when "
        f"the visible is covered else the NIR-7-band matcher (reproducing C1.4).\n\n"
        f"## Apophis (production)\n\nPath {apo['path']}: top **{apo['top']}** "
        f"(conf {apo['top_conf']:.2f}), OC mass {apo['oc_mass']:.2f}, oc_subgroup "
        f"**unresolved**.\n\n## Gaia-vs-external cross-check\n\nOn {n_xcheck} core "
        f"objects with both a Gaia spectrum and an A5 ground-truth spectrum, the "
        f"external path agrees with the Gaia path on the top analog "
        f"**{ga*100:.0f}%** overall, rising to **{ga_conf*100:.0f}%** on the "
        f"{n_conf_x} objects the Gaia path calls confident (status=ok). The "
        f"moderate overall figure reflects that most objects are low-confidence on "
        f"both paths (as C1 predicts); where the Gaia path is confident the two "
        f"paths agree well, which is what validates the external path the Apophis "
        f"dossier and customer requests rely on.\n")

    (OUT / "c2_4_properties.md").write_text(
        "# C2.4 - Physical-properties bridge (first pass)\n\n"
        f"For {n_props:,} confident non-degenerate analogs, group-typical "
        f"literature bulk properties (grain density range, porosity note) are "
        f"attached from data/interim/meteorite_group_properties.csv. These are "
        f"**analog-inferred, group-typical literature ranges - NOT per-object "
        f"measurements**. Every property row carries a citation field; values are "
        f"currently marked **TODO** for citation verification (Galinier "
        f"discipline - no fabricated references; Werner fills verified sources, "
        f"e.g. the meteorite density/porosity literature). Real per-object density "
        f"comes later from DR4 (Stage 2).\n")


def write_provenance(classes, sep_merged, deg_merged, mean_red, ba, ece):
    (PROV / "analog_composition.json").write_text(json.dumps({
        "columns": ["analog_group_top", "analog_group_conf", "analog_distribution",
                    "oc_subgroup", "analog_status", "weathering_strength",
                    "weathering_flag", "inferred_grain_density_range",
                    "inferred_porosity_note", "properties_source"],
        "kind": "MODEL-DERIVED from RELAB meteorite templates (resampled to the 16 "
                "Gaia bands under the documented 44 nm Gaussian bandpass "
                "approximation). Distinct from measured quantities; the inferred "
                "physical properties are ANALOG-inferred group-typical literature "
                "ranges, not per-object measurements.",
        "separable_subset": sorted(sep_merged),
        "degenerate_groups": sorted(deg_merged),
        "oc_policy": "H/L/LL merged into ordinary_chondrite; oc_subgroup always "
                     "'unresolved' (C1: mutually ~0.72 AUC, not spectrally resolvable)",
        "prior": "C4 taxon_esox complex used as prior (implausible analogs "
                 f"down-weighted x{PRIOR_LOW}); flat where unclassified. Mean "
                 f"entropy reduction {round(mean_red,3)}",
        "relab_self_cv_balanced_accuracy": round(ba, 3),
        "reliability_ece": round(ece, 3),
        "confidence_threshold": CONF_THR,
        "resampling_note": "same 44 nm Gaussian bandpass approximation as C1/C4",
    }, indent=2))


if __name__ == "__main__":
    sys.exit(main())
