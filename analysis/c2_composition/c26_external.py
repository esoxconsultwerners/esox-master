"""C2.6 - migrated external-spectrum composition path.

Retires the RandomForest external matcher (a template-density catch-all, see
diag/C25_MATCHER_FINDING.md) in favour of the C2.5 per-group aggregate-likelihood
matcher + RELAB-manifold gate, for BOTH the full-band and NIR external paths.
NIR ground-truth spectra were measured to be 93.7% inside the RELAB manifold
(nir/c26_1_nir_manifold.md), so confident NIR analogs are defensible; the gate
still prunes any out-of-manifold spectrum exactly as it does the Gaia core - no
special pleading for the external path.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "analysis" / "c1_degeneration"))
sys.path.insert(0, str(HERE))
import c1_common as cc
import c2_composition as c2
import c25_matcher as m5

NIR = c2.NIR_BANDS
NREF = c2.NIR_REF
_NRI = NIR.index(NREF)
_NFEAT = [b for b in NIR if b != NREF]
_NFIDX = [NIR.index(b) for b in _NFEAT]
OC_LABEL = c2.OC_LABEL
CONF_THR = 0.55

_STATE = None


def nir_features(V7):
    V = np.atleast_2d(np.asarray(V7, float))
    ref = V[:, _NRI]
    return np.log(np.clip((V / ref[:, None])[:, _NFIDX], 1e-4, None))


def _fit_with_q95(F, y, groups):
    mm = m5.GroupGaussianMatcher().fit(F, y)
    mm.fit_calibration(F, y, groups)
    ci = {c: i for i, c in enumerate(mm.classes_)}
    d2 = np.array([(F[i] - mm.mu_[ci[y[i]]]) @ mm.prec_ @ (F[i] - mm.mu_[ci[y[i]]])
                   for i in range(len(y))])
    return mm, float(np.quantile(d2, 0.95))


def _build():
    global _STATE
    if _STATE is not None:
        return _STATE
    _, meta, _ = cc.load_relab_resampled()
    V16 = c2.relab_features(meta, cc.BANDS, cc.REF_NM)
    ok = np.isfinite(V16).all(1)
    yf = meta[ok]["relab_group"].map(c2.merge_label).to_numpy()
    gf = meta[ok]["meteorite_name"].to_numpy()
    full, q95f = _fit_with_q95(m5.features_from_16(V16[ok]), yf, gf)
    V7 = c2.relab_features(meta, NIR, NREF)
    okn = np.isfinite(V7).all(1)
    yn = meta[okn]["relab_group"].map(c2.merge_label).to_numpy()
    gn = meta[okn]["meteorite_name"].to_numpy()
    nir, q95n = _fit_with_q95(nir_features(V7[okn]), yn, gn)
    sep, deg, sepm = c2.separable_sets()
    _STATE = dict(full=full, q95f=q95f, nir=nir, q95n=q95n, sepm=sepm,
                  degm={c2.merge_label(x) for x in deg})
    return _STATE


def match_external_spectrum(wl_nm, refl, complex_label=None):
    st = _build()
    wl_nm = np.asarray(wl_nm, float)
    refl = np.asarray(refl, float)
    v16 = cc.gaussian_resample(wl_nm, refl, cc.BANDS)
    ri = cc.BANDS.index(cc.REF_NM)
    if np.isfinite(v16[ri]) and v16[ri] and np.isfinite(v16).sum() >= 10:
        mm, q95, F1, path = st["full"], st["q95f"], m5.features_from_16(v16)[0], "full-13band"
    else:
        vn = cc.gaussian_resample(wl_nm, refl, NIR)
        mm, q95, F1, path = st["nir"], st["q95n"], nir_features(vn)[0], "nir-6band"
    raw = mm.predict_proba(F1[None, :])[0]
    post, _ = c2.apply_prior(raw, mm.classes_, complex_label)
    ti = int(post.argmax())
    top = mm.classes_[ti]
    d = F1 - mm.mu_[ti]
    d2 = float(d @ mm.prec_ @ d)
    within = bool(d2 <= q95)
    conf = float(mm.calibrate([float(post.max())])[0])
    if top in st["degm"]:
        status, confident = "degenerate", "degenerate"
    elif within and conf >= CONF_THR:
        status, confident = "ok", top
    else:
        status, confident = "unresolved", "unresolved"
    dist = {str(c): float(p) for c, p in zip(mm.classes_, post) if c in st["sepm"]}
    s = sum(dist.values())
    if s > 0:
        dist = {c: round(v / s, 4) for c, v in dist.items()}
    oc_mass = float(sum(p for c, p in zip(mm.classes_, post) if c == OC_LABEL))
    return {"path": path, "top": top, "top_conf": round(conf, 4),
            "oc_mass": round(oc_mass, 4), "manifold_d2": round(d2, 3),
            "manifold_q95": round(q95, 3), "within_manifold": within,
            "analog_status": status, "confident_analog": confident,
            "oc_subgroup": "unresolved" if top == OC_LABEL else None,
            "matcher": "c2.5-per-group", "distribution": dist}
