"""C1 degeneration analysis - shared utilities.

Band grid, RELAB->Gaia resampling (guardrail 6), identical 550 nm normalization
(guardrail 3), realistic per-band Gaia noise (guardrail 2), viability filtering
(guardrail 4) and group_kind handling (guardrail 5). No modelling here - just
the data plumbing every stage shares.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
I = ROOT / "data" / "interim"
GASP_CAT = Path.home() / "gasp" / "data" / "final" / "gasp_catalog_v2.parquet"

# GASP / Gaia DR3 16-band grid (nm), 44 nm spacing, normalized at 550 nm.
BANDS = [374, 418, 462, 506, 550, 594, 638, 682, 726, 770, 814, 858, 902, 946, 990, 1034]
REF_NM = 550
BAND_FWHM_NM = 44.0          # effective resolution proxy = sample spacing
BAND_SIGMA = BAND_FWHM_NM / 2.3548
VIABLE_MIN_SPECTRA = 8
VIABLE_MIN_METEORITES = 3


# Bus-DeMeo / Mahlke class -> taxonomic complex (guardrail: complex vs class).
COMPLEX_MAP = {
    "S": "S", "Sa": "S", "Sq": "S", "Sr": "S", "Sv": "S", "Sw": "S", "Q": "S",
    "C": "C", "Cb": "C", "Cg": "C", "Cgh": "C", "Ch": "C", "B": "C",
    "X": "X", "Xc": "X", "Xe": "X", "Xk": "X", "Xn": "X", "M": "X", "P": "X",
    "Pk": "X", "Mk": "X", "E": "X",
    "V": "V", "D": "D", "A": "A",
    "K": "K/L", "L": "K/L", "Z": "other", "T": "other", "O": "other", "R": "other",
}


def to_complex(cls):
    if cls is None or (isinstance(cls, float)):
        return None
    return COMPLEX_MAP.get(str(cls).strip(), "other")


def group_kind_map():
    gm = pd.read_csv(I / "relab_group_map.csv")
    return dict(zip(gm["canonical"], gm["group_kind"]))


def gasp_band_errors():
    """Per-band absolute error samples from the GASP core (for noise injection).
    Returns dict band -> np.array of error magnitudes (guardrail 2)."""
    g = pd.read_parquet(GASP_CAT, columns=[f"err_{b}" for b in BANDS])
    out = {}
    for b in BANDS:
        e = g[f"err_{b}"].to_numpy()
        out[b] = e[np.isfinite(e) & (e > 0)]
    return out


def gaussian_resample(wl_nm, refl, centers=BANDS, sigma=BAND_SIGMA):
    """Integrate a high-res spectrum onto Gaia band centers with a Gaussian
    bandpass of FWHM = band spacing (guardrail 6). A band is NaN if the source
    spectrum has no sample within +/- sigma of the center (no extrapolation).

    METHOD NOTE: the published Gaia asteroid passband response curves are not
    in the local dataset; the Gaussian bandpass at the 44 nm sampling interval
    is an explicit, documented approximation of the effective resolution. If
    the true response becomes available it drops in here unchanged elsewhere.
    """
    wl_nm = np.asarray(wl_nm, float)
    refl = np.asarray(refl, float)
    order = np.argsort(wl_nm)
    wl_nm, refl = wl_nm[order], refl[order]
    out = np.full(len(centers), np.nan)
    for k, c in enumerate(centers):
        w = np.exp(-0.5 * ((wl_nm - c) / sigma) ** 2)
        near = np.abs(wl_nm - c) <= sigma
        if not near.any() or w.sum() <= 0:
            continue
        out[k] = float(np.sum(w * refl) / np.sum(w))
    return out


def normalize_at_ref(vec, centers=BANDS, ref=REF_NM):
    """Divide by the reference-band value (identical to GASP's 550 nm norm)."""
    ri = centers.index(ref)
    v = np.asarray(vec, float)
    if not np.isfinite(v[ri]) or v[ri] == 0:
        return np.full_like(v, np.nan)
    return v / v[ri]


def load_relab_resampled(require_full=True, centers=BANDS):
    """Viable primary-group RELAB spectra, resampled to the band grid and
    normalized at 550 nm. Returns (X, meta) where meta has relab_spectrum_id,
    meteorite_name, relab_sample_id, relab_group. Only primary group_kind
    (guardrail 5); viability >=8 spectra & >=3 meteorites (guardrail 4)."""
    S = pd.read_parquet(I / "relab_samples.parquet")
    SP = pd.read_parquet(I / "relab_spectra.parquet")
    kind = group_kind_map()
    S["group_kind"] = S["relab_group"].map(kind)
    prim = S[S["group_kind"] == "primary"].copy()
    agg = prim.groupby("relab_group").agg(
        n=("relab_spectrum_id", "size"),
        m=("meteorite_name", lambda x: x.str.lower().nunique()))
    viable = set(agg[(agg["n"] >= VIABLE_MIN_SPECTRA)
                     & (agg["m"] >= VIABLE_MIN_METEORITES)].index)
    prim = prim[prim["relab_group"].isin(viable)]

    SP = SP[SP["relab_spectrum_id"].isin(set(prim["relab_spectrum_id"]))].copy()
    SP["wl_nm"] = SP["wavelength_um"] * 1000.0
    rows, meta = [], []
    keep = prim.set_index("relab_spectrum_id")
    for sid, grp in SP.groupby("relab_spectrum_id"):
        vec = gaussian_resample(grp["wl_nm"].to_numpy(), grp["reflectance"].to_numpy(),
                                centers)
        vec = normalize_at_ref(vec, centers)
        if require_full and not np.isfinite(vec).all():
            continue
        r = keep.loc[sid]
        rows.append(vec)
        meta.append({"relab_spectrum_id": sid,
                     "meteorite_name": str(r["meteorite_name"]).lower().strip(),
                     "relab_sample_id": r["relab_sample_id"],
                     "relab_group": r["relab_group"]})
    X = np.vstack(rows)
    return X, pd.DataFrame(meta), sorted(viable)


def inject_gaia_noise(X, band_errors, rng, centers=BANDS, ref=REF_NM):
    """Add realistic per-band Gaussian noise drawn from the GASP error
    distribution (guardrail 2). The 550 nm anchor stays fixed (its error is
    definitional zero after normalization). Returns a noisy copy."""
    Xn = X.copy()
    ri = centers.index(ref)
    for k, b in enumerate(centers):
        if k == ri:
            continue
        sig = rng.choice(band_errors[b], size=len(X))
        Xn[:, k] = X[:, k] + rng.normal(0.0, 1.0, size=len(X)) * sig
    return Xn
