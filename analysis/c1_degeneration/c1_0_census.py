#!/usr/bin/env python3
"""C1.0 - Census & band grid. Cheap first stage with the C1.0 decision gate."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c1_common as cc

OUT = Path(__file__).resolve().parent
GASP = pd.read_parquet(cc.GASP_CAT)


def relab_viable():
    S = pd.read_parquet(cc.I / "relab_samples.parquet")
    S["group_kind"] = S["relab_group"].map(cc.group_kind_map())
    rows = []
    for kind in ("primary", "soft", "coarse"):
        sub = S[S["group_kind"] == kind]
        agg = sub.groupby("relab_group").agg(
            n_spec=("relab_spectrum_id", "size"),
            n_met=("meteorite_name", lambda x: x.str.lower().nunique()))
        agg["group_kind"] = kind
        agg["viable"] = (agg["n_spec"] >= cc.VIABLE_MIN_SPECTRA) & (agg["n_met"] >= cc.VIABLE_MIN_METEORITES)
        rows.append(agg)
    return pd.concat(rows).reset_index()


def taxonomy_census():
    mah = GASP[["number_mp", "taxonomy"]].dropna(subset=["taxonomy"])
    mah_cls = mah["taxonomy"].value_counts()
    mah_cx = mah["taxonomy"].map(cc.to_complex).value_counts()
    tw = pd.read_parquet(cc.I / "taxonomy_wide.parquet")
    core = set(GASP["number_mp"].astype(int))
    twc = tw[tw["number_mp"].isin(core)]
    pds = {c: int(twc[c].notna().sum()) for c in
           ["taxon_tholen", "taxon_bus", "taxon_bus_demeo", "taxon_s3os2"]}
    return mah, mah_cls, mah_cx, pds


def main():
    va = relab_viable()
    prim = va[va["group_kind"] == "primary"]
    prim_viable = prim[prim["viable"]].sort_values("n_spec", ascending=False)
    n_prim_viable = len(prim_viable)

    mah, mah_cls, mah_cx, pds = taxonomy_census()
    n_labeled = len(mah)

    gate_pass = (n_prim_viable >= 6) and (n_labeled >= 150)

    lines = ["# C1.0 - Census & band grid", "",
             "## GASP / Gaia DR3 16-band grid",
             "", f"Centers (nm): {cc.BANDS}", f"Spacing: 44 nm uniform. "
             f"Normalization reference: **{cc.REF_NM} nm** (refl_550 == 1.0 for all "
             "GASP spectra; RELAB is normalized identically in every stage).", "",
             "## RELAB primary-group viability (>=8 spectra AND >=3 meteorites)",
             "", f"**{n_prim_viable} primary groups survive** the viability filter.", "",
             "| group | n_spectra | n_meteorites |", "|---|---:|---:|"]
    for _, r in prim_viable.iterrows():
        lines.append(f"| {r['relab_group']} | {int(r['n_spec'])} | {int(r['n_met'])} |")
    nonv = prim[~prim["viable"]]
    lines += ["", f"Primary groups excluded (below threshold): "
              f"{', '.join(f'{r.relab_group}({int(r.n_spec)})' for r in nonv.itertuples())}",
              "", "soft groups (L-LL, H-L) and coarse groups (C-ungrouped, iron, "
              "OC-ung, E) are reported separately in later stages and never counted "
              "as clean separability wins (guardrail 5).", ""]

    lines += ["## Taxonomy arm training set (GASP core, 19,190 objects)", "",
              f"GASP-Mahlke labels in core: **{n_labeled}** objects.", "",
              "By complex:", "",
              "| complex | n |", "|---|---:|"]
    for cx, n in mah_cx.items():
        lines.append(f"| {cx} | {int(n)} |")
    lines += ["", "By class (top):", "",
              "| class | n |", "|---|---:|"]
    for cl, n in mah_cls.head(18).items():
        lines.append(f"| {cl} | {int(n)} |")
    lines += ["", "Independent validation labels in core (PDS, A5b):",
              "", f"- Bus-DeMeo: **{pds['taxon_bus_demeo']}**  |  Bus: {pds['taxon_bus']}  "
              f"|  Tholen: {pds['taxon_tholen']}  |  S3OS2: {pds['taxon_s3os2']}",
              "", "(The ~371-object PDS Bus-DeMeo set is mostly outside the GASP "
              f"spectral core; only {pds['taxon_bus_demeo']} overlap and serve as the "
              "held-out validation set.)", ""]

    lines += ["## Decision gate C1.0", "",
              f"- primary RELAB groups viable: {n_prim_viable} (need >= 6) -> "
              f"{'PASS' if n_prim_viable >= 6 else 'FAIL'}",
              f"- labeled GASP-core objects (Mahlke): {n_labeled} (need >= ~150) -> "
              f"{'PASS' if n_labeled >= 150 else 'FAIL'}", "",
              f"**GATE {'PASSED - proceed to C1.1' if gate_pass else 'FAILED - STOP and renegotiate scope'}.**",
              ""]
    (OUT / "c1_0_census.md").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("\n[c1.0] gate", "PASS" if gate_pass else "FAIL")
    return 0 if gate_pass else 2


if __name__ == "__main__":
    sys.exit(main())
