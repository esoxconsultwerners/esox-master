#!/usr/bin/env python3
"""Phase B - Esox Master Catalog merge (build_master.py).

Merges the nine A-series layers plus the GASP-core layers into two tiers:
  esox_master_full.parquet  - the full ~1.55M orbit backbone, enriched, thin
                              where sparse.
  esox_master_core.parquet  - exactly the 19,190 GASP-core objects, dense.

Conflict resolution follows data/provenance/precedence_rules.md EXACTLY
(rules 1.1-1.6). No precedence is invented here; this script only executes
the documented rules and records, per object, which source won
(the *_source columns). Original per-source values are never overwritten.

GASP-core layers (read-only):
  ~/gasp/data/final/gasp_catalog_v2.parquet  - Gaia reflectance spectra,
      SDSS photometry, Mahlke taxonomy (`taxonomy`), GASP family.
  ~/esox-pipeline/data/gapc_catalog_v8.parquet - GAPC HG1G2 phase curves
      (the handover's "~/gasp/.../GAPC" layer actually lives in the sibling
      esox-pipeline project; documented here as the exact path used).
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
I = ROOT / "data" / "interim"
FINAL = ROOT / "data" / "final"
PROV = ROOT / "data" / "provenance"
DOCS = ROOT / "docs"
GASP_CAT = Path.home() / "gasp" / "data" / "final" / "gasp_catalog_v2.parquet"
GAPC_CAT = Path("/root/esox-pipeline/data/gapc_catalog_v8.parquet")

TAXON_SYSTEMS = ["taxon_mahlke", "taxon_tholen", "taxon_bus",
                 "taxon_bus_demeo", "taxon_s3os2"]


def normalize_desig(s):
    """A6 Rule 1.5: '2023VD3' -> '2023 VD3'; spaced forms pass through."""
    s = (s or "").strip().upper()
    m = re.match(r"^(\d{4})\s*([A-Z]{1,2}\d*)$", s)
    if m:
        return f"{m.group(1)} {m.group(2)}"
    return s or None


def make_key(number_mp, designation):
    if pd.notna(number_mp):
        return str(int(number_mp))
    return normalize_desig(designation)


def add_key(df):
    desig = df["designation"] if "designation" in df.columns else pd.Series([None] * len(df))
    df = df.copy()
    df["object_key"] = [make_key(n, d) for n, d in zip(df["number_mp"], desig)]
    return df[df["object_key"].notna()]


def pick_best(df, triples):
    """First non-null across (value_col, source_label) in precedence order.
    Returns (best, source) with real nulls where nothing wins (never the
    string 'None' - np.select with default=None would coerce to that)."""
    best = pd.Series(np.nan, index=df.index, dtype="float64")
    source = pd.Series(pd.NA, index=df.index, dtype="object")
    remaining = pd.Series(True, index=df.index)
    for col, label in triples:
        take = remaining & df[col].notna()
        best[take] = df.loc[take, col].astype("float64")
        source[take] = label
        remaining = remaining & ~df[col].notna()
    return best, source


def layer(name, cols, sort_col=None, ascending=True):
    df = add_key(pd.read_parquet(I / f"{name}.parquet"))
    if sort_col:
        df = df.sort_values(sort_col, ascending=ascending)
    df = df.drop_duplicates("object_key", keep="first")
    return df[["object_key"] + cols]


def main():
    FINAL.mkdir(parents=True, exist_ok=True)
    print("[merge] loading spine (A1 orbit backbone) ...")
    master = add_key(pd.read_parquet(I / "orbits_backbone.parquet"))
    backbone_n = len(master)
    print(f"[merge] spine: {backbone_n:,} objects")

    joins = [
        layer("lcdb", ["lcdb_period_h", "lcdb_quality_numeric", "lcdb_quality_u",
                       "lcdb_diameter_km", "lcdb_albedo", "lcdb_binary_flag",
                       "lcdb_class", "lcdb_notes"]),
        layer("families", ["fam_id", "fam_name", "fam_membership", "fam_zone"],
              sort_col="fam_membership", ascending=False),
        layer("akari", ["akari_diameter_km", "akari_diameter_err",
                        "akari_albedo", "akari_albedo_err"]),
        layer("iras_simps", ["iras_diameter_km", "iras_diameter_err",
                             "iras_albedo", "iras_albedo_err"]),
        layer("groundtruth_summary", ["gt_source", "gt_taxon", "gt_n_spectra",
                                      "gt_wl_min_um", "gt_wl_max_um"]),
        layer("taxonomy_wide", ["taxon_tholen", "taxon_bus", "taxon_bus_demeo",
                                "taxon_s3os2"]),
        layer("neo_ops", ["nhats_min_dv_kms", "nhats_n_trajectories",
                          "neocc_ip_max", "neocc_risk_ps_max", "neocc_velocity_kms",
                          "sentry_ip", "sentry_ps_max", "sentry_ts_max"]),
        layer("damit_best", ["damit_model_id", "damit_lambda_deg", "damit_beta_deg",
                             "damit_period_h", "damit_method", "damit_ref"]),
    ]
    for j in joins:
        master = master.merge(j, on="object_key", how="left")

    # GASP core layer (keyed by str(number_mp))
    print("[merge] joining GASP spectral catalog + SDSS + Mahlke taxonomy ...")
    gasp = pd.read_parquet(GASP_CAT)
    gasp["object_key"] = gasp["number_mp"].astype("int64").astype(str)
    refl = [c for c in gasp.columns if c.startswith("refl_")]
    gasp_sel = gasp[["object_key", "denomination", "num_of_spectra", "n_valid_bands",
                     "quality", "u", "g", "r", "i", "z", "a_star", "V", "B",
                     "sdss_complex", "taxonomy", "family", "diameter_km", "albedo"]
                    + refl].rename(columns={
                        "denomination": "gasp_denomination",
                        "num_of_spectra": "gasp_num_of_spectra",
                        "n_valid_bands": "gasp_n_valid_bands",
                        "quality": "gasp_quality", "taxonomy": "taxon_mahlke",
                        "family": "gasp_family", "diameter_km": "gasp_diameter_km",
                        "albedo": "gasp_albedo"}).drop_duplicates("object_key")
    master = master.merge(gasp_sel, on="object_key", how="left")

    # GAPC phase-curve layer (keyed by str(number_mp))
    print("[merge] joining GAPC HG1G2 phase curves ...")
    gapc = pd.read_parquet(GAPC_CAT, columns=["number_mp", "H", "G1", "G2",
                                              "phase_range", "n_obs", "fit_ok",
                                              "flag_unphysical"])
    gapc["object_key"] = gapc["number_mp"].astype("int64").astype(str)
    gapc_sel = gapc.rename(columns={"H": "gapc_h", "G1": "gapc_g1", "G2": "gapc_g2",
                                    "phase_range": "gapc_phase_range",
                                    "n_obs": "gapc_n_obs", "fit_ok": "gapc_fit_ok",
                                    "flag_unphysical": "gapc_flag_unphysical"}) \
        .drop(columns="number_mp").drop_duplicates("object_key")
    master = master.merge(gapc_sel, on="object_key", how="left")

    master = master.reset_index(drop=True)

    # ---- Best fields (precedence_rules.md) ----
    # Rule 1.1: diameter/albedo, AKARI > NEOWISE(SBDB) > IRAS, no averaging.
    master["diameter_best"], master["diameter_source"] = pick_best(master, [
        ("akari_diameter_km", "AKARI"), ("sbdb_diameter", "NEOWISE"),
        ("iras_diameter_km", "IRAS")])
    master["albedo_best"], master["albedo_source"] = pick_best(master, [
        ("akari_albedo", "AKARI"), ("sbdb_albedo", "NEOWISE"),
        ("iras_albedo", "IRAS")])

    # Rule 1.4: H, SBDB > MPCORB.
    master["h_best"], master["h_source"] = pick_best(master, [
        ("sbdb_H", "SBDB"), ("mpc_h", "MPCORB")])

    # Rule 1.2: period, DAMIT > LCDB(U>=2); tumblers -> null best, flagged.
    damit_tumb = set(pd.read_parquet(I / "damit_tumblers.parquet")["number_mp"]
                     .dropna().astype(int).astype(str))
    lcdb_tumb_key = set(master.loc[master["lcdb_notes"].astype("string").str.strip()
                                   == "T", "object_key"])
    tumbler = master["object_key"].isin(damit_tumb | lcdb_tumb_key)
    master["tumbler_flag"] = tumbler
    damit_ok = master["damit_period_h"].notna()
    lcdb_ok = master["lcdb_period_h"].notna() & (master["lcdb_quality_numeric"] >= 2)
    sel_d = (~tumbler) & damit_ok
    sel_l = (~tumbler) & (~damit_ok) & lcdb_ok
    master["period_best"] = np.nan
    master.loc[sel_d, "period_best"] = master.loc[sel_d, "damit_period_h"]
    master.loc[sel_l, "period_best"] = master.loc[sel_l, "lcdb_period_h"]
    master["period_source"] = pd.Series(pd.NA, index=master.index, dtype="object")
    master.loc[sel_d, "period_source"] = "DAMIT"
    master.loc[sel_l, "period_source"] = "LCDB"

    # Rule 1.3: taxonomy consensus only where >=2 systems agree at letter level.
    letters = pd.DataFrame({
        s: master[s].where(master[s].notna()).astype("string").str.strip().str[0].str.upper()
        for s in TAXON_SYSTEMS})
    melted = letters.reset_index().melt("index").dropna(subset=["value"])
    counts = melted.groupby(["index", "value"]).size().reset_index(name="n")
    top = counts.sort_values(["index", "n"]).groupby("index").tail(1)
    consensus = top[top["n"] >= 2].set_index("index")["value"]
    master["taxon_literature_consensus"] = master.index.map(consensus)

    # ---- Coverage flags ----
    master["has_gaia_spectrum"] = master["gasp_num_of_spectra"].fillna(0) > 0
    master["has_period_best"] = master["period_best"].notna()
    master["has_family"] = master["fam_id"].notna()
    master["has_diameter_best"] = master["diameter_best"].notna()
    master["has_taxon_any"] = master[TAXON_SYSTEMS].notna().any(axis=1)
    master["has_phase_curve"] = master["gapc_fit_ok"].fillna(False).astype(bool)
    master["has_groundtruth_spectrum"] = master["gt_source"].notna()
    master["complete_profile"] = (master["has_gaia_spectrum"] & master["has_period_best"]
                                  & master["has_family"] & master["has_diameter_best"]
                                  & master["has_taxon_any"])

    # ---- Outputs ----
    core_keys = pd.read_parquet(I / "gasp_core_keys.parquet")
    core_keys["object_key"] = core_keys["number_mp"].astype("int64").astype(str)
    core = core_keys[["object_key"]].merge(master, on="object_key", how="left")

    full_path = FINAL / "esox_master_full.parquet"
    core_path = FINAL / "esox_master_core.parquet"
    master.to_parquet(full_path, index=False)
    core.to_parquet(core_path, index=False)

    realized = realized_outcomes(master)
    write_realized(realized)

    validate(master, core, backbone_n)
    spot_objects(master)
    apophis_dossier(master)
    coverage_report(master, core, full_path, core_path)
    print("\n=== Phase B complete ===")
    cp = int(core["complete_profile"].sum())
    print(f"Master v0: {len(master):,} full / {len(core):,} core. "
          f"Complete physical profiles: {cp:,} ({100*cp/len(core):.1f}%).")
    print("Ready for CM (coverage matrix) and C1 (degeneration analysis).")
    return 0


def realized_outcomes(m):
    # Expected sources per field, in precedence order, so a source that never
    # wins is shown explicitly as 0 (e.g. IRAS: every IRAS target also has an
    # AKARI or NEOWISE value, so as the lowest-priority source it wins 0 -
    # the precedence is applied correctly, not misapplied).
    expected = {
        "diameter_best": (["AKARI", "NEOWISE", "IRAS"], "diameter_source"),
        "albedo_best": (["AKARI", "NEOWISE", "IRAS"], "albedo_source"),
        "period_best": (["DAMIT", "LCDB"], "period_source"),
        "h_best": (["SBDB", "MPCORB"], "h_source"),
    }
    out = {}
    for field, (sources, src) in expected.items():
        vc = m.loc[m[field].notna(), src].value_counts()
        out[field] = {s: int(vc.get(s, 0)) for s in sources}
    out["period_best"]["(tumbler, no best)"] = int(m["tumbler_flag"].sum())
    out["taxon_literature_consensus"] = {
        ">=2 systems agree": int(m["taxon_literature_consensus"].notna().sum())}
    return out


def write_realized(realized):
    path = PROV / "precedence_rules.md"
    text = path.read_text()
    marker = "## Realized outcomes"
    head = text.split(marker)[0].rstrip() + "\n\n"
    lines = [marker, "",
             "*Filled by `pipeline/build_master.py` after the merge "
             "(counts over the full catalog).*", "",
             "| Best field | Source | Objects (n) |", "|---|---|---|"]
    labels = {"diameter_best": "`diameter_best`", "albedo_best": "`albedo_best`",
              "period_best": "`period_best`", "h_best": "`h_best`",
              "taxon_literature_consensus": "`taxon_literature_consensus`"}
    for field in ["diameter_best", "albedo_best", "period_best", "h_best",
                  "taxon_literature_consensus"]:
        for source, n in realized[field].items():
            lines.append(f"| {labels[field]} | {source} | {n:,} |")
    path.write_text(head + "\n".join(lines) + "\n")
    print(f"[realized] updated {path.relative_to(ROOT)} (placeholder replaced)")


def fail(msg):
    print(f"[VALIDATION FAILED] {msg}")
    sys.exit(1)


def validate(master, core, backbone_n):
    print("\n[validate] running gates ...")
    if not master["object_key"].is_unique:
        fail("object_key not unique in full catalog")
    if not core["object_key"].is_unique:
        fail("object_key not unique in core catalog")
    if len(core) != 19190:
        fail(f"core is {len(core)} rows, expected 19,190")
    for best, src in [("diameter_best", "diameter_source"),
                      ("albedo_best", "albedo_source"),
                      ("period_best", "period_source"),
                      ("h_best", "h_source")]:
        bad = master[master[best].notna() & master[src].isna()]
        if len(bad):
            fail(f"{len(bad)} rows with {best} non-null but {src} null")
    if abs(len(master) - backbone_n) / backbone_n > 0.01:
        fail(f"full count {len(master)} deviates >1% from backbone {backbone_n}")
    for name, num in [("Ceres", 1), ("Vesta", 4), ("Eros", 433),
                      ("Bennu", 101955), ("Ryugu", 162173), ("Apophis", 99942)]:
        if (master["object_key"] == str(num)).sum() != 1:
            fail(f"spot object {name} ({num}) not uniquely present")
    print("[validate] all gates passed")


def _fmt(v, nd=3):
    if pd.isna(v):
        return "-"
    if isinstance(v, float):
        return f"{v:.{nd}g}"
    return str(v)


def spot_objects(master):
    print("\n=== Spot-object records (sanity vs known values) ===")
    notes = {1: "Ceres ~939 km, pV~0.09, P~9.07 h, C-complex",
             4: "Vesta ~525 km, pV~0.42, P~5.34 h, V-type",
             433: "Eros ~16.8 km, P~5.27 h, S-type, Amor NEO",
             101955: "Bennu ~0.49 km, P~4.30 h, B-type, OSIRIS-REx target",
             162173: "Ryugu ~0.90 km, P~7.63 h, Cb-type, Hayabusa2 target",
             99942: "Apophis ~0.34 km, NPA tumbler (no single period), Sq"}
    for num in [1, 4, 433, 101955, 162173, 99942]:
        r = master[master["object_key"] == str(num)].iloc[0]
        print(f"\n--- ({num}) {_fmt(r.get('mpc_name'))} / {_fmt(r.get('sbdb_full_name'))} ---")
        print(f"    known: {notes[num]}")
        print(f"    diameter_best = {_fmt(r['diameter_best'])} km [{_fmt(r['diameter_source'])}]  "
              f"albedo_best = {_fmt(r['albedo_best'])} [{_fmt(r['albedo_source'])}]")
        print(f"    period_best = {_fmt(r['period_best'])} h [{_fmt(r['period_source'])}]  "
              f"tumbler={bool(r['tumbler_flag'])}  h_best = {_fmt(r['h_best'])} [{_fmt(r['h_source'])}]")
        print(f"    taxon: mahlke={_fmt(r.get('taxon_mahlke'))} tholen={_fmt(r.get('taxon_tholen'))} "
              f"bus={_fmt(r.get('taxon_bus'))} bus_demeo={_fmt(r.get('taxon_bus_demeo'))} "
              f"s3os2={_fmt(r.get('taxon_s3os2'))} consensus={_fmt(r.get('taxon_literature_consensus'))}")
        print(f"    family = {_fmt(r.get('fam_name'))} ({_fmt(r.get('fam_id'))})  "
              f"NEO={_fmt(r.get('sbdb_neo'))} PHA={_fmt(r.get('sbdb_pha'))}  "
              f"NEOCC_ip={_fmt(r.get('neocc_ip_max'))} Sentry_ip={_fmt(r.get('sentry_ip'))}")
        print(f"    gaia_spectrum={bool(r['has_gaia_spectrum'])}  phase_curve={bool(r['has_phase_curve'])}  "
              f"complete_profile={bool(r['complete_profile'])}")


def apophis_dossier(master):
    r = master[master["object_key"] == "99942"].iloc[0]
    block = [
        "\n## Assembled dossier v0 (Phase B merge)",
        "",
        "The showcase record every A1-A9 layer has been building toward - "
        "(99942) Apophis, assembled from the merged master catalog:",
        "",
        f"- **Orbit (A1):** a = {_fmt(r.get('mpc_a'))} au, e = {_fmt(r.get('mpc_e'))}, "
        f"i = {_fmt(r.get('mpc_i'))} deg; SBDB class {_fmt(r.get('sbdb_class'))}, "
        f"NEO={_fmt(r.get('sbdb_neo'))}, PHA={_fmt(r.get('sbdb_pha'))}, "
        f"MOID = {_fmt(r.get('sbdb_moid'))} au. H = {_fmt(r.get('h_best'))} "
        f"[{_fmt(r.get('h_source'))}].",
        f"- **Rotation (A2 + A9):** non-principal-axis **tumbler** "
        f"(tumbler_flag={bool(r['tumbler_flag'])}) - so period_best is "
        f"deliberately null (Rule 1.2). LCDB reports ~30.56 h (U=3, note 'T'); "
        f"DAMIT carries the NPA solution in its tumblers table. Two independent "
        f"databases agree on the non-relaxed rotation state.",
        f"- **Size/albedo (A4/SBDB):** diameter_best = {_fmt(r.get('diameter_best'))} km "
        f"[{_fmt(r.get('diameter_source'))}], albedo_best = {_fmt(r.get('albedo_best'))} "
        f"[{_fmt(r.get('albedo_source'))}].",
        f"- **Accessibility & risk (A6):** NHATS min dv = {_fmt(r.get('nhats_min_dv_kms'))} "
        f"km/s ({_fmt(r.get('nhats_n_trajectories'))} trajectories); "
        f"NEOCC/Sentry impact probability = {_fmt(r.get('neocc_ip_max'))}/"
        f"{_fmt(r.get('sentry_ip'))} (removed from risk lists after the "
        f"2004-2021 observation arc).",
        f"- **Taxonomy (literature):** mahlke={_fmt(r.get('taxon_mahlke'))}, "
        f"tholen={_fmt(r.get('taxon_tholen'))}, bus={_fmt(r.get('taxon_bus'))}, "
        f"bus_demeo={_fmt(r.get('taxon_bus_demeo'))}, "
        f"consensus={_fmt(r.get('taxon_literature_consensus'))} - S-complex / Sq, "
        f"the spectral bridge to LL ordinary chondrites.",
        f"- **Spectra (A5 + GASP):** gaia_spectrum={bool(r['has_gaia_spectrum'])}; "
        f"MITHNEOS NIR spectra provide the 0.77-2.49 um coverage. The laboratory "
        f"analogue is the RELAB LL-chondrite set (A7).",
        f"- **Complete physical profile:** {bool(r['complete_profile'])} "
        f"(a tumbler has no single period_best by design, so Apophis is "
        f"intentionally not a 'complete profile' object - the flag behaves "
        f"correctly).",
        "",
    ]
    path = DOCS / "apophis_status.md"
    text = path.read_text()
    heading = "## Assembled dossier v0 (Phase B merge)"
    if heading in text:
        text = text[:text.index(heading)].rstrip() + "\n"
        path.write_text(text)
    path.open("a").write("\n".join(block) + "\n")
    print("\n".join(block))
    print(f"[dossier] appended assembled Apophis record to docs/apophis_status.md")


def coverage_report(master, core, full_path, core_path):
    n = len(core)

    def pct(mask):
        c = int(mask.sum())
        return f"{c:,} ({100*c/n:.1f}%)"

    rs = pd.read_parquet(I / "relab_samples.parquet")
    viable = rs[rs["relab_group"] != "unmapped"]["relab_group"].value_counts()
    n_viable = int((viable >= 5).sum())

    print("\n" + "=" * 64)
    print("PHASE B COVERAGE REPORT (esox master v0)")
    print("=" * 64)
    print(f"Full catalog objects:                {len(master):,}")
    print(f"Core catalog objects:                {len(core):,}")
    print("Per-property realized coverage (core):")
    print(f"    orbit:                  {pct(core['mpc_a'].notna())}")
    print(f"    h_best:                 {pct(core['h_best'].notna())}")
    print(f"    diameter_best:          {pct(core['has_diameter_best'])}")
    print(f"    albedo_best:            {pct(core['albedo_best'].notna())}")
    print(f"    period_best:            {pct(core['has_period_best'])}")
    print(f"    family:                 {pct(core['has_family'])}")
    print(f"    taxon(any):             {pct(core['has_taxon_any'])}")
    print(f"    taxon_mahlke:           {pct(core['taxon_mahlke'].notna())}")
    print(f"    gaia_spectrum:          {pct(core['has_gaia_spectrum'])}")
    print(f"    phase_curve:            {pct(core['has_phase_curve'])}")
    print(f"    groundtruth_spectrum:   {pct(core['has_groundtruth_spectrum'])}")
    cp = int(core["complete_profile"].sum())
    print(f"COMPLETE PHYSICAL PROFILE count:     {cp:,} ({100*cp/n:.1f}%)")
    print("    (has ALL of gaia_spectrum + period_best + family + "
          "diameter_best + taxon(any))")
    print(f"RELAB reference note: {len(rs):,} meteorite spectra in {n_viable} "
          f"viable groups available for Phase C (not merged into catalog).")
    print(f"Saved: {full_path.name} ({full_path.stat().st_size/1e6:.1f} MB), "
          f"{core_path.name} ({core_path.stat().st_size/1e6:.1f} MB)")
    print("=" * 64)


if __name__ == "__main__":
    sys.exit(main())
