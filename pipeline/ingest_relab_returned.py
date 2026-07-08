#!/usr/bin/env python3
"""AP_A7b: RELAB returned-asteroid-sample reference set.

AP_A7's meteorite filter correctly excluded the "Returned Asteroid Sample"
spectra - they are not meteorites. But returned asteroid material is the only
laboratory-measured, genuinely space-weathered asteroid surface in existence,
and is the ground truth for the Phase C2 space-weathering model. This package
recovers it as a separate, clearly labelled reference set - NOT part of the
meteorite library, NOT merged into the asteroid catalog.

Reuses the already-cached RELAB PDS4 labels under data/raw/relab/labels/
(no re-crawl; the label set is read straight off disk). Only the 31 data
tables for the kept spectra are downloaded.

What is actually in this bundle (verified by parsing all 23,606 cached
labels and by raw grep of every label XML):
  - specimen_type == "Returned Asteroid Sample": 31 spectra, ALL Ryugu
    (sample C0002, Hayabusa2 / JAXA; collection location "Asteroid Ryugu").
  - Itokawa: absent from this PDS4 bundle (0 mentions anywhere).
  - Bennu: absent (0 mentions).
  - 1,510 "Returned Lunar Sample" spectra exist but are DELIBERATELY EXCLUDED
    here: the Moon is not an asteroid, so Apollo/Luna material is out of scope
    for an asteroid space-weathering reference set (documented choice).

Quality rules follow A7 (reflectance QA, no mis-scaled runs), but the
grain-size-documented filter is NOT applied - these samples are precious
regardless of grain-size metadata; grain size is flagged where the label
gives it (structured specimen_min_size / specimen_max_size fields). Archive
values stay pure: no resampling/normalization, only nm->um unit conversion.
"""

import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ingest_relab import (BASE, RAW, INTERIM, PROVENANCE, fetch, parse_tab,
                          texture_of, sha256_file)

# body substring -> (parent_body, mission). Only Ryugu is present in this
# bundle; the others are kept so a future bundle refresh classifies them.
BODY_MISSION = {
    "ryugu": ("Ryugu", "Hayabusa2"),
    "itokawa": ("Itokawa", "Hayabusa"),
    "bennu": ("Bennu", "OSIRIS-REx"),
}


def txt(el):
    return (el.text or "").strip() if el is not None else ""


def parse_returned_label(path):
    """Richer parse than A7's: also pulls the structured grain-size,
    collection-location and provenance fields present in these labels.
    """
    tags = {}
    for el in ET.parse(path).getroot().iter():
        tags.setdefault(el.tag.split("}")[-1], []).append((el.text or "").strip())

    def one(name, default=""):
        v = [x for x in tags.get(name, []) if x != ""]
        return v[0] if v else (tags.get(name, [default])[0] if tags.get(name) else default)

    def num(name):
        v = one(name, "")
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan

    unit = (one("spectral_range_unit_name") or "nm").lower()
    gmin = num("specimen_min_size")
    gmax = num("specimen_max_size")
    return {
        "relab_spectrum_id": path.stem,
        "subdir": path.parent.name,
        "relab_sample_id": one("specimen_id"),
        "sample_name": one("specimen_name"),
        "description": one("specimen_description"),
        "specimen_type": one("specimen_type"),
        "grain_size_min_um": gmin,
        "grain_size_max_um": gmax,
        "relab_texture": texture_of("|".join(tags.get("material_subtype", []))),
        "measurement_type_raw": one("measurement_type"),
        "geometry": one("measurement_geometry_type"),
        "instrument": one("instrument_name"),
        "collection_location": one("specimen_collection_location"),
        "owner_name": one("specimen_owner_name"),
        "measurement_date": one("measurement_date_time"),
        "range_unit": unit,
    }


def body_mission_of(row):
    blob = f"{row['sample_name']} {row['collection_location']}".lower()
    for key, (body, mission) in BODY_MISSION.items():
        if key in blob:
            return body, mission
    return "unknown", "unknown"


def main():
    labels_dir = RAW / "labels"
    if not labels_dir.exists():
        print(f"[a7b] {labels_dir} missing - AP_A7 must have cached labels first")
        return 1
    labs = sorted(labels_dir.glob("*/*.xml"))
    print(f"[a7b] reading {len(labs):,} cached labels (no re-crawl) ...")
    meta = pd.DataFrame([parse_returned_label(p) for p in labs])

    ras = meta[meta["specimen_type"] == "Returned Asteroid Sample"].copy()
    print(f"[a7b] Returned Asteroid Sample spectra: {len(ras)}")
    bm = ras.apply(body_mission_of, axis=1)
    ras["parent_body"] = [b for b, _ in bm]
    ras["mission"] = [m for _, m in bm]

    session = requests.Session()
    print(f"[a7b] downloading {len(ras)} data tables ...")
    spectra_rows, npts, wlmin, wlmax = [], [], [], []
    for row in ras.itertuples():
        rel = f"data_reflectance/{row.subdir}/{row.relab_spectrum_id}.tab"
        dest = RAW / "spectra" / row.subdir / f"{row.relab_spectrum_id}.tab"
        fetch(session, f"{BASE}/{rel}", dest)
        wl, rf = parse_tab(dest)
        factor = {"nm": 1e-3}.get(row.range_unit, 1.0)
        wl_um = [w * factor for w in wl]
        spectra_rows.extend((row.relab_spectrum_id, w, r) for w, r in zip(wl_um, rf))
        npts.append(len(wl_um))
        wlmin.append(min(wl_um) if wl_um else np.nan)
        wlmax.append(max(wl_um) if wl_um else np.nan)
    ras["n_points"], ras["wl_min_um"], ras["wl_max_um"] = npts, wlmin, wlmax

    spectra = pd.DataFrame(spectra_rows,
                           columns=["relab_spectrum_id", "wavelength_um", "reflectance"])
    # Same QA as A7: whole spectra with reflectance >= 1.5 are specular/
    # mis-scaled, never diffuse templates - exclude the spectrum, never edit.
    bad = set(spectra.loc[spectra["reflectance"] >= 1.5, "relab_spectrum_id"].unique())
    if bad:
        print(f"[a7b] QA-excluded {len(bad)} spectra (reflectance >= 1.5): {sorted(bad)}")
        ras = ras[~ras["relab_spectrum_id"].isin(bad)]
        spectra = spectra[~spectra["relab_spectrum_id"].isin(bad)]
    ras = ras[ras["relab_spectrum_id"].isin(set(spectra["relab_spectrum_id"]))]

    ras["measurement_type"] = (ras["measurement_type_raw"] + "/"
                               + ras["geometry"].replace("", "unspecified"))

    samples = ras[[
        "relab_spectrum_id", "relab_sample_id", "sample_name", "parent_body",
        "mission", "grain_size_min_um", "grain_size_max_um", "relab_texture",
        "wl_min_um", "wl_max_um", "n_points", "measurement_type",
        "collection_location", "owner_name", "measurement_date",
    ]].reset_index(drop=True)

    INTERIM.mkdir(parents=True, exist_ok=True)
    s_out = INTERIM / "relab_returned_samples.parquet"
    sp_out = INTERIM / "relab_returned_spectra.parquet"
    samples.to_parquet(s_out, index=False)
    spectra.reset_index(drop=True).to_parquet(sp_out, index=False)

    relab_prov = json.loads((PROVENANCE / "relab.json").read_text())
    grain_doc = (samples["grain_size_min_um"].notna()
                 | samples["grain_size_max_um"].notna()).mean()
    prov = {
        "source": relab_prov["source"],
        "url": relab_prov["url"],
        "urn": relab_prov["urn"],
        "version": relab_prov["version"],
        "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "derived_from": "data/provenance/relab.json (same PDS4 bundle, no re-crawl)",
        "selection": ("specimen_type == 'Returned Asteroid Sample'. "
                      "Returned Lunar Samples (1,510) deliberately excluded: the "
                      "Moon is not an asteroid (out of scope for an asteroid "
                      "space-weathering reference)."),
        "bodies_present": sorted(samples["parent_body"].unique().tolist()),
        "bodies_absent_in_bundle": ["Itokawa", "Bennu"],
        "row_count": int(len(samples)),
        "spectra_points": int(len(spectra)),
        "grain_size_documented_fraction": round(float(grain_doc), 3),
        "qa_reflectance_excluded": sorted(bad),
        "sample_note": ("All returned-asteroid spectra in this bundle are Ryugu "
                        "sample C0002 (Hayabusa2 / JAXA), measured at the Brown "
                        "University RELAB facility (owner/requestor Ralph E. "
                        "Milliken). Collection location per label e.g. 'Asteroid "
                        "Ryugu SCI crater' (artificial-impact subsurface material). "
                        "The PDS4 labels carry no per-sample paper DOI; sample "
                        "origin traces to the Hayabusa2 return (Yada et al. 2022, "
                        "Nat. Astron. 6, 214; Pilorget et al. 2022)."),
        "citation": relab_prov["citation"],
        "license_terms": relab_prov["license_terms"],
    }
    (PROVENANCE / "relab_returned.json").write_text(json.dumps(prov, indent=2))

    print("\n=== Coverage report (AP_A7b) ===")
    print(f"Returned-sample spectra kept:        {len(samples)}")
    print("Per body / per mission:")
    tab = samples.groupby(["parent_body", "mission"]).size()
    for (body, mission), n in tab.items():
        print(f"    {body:<10} {mission:<12} {n}")
    print(f"Wavelength coverage span:            {samples['wl_min_um'].min():.3f}"
          f"-{samples['wl_max_um'].max():.3f} um")
    visnir = samples[(samples["wl_min_um"] <= 0.8) & (samples["wl_max_um"] >= 2.4)]
    print(f"  spectra spanning VIS-NIR (<=0.8 & >=2.4 um): {len(visnir)} "
          f"(the C2 space-weathering window)")
    print(f"Grain-size documented fraction:      {grain_doc*100:.1f}%")
    print(f"Saved: {s_out.name} ({s_out.stat().st_size/1024:.1f} KB), "
          f"{sp_out.name} ({sp_out.stat().st_size/1024:.1f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
