#!/usr/bin/env python3
"""AP_A7: RELAB meteorite spectral library.

Source: PDS4 bundle urn:nasa:pds:relab at the PDS Geosciences Node
(https://pds-geosciences.wustl.edu/speclib/urn-nasa-pds-relab/).
The legacy Brown RelabDB zip (Sample_Catalogue/Spectra_Catalogue) is no
longer distributed (sites.brown.edu links now point to the PDS
spectral library), so the PDS4 bundle is both the primary and the only
bulk-accessible source. It has no monolithic catalog table: per-spectrum
metadata lives in the PDS4 label of each product. This script therefore
crawls the labels once (file list taken from the bundle's MD5 manifest,
everything cached under data/raw/relab/), extracts metadata, applies the
filter funnel and downloads the data tables only for kept spectra.

Metadata sources inside each label (speclib dictionary):
  specimen_type          -> meteorite test (criterion a)
  measurement_type       -> reflectance test (criterion b)
  spectral_range_min/max -> Gaia range test (criterion c), unit-aware
  material_subtype       -> relab_texture (slab | chip | particulate)
  specimen_name/descr.   -> grain size (parsed from text like
                            "0-125 um", "<45 um"; the PDS4 labels have
                            no structured grain-size field) and
                            meteorite class token (e.g. "CM2", "L5",
                            "Howardite" - no structured class field
                            either); both rules documented here.

Group mapping: extracted raw class tokens are mapped to canonical
groups via data/interim/relab_group_map.csv. If that file already
exists it WINS (it is a reviewed artifact - edits by the reviewer
persist across re-runs); otherwise it is generated from the built-in
rule table below and written for review. Unmappable raw values get
relab_group = "unmapped", never dropped.

Archive values stay pure: no resampling, smoothing or normalization.
The only transformation is wavelength unit conversion to micrometers.
"""

import hashlib
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "relab"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"
DOCS = ROOT / "docs"

BASE = "https://pds-geosciences.wustl.edu/speclib/urn-nasa-pds-relab"
PDS_URN = "urn:nasa:pds:relab"
UA = {"User-Agent": "esox-master-catalog/A7 (vienna@esoxconsult.com)"}

TODAY = datetime.now(timezone.utc).strftime("%Y%m%d")

SPECLIB_NS = "{http://pds.nasa.gov/pds4/speclib/v1}"

GROUP_RULES = [
    (r"^CI\d?(\.\d)?$", "CI"), (r"^CM\d?(\.\d)?$", "CM"),
    (r"^CO\d?(\.\d)?$", "CO"), (r"^CV\d?(\.\d)?$", "CV"),
    (r"^CK\d?(\.\d)?$", "CK"), (r"^CR\d?(\.\d)?$", "CR"),
    (r"^CH\d?$", "CH"), (r"^CB\d?[ab]?$", "CB"),
    (r"^C\d?[ -]?(UNG|UNGROUPED)$", "C-ungrouped"),
    (r"^TAGISH LAKE$", "C-ungrouped"),
    (r"^LL\d(\.\d)?$", "LL"), (r"^LL$", "LL"),
    (r"^L\d(\.\d)?$", "L"), (r"^L$", "L"),
    (r"^H\d(\.\d)?$", "H"), (r"^H$", "H"),
    (r"^L/LL\d?$", "LL"), (r"^H/L\d?$", "L"),
    (r"^EH\d?$", "EH"), (r"^EL\d?$", "EL"), (r"^E\d?$", "E"),
    (r"^R\d(\.\d)?$", "R"), (r"^R$", "R"),
    (r"^HOWARDITE$", "howardite"), (r"^EUCRITE$", "eucrite"),
    (r"^DIOGENITE$", "diogenite"), (r"^AUBRITE$", "aubrite"),
    (r"^UREILITE$", "ureilite"), (r"^ANGRITE$", "angrite"),
    (r"^BRACHINITE$", "brachinite"),
    (r"^(LODRANITE|ACAPULCOITE)$", "lodranite-acapulcoite"),
    (r"^WINONAITE$", "winonaite"),
    (r"^MESOSIDERITE$", "mesosiderite"), (r"^PALLASITE$", "pallasite"),
    (r"^(IRON|I{1,3}[A-G]B?|IV[AB]|OCTAHEDRITE|HEXAHEDRITE|ATAXITE)$", "iron"),
    (r"^LUNAR( METEORITE)?$", "lunar-meteorite"),
    (r"^(MARTIAN|SNC|SHERGOTTITE|NAKHLITE|CHASSIGNITE)$", "martian"),
]

CLASS_TOKEN_RE = re.compile(
    r"\b(CI\d?(?:\.\d)?|CM\d?(?:\.\d)?|CO\d(?:\.\d)?|CV\d(?:\.\d)?"
    r"|CK\d(?:\.\d)?|CR\d(?:\.\d)?|CH\d?|CB\d?[ab]?"
    r"|C\d?[ -]?ung(?:rouped)?"
    r"|LL\d(?:\.\d)?|L/LL\d?|H/L\d?"
    r"|EH\d?|EL\d?"
    r"|howardite|eucrite|diogenite|aubrite|ureilite|angrite|brachinite"
    r"|lodranite|acapulcoite|winonaite|mesosiderite|pallasite"
    r"|shergottite|nakhlite|chassignite|SNC"
    r"|octahedrite|hexahedrite|ataxite"
    r"|I{1,3}[A-G]B?(?![a-z])|IV[AB](?![a-z])"
    r"|tagish lake"
    r"|[HLER]\d(?:\.\d)?)\b",
    re.IGNORECASE)

ORDINARY_HINT_RE = re.compile(
    r"(?:ordinary\s+)?chondrite\s*\(?\s*(LL|L|H|EH|EL|E|R)\s*\d?\s*\)?",
    re.IGNORECASE)

GRAIN_RANGE_RE = re.compile(
    r"(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)\s*(?:um|µm|μm|micron)", re.I)
GRAIN_LT_RE = re.compile(r"<\s*(\d+(?:\.\d+)?)\s*(?:um|µm|μm|micron)", re.I)
GRAIN_GT_RE = re.compile(r">\s*(\d+(?:\.\d+)?)\s*(?:um|µm|μm|micron)", re.I)


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def fetch(session, url, dest, retries=3):
    if dest.exists() and dest.stat().st_size > 0:
        return dest, False
    for attempt in range(1, retries + 1):
        try:
            r = session.get(url, headers=UA, timeout=(15, 120))
            r.raise_for_status()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(r.content)
            return dest, True
        except requests.RequestException as e:
            if attempt == retries:
                raise
            time.sleep(3 * attempt)


def label_list(session):
    """All reflectance label paths (subdir/name.xml) from the MD5 manifest."""
    manifest, _ = fetch(session, f"{BASE}/urn-nasa-pds-relab.md5",
                        RAW / f"urn-nasa-pds-relab_{TODAY}.md5")
    paths = []
    for line in manifest.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        rel = parts[1].strip().replace("\\", "/").lstrip("/")
        if re.match(r"^data_reflectance/(bdr|ftir)\w*/[^/]+\.xml$", rel):
            paths.append(rel)
    return sorted(paths)


def crawl_labels(session, paths):
    labels_dir = RAW / "labels"
    n_new = 0
    for i, rel in enumerate(paths, 1):
        dest = labels_dir / Path(rel).parent.name / Path(rel).name
        _, downloaded = fetch(session, f"{BASE}/{rel}", dest)
        if downloaded:
            n_new += 1
            time.sleep(0.005)
        if i % 2000 == 0:
            print(f"[relab]   labels {i:,}/{len(paths):,} ({n_new:,} new)")
    print(f"[relab] label crawl done: {len(paths):,} labels ({n_new:,} downloaded)")
    return labels_dir


def txt(el):
    return (el.text or "").strip() if el is not None else ""


def parse_label(path):
    root = ET.parse(path).getroot()
    tags = {}
    for el in root.iter():
        tag = el.tag.split("}")[-1]
        tags.setdefault(tag, []).append((el.text or "").strip())

    def one(name, default=""):
        v = tags.get(name, [])
        return v[0] if v else default

    def all_(name):
        return [v for v in tags.get(name, []) if v]

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan

    unit = (one("spectral_range_unit_name") or "nm").lower()
    factor = {"nm": 1e-3, "nanometer": 1e-3, "um": 1.0, "micrometer": 1.0,
              "micron": 1.0}.get(unit)
    rmins = [num(v) for v in all_("spectral_range_min")]
    rmaxs = [num(v) for v in all_("spectral_range_max")]
    return {
        "relab_spectrum_id": path.stem,
        "subdir": path.parent.name,
        "relab_sample_id": one("specimen_id"),
        "relab_sample_name": one("specimen_name"),
        "description": one("specimen_description"),
        "specimen_type": one("specimen_type"),
        "material_subtypes": "|".join(all_("material_subtype")),
        "measurement_type_raw": one("measurement_type"),
        "geometry": one("measurement_geometry_type"),
        "range_unit": unit,
        "range_min_um": min(rmins) * factor if rmins and factor else np.nan,
        "range_max_um": max(rmaxs) * factor if rmaxs and factor else np.nan,
        "incidence_deg": num(one("incidence_angle")),
        "emission_deg": num(one("emission_angle")),
    }


def texture_of(subtypes):
    s = subtypes.lower()
    if "slab" in s:
        return "slab"
    if "chip" in s:
        return "chip"
    if "particulate" in s:
        return "particulate"
    return "unknown"


def grain_size_of(text):
    m = GRAIN_RANGE_RE.search(text)
    if m:
        # a few archive names write the range reversed ("150-125 um");
        # normalize so min <= max
        a, b = float(m.group(1)), float(m.group(2))
        return min(a, b), max(a, b)
    m = GRAIN_LT_RE.search(text)
    if m:
        return np.nan, float(m.group(1))
    m = GRAIN_GT_RE.search(text)
    if m:
        return float(m.group(1)), np.nan
    return np.nan, np.nan


def raw_class_of(row):
    """Extract a verbatim class token from type/name/description."""
    stype = row["specimen_type"].lower()
    if "lunar" in stype:
        return "Lunar Meteorite"
    if "martian" in stype:
        return "Martian"
    text = f"{row['relab_sample_name']} {row['description']}"
    m = ORDINARY_HINT_RE.search(text)
    if m:
        hint = m.group(1).upper()
        m2 = re.search(rf"\b{hint}\s?(\d(?:\.\d)?)\b", text, re.I)
        return f"{hint}{m2.group(1)}" if m2 else hint
    m = CLASS_TOKEN_RE.search(text)
    if m:
        return m.group(1)
    return None


def canonical_group(raw):
    if raw is None:
        return "unmapped"
    token = raw.strip().upper()
    for pattern, canon in GROUP_RULES:
        if re.match(pattern, token, re.IGNORECASE):
            return canon
    return "unmapped"


def load_or_build_group_map(raws):
    """The mapping CSV is a reviewed artifact: if present it wins."""
    map_path = INTERIM / "relab_group_map.csv"
    if map_path.exists():
        gm = pd.read_csv(map_path)
        print(f"[relab] using existing (reviewed) {map_path.name}: {len(gm)} entries")
        return dict(zip(gm["raw"], gm["canonical"])), map_path
    rows = sorted({(r, canonical_group(r)) for r in raws if r is not None})
    gm = pd.DataFrame(rows, columns=["raw", "canonical"])
    gm.to_csv(map_path, index=False)
    print(f"[relab] wrote {map_path.name}: {len(gm)} raw values "
          "(reviewed artifact - edits persist)")
    return dict(zip(gm["raw"], gm["canonical"])), map_path


def parse_tab(path):
    wl, rf = [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        toks = line.split()
        if len(toks) < 2:
            continue
        try:
            w, r = float(toks[0]), float(toks[1])
        except ValueError:
            continue
        wl.append(w)
        rf.append(r)
    return wl, rf


def main():
    for d in (RAW, INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)
    session = requests.Session()

    bundle_xml, _ = fetch(session, f"{BASE}/bundle_relab.xml",
                          RAW / f"bundle_relab_{TODAY}.xml")
    broot = ET.parse(bundle_xml).getroot()
    version = None
    for el in broot.iter():
        if el.tag.split("}")[-1] == "version_id":
            version = el.text.strip()
            break
    print(f"[relab] bundle {PDS_URN} version_id: {version}")

    paths = label_list(session)
    print(f"[relab] manifest lists {len(paths):,} reflectance labels")
    labels_dir = crawl_labels(session, paths)

    print("[relab] parsing labels ...")
    meta = pd.DataFrame([parse_label(labels_dir / Path(p).parent.name / Path(p).name)
                         for p in paths])
    total = len(meta)

    type_counts = meta["specimen_type"].value_counts()
    print("[relab] specimen_type counts (top 12):")
    for t, c in type_counts.head(12).items():
        print(f"    {t:<28} {c:,}")

    is_met = meta["specimen_type"].str.lower().str.contains("meteorite")
    met = meta[is_met].copy()

    is_refl = met["measurement_type_raw"].str.lower().eq("reflectance")
    refl = met[is_refl].copy()

    in_range = (refl["range_min_um"] <= 0.40) & (refl["range_max_um"] >= 1.00)
    gaia = refl[in_range].copy()

    gaia["relab_texture"] = gaia["material_subtypes"].map(texture_of)
    sizes = gaia.apply(
        lambda r: grain_size_of(f"{r['relab_sample_name']} {r['description']}"),
        axis=1)
    gaia["grain_size_min_um"] = [s[0] for s in sizes]
    gaia["grain_size_max_um"] = [s[1] for s in sizes]
    documented = (gaia["grain_size_min_um"].notna()
                  | gaia["grain_size_max_um"].notna()
                  | gaia["relab_texture"].isin(["slab", "chip"]))
    kept = gaia[documented].copy()

    print(f"[relab] funnel: total {total:,} -> meteorite {len(met):,} "
          f"-> reflectance {len(refl):,} -> Gaia-range {len(gaia):,} "
          f"-> grain/texture documented {len(kept):,}")

    kept["relab_group_raw"] = kept.apply(raw_class_of, axis=1)
    gmap, map_path = load_or_build_group_map(kept["relab_group_raw"].dropna().unique())
    kept["relab_group"] = kept["relab_group_raw"].map(
        lambda r: gmap.get(r, canonical_group(r)) if r is not None else "unmapped")
    kept["relab_group"] = kept["relab_group"].fillna("unmapped")

    kept["meteorite_name"] = (
        kept["relab_sample_name"]
        .str.replace(GRAIN_RANGE_RE, "", regex=True)
        .str.replace(GRAIN_LT_RE, "", regex=True)
        .str.replace(GRAIN_GT_RE, "", regex=True)
        .str.strip(" ,;-"))
    kept["measurement_type"] = (kept["measurement_type_raw"] + "/"
                                + kept["geometry"].replace("", "unspecified"))

    print(f"[relab] downloading {len(kept):,} kept data tables ...")
    spectra_rows = []
    n_pts_col, wlmin_col, wlmax_col = [], [], []
    for i, row in enumerate(kept.itertuples(), 1):
        rel = f"data_reflectance/{row.subdir}/{row.relab_spectrum_id}.tab"
        dest = RAW / "spectra" / row.subdir / f"{row.relab_spectrum_id}.tab"
        fetch(session, f"{BASE}/{rel}", dest)
        wl, rf = parse_tab(dest)
        factor = {"nm": 1e-3}.get(row.range_unit, 1.0)
        wl_um = [w * factor for w in wl]
        spectra_rows.extend(
            (row.relab_spectrum_id, w, r) for w, r in zip(wl_um, rf))
        n_pts_col.append(len(wl_um))
        wlmin_col.append(min(wl_um) if wl_um else np.nan)
        wlmax_col.append(max(wl_um) if wl_um else np.nan)
        if i % 500 == 0:
            print(f"[relab]   spectra {i:,}/{len(kept):,}")
    kept["n_points"] = n_pts_col
    kept["wl_min_um"] = wlmin_col
    kept["wl_max_um"] = wlmax_col

    spectra_all = pd.DataFrame(spectra_rows, columns=[
        "relab_spectrum_id", "wavelength_um", "reflectance"])
    # QA: spectra containing reflectance >= 1.5 are not diffuse
    # reflectance templates (specular slab/mirror geometry, e.g. the
    # Mundrabilla "mirror surface" run) or are mis-scaled in the
    # archive; exclude the whole spectrum, never edit values.
    bad_ids = set(spectra_all.loc[spectra_all["reflectance"] >= 1.5,
                                  "relab_spectrum_id"].unique())
    if bad_ids:
        print(f"[relab] QA-excluded {len(bad_ids)} spectra with reflectance "
              f">= 1.5 (specular/mis-scaled): {sorted(bad_ids)}")
        kept = kept[~kept["relab_spectrum_id"].isin(bad_ids)]
        spectra_rows = [r for r in spectra_rows if r[0] not in bad_ids]

    samples = kept[[
        "relab_spectrum_id", "relab_sample_id", "relab_sample_name",
        "meteorite_name", "relab_group_raw", "relab_group", "relab_texture",
        "grain_size_min_um", "grain_size_max_um", "wl_min_um", "wl_max_um",
        "n_points", "measurement_type", "incidence_deg", "emission_deg",
    ]].reset_index(drop=True)
    spectra = pd.DataFrame(spectra_rows, columns=[
        "relab_spectrum_id", "wavelength_um", "reflectance"])

    samples_out = INTERIM / "relab_samples.parquet"
    spectra_out = INTERIM / "relab_spectra.parquet"
    samples.to_parquet(samples_out, engine="pyarrow", compression="snappy", index=False)
    spectra.to_parquet(spectra_out, engine="pyarrow", compression="snappy", index=False)

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    group_counts = samples["relab_group"].value_counts()
    viable = group_counts[(group_counts >= 5) & (group_counts.index != "unmapped")]
    n_ll = int((samples["relab_group"] == "LL").sum())

    PROVENANCE.joinpath("relab.json").write_text(json.dumps({
        "source": "RELAB Spectral Library, PDS4 bundle at the PDS Geosciences Node",
        "source_note": ("legacy Brown RelabDB zip no longer distributed; PDS4 "
                        "bundle is the primary and only bulk source. Metadata "
                        "extracted from per-product PDS4 labels (no monolithic "
                        "catalog exists in the bundle); grain size and meteorite "
                        "class parsed from specimen name/description text with "
                        "documented rules (no structured fields available)."),
        "url": BASE + "/",
        "urn": PDS_URN,
        "version": version,
        "date": retrieved,
        "manifest_sha256": sha256_file(RAW / f"urn-nasa-pds-relab_{TODAY}.md5"),
        "labels_crawled": int(total),
        "row_count": int(len(samples)),
        "spectra_points": int(len(spectra)),
        "funnel": {
            "all_reflectance_products": int(total),
            "meteorite": int(len(met)),
            "reflectance_measurement": int(len(refl)),
            "gaia_range_040_100um": int(len(gaia)),
            "grain_or_texture_documented": int(len(kept) + len(bad_ids)),
            "qa_reflectance_excluded": int(len(bad_ids)),
            "kept": int(len(kept)),
        },
        "citation": ("Pieters, C.M., Hiroi, T. (2004), RELAB (Reflectance "
                     "Experiment Laboratory): A NASA Multiuser Spectroscopy "
                     "Facility, LPSC XXXV, #1720; Milliken, R. et al., RELAB "
                     f"Spectral Library Bundle {PDS_URN} V{version}, PDS "
                     "Geosciences Node; PDS Spectral Library: "
                     f"doi:10.3847/PSJ/ad5af7; retrieved {retrieved}"),
        "license_terms": ("NASA PDS open data (publicly released, no access "
                          "restrictions). RELAB is a NASA-supported multiuser "
                          "facility; publications using these data customarily "
                          "acknowledge the RELAB facility at Brown University "
                          "and cite Pieters & Hiroi (2004) and the PDS "
                          "Spectral Library."),
    }, indent=2))
    print("[provenance] wrote data/provenance/relab.json")

    print()
    print("=== Coverage report (AP_A7) ===")
    print(f"RELAB spectra total (all holdings):     {total:,}")
    print(f"Funnel: meteorite {len(met):,} -> reflectance {len(refl):,}")
    print(f"        -> Gaia-range {len(gaia):,} -> grain/texture documented "
          f"{len(kept):,}  = KEPT")
    print(f"Kept spectra / distinct meteorites:     {len(samples):,} / "
          f"{samples['meteorite_name'].str.lower().nunique():,}")
    n_groups = int((group_counts.index != "unmapped").sum())
    print(f"Groups covered / unmapped raw values:   {n_groups} / "
          f"{int((samples['relab_group'] == 'unmapped').sum()):,} spectra unmapped")
    print("Per-group spectrum counts:")
    for g, c in group_counts.items():
        print(f"    {g:<24} {c:,}")
    print(f"Groups with >= 5 spectra:               {len(viable)}")
    print(f"LL chondrite spectra:                   {n_ll}")
    print(f"Saved: data/interim/relab_samples.parquet  "
          f"{samples_out.stat().st_size/1e6:.2f} MB")
    print(f"Saved: data/interim/relab_spectra.parquet  "
          f"{spectra_out.stat().st_size/1e6:.2f} MB")
    print(f"Saved: data/interim/relab_group_map.csv")

    status = DOCS / "apophis_status.md"
    if status.exists() and "| A7 |" not in status.read_text():
        lines = status.read_text().splitlines()
        last_row = max(i for i, l in enumerate(lines) if l.startswith("|"))
        lines.insert(last_row + 1,
                     f"| A7 | RELAB (lab reference set) | **{n_ll} LL-chondrite spectra** "
                     f"available as the laboratory counterpart to Apophis' literature "
                     f"Sq classification (Sq ↔ LL connection, Binzel et al. 2009) |")
        status.write_text("\n".join(lines) + "\n")
        print("[docs] appended LL-chondrite line to docs/apophis_status.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
