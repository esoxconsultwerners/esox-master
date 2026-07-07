#!/usr/bin/env python3
"""AP_A5: SMASS/MITHNEOS + ECAS ground-truth spectra.

Three PDS4 bundles from the NASA PDS Small Bodies Node, each downloaded
once as the versioned bundle zip and read directly from the cached zip.
Every table layout is parsed from its PDS4 XML label, never guessed.

  smass2:   gbo.ast.smass2.spectra::1.0  (Bus & Binzel 2002, VIS)
  mithneos: gbo.mithneos.spectra_2000-2021::1.0  (Binzel et al. 2019 +
            Marsset et al. 2022 collections, NIR)
  ecas:     gbo.ast.ecas.phot::1.0  (Zellner, Tholen, Tedesco 1985,
            eight-band photometry -> seven mean color indices vs v)

No resampling, no renormalization - archive values as published.
Taxonomy: none of the three bundles distributes per-object classes
(Tholen/Bus-DeMeo classes live in the separate PDS Asteroid Taxonomy
dataset), so gt_taxon is null throughout this package.

The long table carries gt_spectrum_id (source file stem) so that
multiple observation epochs of the same object stay separable.
"""

import csv
import hashlib
import io
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"

PDS4_BASE = "https://sbnarchive.psi.edu/pds4/non_mission"

SOURCES = {
    "smass2": {
        "zip_url": f"{PDS4_BASE}/gbo.ast.smass2.spectra.zip",
        "urn": "urn:nasa:pds:gbo.ast.smass2.spectra::1.0",
        "doi": "10.26033/fj1d-vb37",
        "citation": ("Bus, S.J., Binzel, R.P. (2002), Phase II of the Small "
                     "Main-Belt Asteroid Spectroscopic Survey: The Observations, "
                     "Icarus 158, 106-145; PDS4 bundle "
                     "urn:nasa:pds:gbo.ast.smass2.spectra::1.0, "
                     "NASA Planetary Data System (2020), doi:10.26033/fj1d-vb37"),
    },
    "mithneos": {
        "zip_url": f"{PDS4_BASE}/gbo.ast.mithneos.spectra_2000-2021_V1_0.zip",
        "urn": "urn:nasa:pds:gbo.mithneos.spectra_2000-2021::1.0",
        "doi": "10.26033/1aft-4018",
        "citation": ("Binzel, R.P., et al. (2019), Compositional distributions and "
                     "evolutionary processes for the near-Earth object population, "
                     "Icarus 324, 41-76; Marsset, M., et al. (2022), AJ 163, 165; "
                     "MITHNEOS IRTF Spectra of Asteroids from 2000 to 2021 Bundle V1.0, "
                     "urn:nasa:pds:gbo.mithneos.spectra_2000-2021::1.0, "
                     "NASA Planetary Data System (2023), doi:10.26033/1aft-4018"),
    },
    "ecas": {
        "zip_url": f"{PDS4_BASE}/gbo.ast.ecas.phot.zip",
        "urn": "urn:nasa:pds:gbo.ast.ecas.phot::1.0",
        "doi": "10.26033/dhpe-7k98",
        "citation": ("Zellner, B., Tholen, D.J., Tedesco, E.F. (1985), The "
                     "Eight-Color Asteroid Survey: Results for 589 minor planets, "
                     "Icarus 61, 355-416; PDS4 bundle "
                     "urn:nasa:pds:gbo.ast.ecas.phot::1.0, "
                     "NASA Planetary Data System (2020), doi:10.26033/dhpe-7k98"),
    },
}

TODAY = datetime.now(timezone.utc).strftime("%Y%m%d")


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[gt] using cached {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    print(f"[gt] downloading {url}")
    tmp = dest.with_suffix(".part")
    with requests.get(url, headers={"User-Agent": "esox-master-catalog/A5"},
                      stream=True, timeout=(30, 1800)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[gt] saved {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


def strip_ns(tag):
    return tag.split("}")[-1]


def parse_char_label(xml_bytes):
    root = ET.fromstring(xml_bytes)
    fields, records = [], None
    for el in root.iter():
        if strip_ns(el.tag) == "Table_Character":
            for sub in el.iter():
                st = strip_ns(sub.tag)
                if st == "records":
                    records = int(sub.text)
                elif st == "Field_Character":
                    name = loc = length = None
                    for f in sub:
                        ft = strip_ns(f.tag)
                        if ft == "name":
                            name = f.text.strip()
                        elif ft == "field_location":
                            loc = int(f.text)
                        elif ft == "field_length":
                            length = int(f.text)
                    fields.append((name, loc, length))
    if not fields:
        raise ValueError("no Table_Character spec in label")
    return fields, records


def parse_delim_label(xml_bytes):
    root = ET.fromstring(xml_bytes)
    names, records = [], None
    for el in root.iter():
        if strip_ns(el.tag) == "Table_Delimited":
            for sub in el.iter():
                st = strip_ns(sub.tag)
                if st == "records":
                    records = int(sub.text)
                elif st == "Field_Delimited":
                    for f in sub:
                        if strip_ns(f.tag) == "name":
                            names.append(f.text.strip())
    if not names:
        raise ValueError("no Table_Delimited spec in label")
    return names, records


def to_num(s):
    s = (s or "").strip()
    if not s or s == "-":
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def normalize_desig(token):
    """'1995bm2' -> '1995 BM2'; anything else uppercased as-is."""
    m = re.match(r"^(\d{4})([a-z]{1,2}\d*)$", token)
    if m:
        return f"{m.group(1)} {m.group(2).upper()}"
    return token.upper()


def object_key_from_stem(stem):
    """Spectrum file stem '<objid>_<suffix>' -> (number_mp, designation)."""
    objid = stem.rsplit("_", 1)[0]
    if objid.isdigit():
        return int(objid), None
    return None, normalize_desig(objid)


def bundle_version(zf, names):
    for n in names:
        if re.search(r"/bundle_[^/]*\.xml$", n) or re.match(r"bundle_[^/]*\.xml$", n):
            root = ET.fromstring(zf.read(n))
            for el in root.iter():
                if strip_ns(el.tag) == "version_id":
                    return el.text.strip()
    return None


def ingest_smass2(zf, names):
    long_rows = []
    labels = sorted(n for n in names
                    if re.search(r"/data/data\d+/[^/]+\.xml$", n))
    for label_name in labels:
        tab_name = label_name[:-4] + ".tab"
        if tab_name not in names:
            continue
        fields, _ = parse_char_label(zf.read(label_name))
        idx = {n: (loc, ln) for n, loc, ln in fields}
        for need in ("WAVELENGTH", "REFLECTANCE", "UNCERTAINTY"):
            if need not in idx:
                raise ValueError(f"{label_name}: missing field {need}")
        stem = Path(tab_name).stem
        num, desig = object_key_from_stem(stem)
        for line in zf.read(tab_name).decode("utf-8", "replace").splitlines():
            if not line.strip():
                continue
            vals = {n: to_num(line[loc - 1:loc - 1 + ln]) for n, (loc, ln) in idx.items()}
            long_rows.append((num, desig, "smass2", stem,
                              vals["WAVELENGTH"], vals["REFLECTANCE"],
                              vals["UNCERTAINTY"]))
    n_spec = len({r[3] for r in long_rows})
    print(f"[smass2] {n_spec} spectra, {len(long_rows):,} points")
    return long_rows


def mithneos_obsparams(zf, names):
    """basename -> (number, designation) from the two obs parameter CSVs."""
    key = {}
    for coll in ("binzel", "marsset"):
        label = next(n for n in names
                     if n.endswith(f"observationalparameters_{coll}.xml"))
        csv_name = label[:-4] + ".csv"
        fields, _ = parse_delim_label(zf.read(label))
        idx = {f: i for i, f in enumerate(fields)}
        for r in csv.reader(io.StringIO(zf.read(csv_name).decode("utf-8", "replace"))):
            if not r or not any(c.strip() for c in r):
                continue
            num = r[idx["ASTEROID_NUMBER"]].strip()
            desig = r[idx["ASTEROID_DESIGNATION"]].strip()
            base = Path(r[idx["FILENAME"]].strip()).name
            key[base] = (int(num) if num.isdigit() and int(num) > 0 else None,
                         desig or None)
    return key


def ingest_mithneos(zf, names):
    obsmap = mithneos_obsparams(zf, names)
    long_rows = []
    labels = sorted(n for n in names
                    if re.search(r"/data/(binzel2019|marsset2022)/[^/]+\.xml$", n))
    n_unmapped = 0
    for label_name in labels:
        csv_name = label_name[:-4] + ".csv"
        if csv_name not in names:
            continue
        fields, _ = parse_delim_label(zf.read(label_name))
        idx = {f.lower(): i for i, f in enumerate(fields)}
        for need in ("wavelength", "reflectance", "error"):
            if need not in idx:
                raise ValueError(f"{label_name}: missing field {need}")
        stem = Path(csv_name).stem
        base = Path(csv_name).name
        if base in obsmap:
            num, desig = obsmap[base]
        else:
            n_unmapped += 1
            num, desig = object_key_from_stem(stem)
        for r in csv.reader(io.StringIO(zf.read(csv_name).decode("utf-8", "replace"))):
            if not r or not any(c.strip() for c in r):
                continue
            long_rows.append((num, desig, "mithneos", stem,
                              to_num(r[idx["wavelength"]]),
                              to_num(r[idx["reflectance"]]),
                              to_num(r[idx["error"]])))
    if n_unmapped:
        print(f"[mithneos] {n_unmapped} spectra not in obs param tables "
              "(object key parsed from file name)")
    n_spec = len({r[3] for r in long_rows})
    print(f"[mithneos] {n_spec} spectra, {len(long_rows):,} points")
    return long_rows


ECAS_COLORS = ["S_V", "U_V", "B_V", "V_W", "V_X", "V_P", "V_Z"]


def ingest_ecas(zf, names):
    label_name = next(n for n in names if n.endswith("data/ecasmean.xml"))
    tab_name = label_name[:-4] + ".tab"
    fields, records = parse_char_label(zf.read(label_name))
    idx = {n: (loc, ln) for n, loc, ln in fields}
    need = {"AST_NUMBER", "AST_NAME", "NIGHTS"} | \
        {f"{c}_MEAN" for c in ECAS_COLORS} | {f"{c}_STD_DEV" for c in ECAS_COLORS}
    missing = need - set(idx)
    if missing:
        raise ValueError(f"ecasmean label lacks fields: {missing}")

    rows = []
    for line in zf.read(tab_name).decode("utf-8", "replace").splitlines():
        if not line.strip():
            continue
        def cut(name):
            loc, ln = idx[name]
            return line[loc - 1:loc - 1 + ln].strip()
        num = cut("AST_NUMBER")
        row = {
            "number_mp": int(num) if num.isdigit() and int(num) > 0 else None,
            "designation": None if num.isdigit() and int(num) > 0 else cut("AST_NAME"),
            "ecas_nights": to_num(cut("NIGHTS")),
        }
        for c in ECAS_COLORS:
            row[f"ecas_{c.lower()}"] = to_num(cut(f"{c}_MEAN"))
            row[f"ecas_{c.lower()}_err"] = to_num(cut(f"{c}_STD_DEV"))
        rows.append(row)
    if records is not None and len(rows) != records:
        print(f"[ecas] WARNING parsed {len(rows)} rows, label declares {records}")
    print(f"[ecas] {len(rows)} objects (mean colors table)")
    return pd.DataFrame(rows)


def spectra_summary(long_df, source):
    """One row per object: canonical key is number_mp where known, else
    designation (the same object can arrive with and without a designation
    depending on whether the obs parameter tables cover its file)."""
    sub = long_df[long_df["gt_source"] == source].copy()
    sub["obj_key"] = sub["number_mp"].astype("string").fillna(sub["designation"])
    g = sub.groupby("obj_key", dropna=False)
    out = g.agg(
        number_mp=("number_mp", "first"),
        designation=("designation",
                     lambda s: s.dropna().iloc[0] if s.notna().any() else None),
        gt_wl_min_um=("wavelength_um", "min"),
        gt_wl_max_um=("wavelength_um", "max"),
        gt_n_points=("wavelength_um", "size"),
        gt_n_spectra=("gt_spectrum_id", "nunique"),
    ).reset_index(drop=True)
    out["number_mp"] = out["number_mp"].astype("Int64")
    out["gt_source"] = source
    out["gt_taxon"] = None
    return out


def main():
    for d in (RAW, INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)

    zips, versions = {}, {}
    for src, cfg in SOURCES.items():
        raw_dir = RAW / src
        raw_dir.mkdir(parents=True, exist_ok=True)
        zips[src] = download(cfg["zip_url"],
                             raw_dir / f"{Path(cfg['zip_url']).stem}_{TODAY}.zip")

    smass_zf = zipfile.ZipFile(zips["smass2"])
    mith_zf = zipfile.ZipFile(zips["mithneos"])
    ecas_zf = zipfile.ZipFile(zips["ecas"])
    smass_names = set(smass_zf.namelist())
    mith_names = set(mith_zf.namelist())
    ecas_names = set(ecas_zf.namelist())
    for src, zf, nm in (("smass2", smass_zf, smass_names),
                        ("mithneos", mith_zf, mith_names),
                        ("ecas", ecas_zf, ecas_names)):
        versions[src] = bundle_version(zf, nm)
        print(f"[{src}] bundle version_id: {versions[src]}")

    long_rows = ingest_smass2(smass_zf, smass_names)
    long_rows += ingest_mithneos(mith_zf, mith_names)
    long_df = pd.DataFrame(long_rows, columns=[
        "number_mp", "designation", "gt_source", "gt_spectrum_id",
        "wavelength_um", "reflectance", "reflectance_err"])
    long_df["number_mp"] = long_df["number_mp"].astype("Int64")

    pad = long_df["wavelength_um"] <= 0
    if pad.any():
        print(f"[gt] dropped {int(pad.sum())} padding rows (wavelength 0.0 = "
              "fill value in archived files, not data)")
        long_df = long_df[~pad].reset_index(drop=True)

    n_bad_refl = int(((long_df["reflectance"] <= 0) |
                      (long_df["reflectance"] >= 5)).sum())
    if n_bad_refl:
        print(f"[gt] NOTE {n_bad_refl:,} points with reflectance outside (0, 5) "
              "kept as archived (no cleaning in this package)")

    ecas_df = ingest_ecas(ecas_zf, ecas_names)
    ecas_df["number_mp"] = ecas_df["number_mp"].astype("Int64")
    ecas_df["gt_source"] = "ecas"
    ecas_df["gt_taxon"] = None

    summary = pd.concat([
        spectra_summary(long_df, "smass2"),
        spectra_summary(long_df, "mithneos"),
        ecas_df,
    ], ignore_index=True)
    col_order = (["number_mp", "designation", "gt_source", "gt_taxon",
                  "gt_wl_min_um", "gt_wl_max_um", "gt_n_points", "gt_n_spectra"]
                 + [c for c in summary.columns if c.startswith("ecas_")])
    summary = summary[col_order]

    sum_out = INTERIM / "groundtruth_summary.parquet"
    spec_out = INTERIM / "groundtruth_spectra.parquet"
    summary.to_parquet(sum_out, engine="pyarrow", compression="snappy", index=False)
    long_df.to_parquet(spec_out, engine="pyarrow", compression="snappy", index=False)

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    counts = {
        "smass2": int((summary["gt_source"] == "smass2").sum()),
        "mithneos": int((summary["gt_source"] == "mithneos").sum()),
        "ecas": int((summary["gt_source"] == "ecas").sum()),
    }
    for src, cfg in SOURCES.items():
        PROVENANCE.joinpath(f"{src}.json").write_text(json.dumps({
            "source": f"PDS4 bundle {cfg['urn']}",
            "url": cfg["zip_url"],
            "urn": cfg["urn"],
            "version": versions[src],
            "doi": cfg["doi"],
            "date": retrieved,
            "sha256": sha256_file(zips[src]),
            "row_count": counts[src],
            "spectra_points": int((long_df["gt_source"] == src).sum()) or None,
            "taxonomy_note": ("bundle distributes no per-object taxonomic classes; "
                              "gt_taxon left null (classes live in the separate PDS "
                              "Asteroid Taxonomy dataset)"),
            "citation": cfg["citation"] + f", retrieved {retrieved}",
            "license_note": "NASA PDS open data; cite the bundle DOI when redistributing derived tables.",
        }, indent=2))
    print("[provenance] wrote smass2.json, mithneos.json, ecas.json")

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")

    def numbers(src):
        return set(summary.loc[(summary["gt_source"] == src) &
                               summary["number_mp"].notna(), "number_mp"])

    smass_n, mith_n, ecas_n = numbers("smass2"), numbers("mithneos"), numbers("ecas")
    union_n = smass_n | mith_n | ecas_n
    g_union = int(gasp["number_mp"].isin(union_n).sum())
    visnir = smass_n & mith_n
    g_visnir = int(gasp["number_mp"].isin(visnir).sum())
    n_taxon = int(summary["gt_taxon"].notna().sum())

    apo_m = summary[(summary["gt_source"] == "mithneos") & (summary["number_mp"] == 99942)]
    apo_s = (summary.loc[summary["gt_source"] == "smass2", "number_mp"] == 99942).any()
    apo_e = (summary.loc[summary["gt_source"] == "ecas", "number_mp"] == 99942).any()
    if len(apo_m):
        r = apo_m.iloc[0]
        apo_line = (f"PRESENT in MITHNEOS ({int(r['gt_n_spectra'])} spectra, "
                    f"{int(r['gt_n_points'])} points, "
                    f"{r['gt_wl_min_um']:.3f}-{r['gt_wl_max_um']:.3f} um); "
                    f"SMASSII: {'present' if apo_s else 'absent'}, "
                    f"ECAS: {'present' if apo_e else 'absent'} (both expected absent)")
    else:
        apo_line = "ABSENT in MITHNEOS - unexpected, investigate"

    total_pts = {s: int((long_df["gt_source"] == s).sum()) for s in ("smass2", "mithneos")}
    print()
    print("=== Coverage report (AP_A5) ===")
    print(f"SMASSII objects / spectra points:   {counts['smass2']:,} / {total_pts['smass2']:,}")
    print(f"MITHNEOS objects / spectra points:  {counts['mithneos']:,} / {total_pts['mithneos']:,}")
    print(f"ECAS objects:                       {counts['ecas']:,}")
    union_all_keys = summary.assign(
        key=summary["number_mp"].astype("string").fillna(summary["designation"]))
    print(f"Any ground truth (union):           {union_all_keys['key'].nunique():,}")
    print(f"With taxon label:                   {n_taxon:,} "
          "(no bundle distributes per-object classes; see provenance)")
    print(f"GASP core matched (union):          {g_union:,} / {len(gasp):,} "
          f"({100*g_union/len(gasp):.2f}%)")
    print(f"GASP core with VIS+NIR overlap:     {g_visnir:,}  "
          f"(SMASS visible + MITHNEOS NIR; full overlap set: {len(visnir):,} objects)")
    print(f"Apophis (99942) check:              {apo_line}")
    print(f"Saved: data/interim/groundtruth_summary.parquet  {sum_out.stat().st_size/1e6:.2f} MB")
    print(f"Saved: data/interim/groundtruth_spectra.parquet  {spec_out.stat().st_size/1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
