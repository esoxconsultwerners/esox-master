#!/usr/bin/env python3
"""AP_A5b: per-object taxonomy labels from the PDS Asteroid Taxonomy
compilation.

Source: PDS4 bundle "Asteroid Taxonomy" V1.1 (Neese, Ed.), the PDS4
migration of the PDS3 lineage EAR-A-5-DDR-TAXONOMY (V6.0). Product
taxonomy10.tab carries per-object classes from nine systems: Tholen,
Barucci, Tedesco, Howell, SMASS (Xu), Bus (SMASSII), S3OS2 (Lazzaro;
both the Tholen-like and Bus-like variants) and Bus-DeMeo.
Fixed-width layout parsed from the PDS4 label, never guessed.

Outputs: taxonomy_long.parquet (authoritative, one row per object x
system) and taxonomy_wide.parquet (convenience view). The wide
taxon_s3os2 column uses the Bus-like variant (S3OS2_CLASS_BB) because
the Phase-C target system is Bus-DeMeo; the Tholen-like variant stays
available in the long table as system s3os2_th. The source has one row
per object, so no dedup is needed; a guard drops and reports any
duplicate (object, system) pairs should a future version introduce them.

No harmonization across systems here (Phase C decision).
"""

import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "taxonomy"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"
DOCS = ROOT / "docs"

ZIP_URL = "https://sbnarchive.psi.edu/pds4/non_mission/ast_taxonomy_v1.1.zip"
PDS_URN = "urn:nasa:pds:ast_taxonomy::1.1"
PDS_DOI = "10.26033/e1p3-xm59"
GASP_CATALOG = Path("/root/gasp/data/final/gasp_catalog_v1.parquet")

SYSTEMS = [
    ("tholen", "THOLEN_CLASS", "THOLEN_PARAM"),
    ("barucci", "BARUCCI_CLASS", "BARUCCI_PARAM"),
    ("tedesco", "TEDESCO_CLASS", "TEDESCO_PARAM"),
    ("howell", "HOWELL_CLASS", "HOWELL_PARAM"),
    ("smass", "SMASS_CLASS", "SMASS_PARAM"),
    ("bus", "BUS_CLASS", "BUS_PARAM"),
    ("s3os2_th", "S3OS2_CLASS_TH", None),
    ("s3os2_bb", "S3OS2_CLASS_BB", None),
    ("bus_demeo", "BUS_DEMEO_CLASS", "DEMEO_REF_CODE"),
]

WIDE_COLUMNS = {
    "taxon_tholen": "tholen",
    "taxon_bus": "bus",
    "taxon_bus_demeo": "bus_demeo",
    "taxon_s3os2": "s3os2_bb",
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


def download_bundle():
    RAW.mkdir(parents=True, exist_ok=True)
    dest = RAW / f"ast_taxonomy_v1.1_{TODAY}.zip"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[taxonomy] using cached {dest.name}")
        return dest
    print(f"[taxonomy] downloading {ZIP_URL}")
    tmp = dest.with_suffix(".part")
    with requests.get(ZIP_URL, headers={"User-Agent": "esox-master-catalog/A5b"},
                      stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[taxonomy] saved {dest.name} ({dest.stat().st_size/1e3:.0f} kB)")
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


def bundle_version(zf, names):
    for n in names:
        if re.search(r"(^|/)bundle_[^/]*\.xml$", n):
            root = ET.fromstring(zf.read(n))
            for el in root.iter():
                if strip_ns(el.tag) == "version_id":
                    return el.text.strip()
    return None


def clean(s):
    s = (s or "").strip()
    return None if s in ("", "-") else s


def main():
    for d in (INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)

    bundle = download_bundle()
    zf = zipfile.ZipFile(bundle)
    names = set(zf.namelist())
    version = bundle_version(zf, names)
    print(f"[taxonomy] bundle version_id: {version}")

    label_name = next(n for n in names if n.endswith("data/taxonomy10.xml"))
    tab_name = label_name[:-4] + ".tab"
    fields, records = parse_char_label(zf.read(label_name))
    idx = {n: (loc, ln) for n, loc, ln in fields}
    need = {"AST_NUMBER", "AST_NAME", "PROV_ID"} | \
        {c for _, c, p in SYSTEMS for c in ([c] + ([p] if p else []))}
    missing = need - set(idx)
    if missing:
        raise ValueError(f"taxonomy10 label lacks fields: {missing}")
    print(f"[taxonomy] label: {len(fields)} fields, {records} records declared")

    long_rows = []
    n_lines = 0
    for line in zf.read(tab_name).decode("utf-8", "replace").splitlines():
        if not line.strip():
            continue
        n_lines += 1

        def cut(name):
            loc, ln = idx[name]
            return line[loc - 1:loc - 1 + ln].strip()

        num_s = cut("AST_NUMBER")
        num = int(num_s) if num_s.isdigit() and int(num_s) > 0 else None
        desig = clean(cut("PROV_ID"))
        for system, class_col, param_col in SYSTEMS:
            cls = clean(cut(class_col))
            if cls is None:
                continue
            note = clean(cut(param_col)) if param_col else None
            long_rows.append((num, desig, system, cls, note))
    if records is not None and n_lines != records:
        print(f"[taxonomy] WARNING parsed {n_lines} rows, label declares {records}")

    long_df = pd.DataFrame(long_rows, columns=[
        "number_mp", "designation", "taxon_system", "taxon_class", "taxon_note"])
    long_df["number_mp"] = long_df["number_mp"].astype("Int64")

    key = long_df["number_mp"].astype("string").fillna(long_df["designation"])
    dup = long_df.assign(key=key).duplicated(subset=["key", "taxon_system"])
    if dup.any():
        print(f"[taxonomy] dropping {int(dup.sum())} duplicate (object, system) rows "
              "(keeping first occurrence)")
        long_df = long_df[~dup.values].reset_index(drop=True)

    wide = long_df.assign(key=key[~dup.values] if dup.any() else key)
    base = wide[["key", "number_mp", "designation"]].drop_duplicates("key")
    wide_df = base.copy()
    for col, system in WIDE_COLUMNS.items():
        sub = wide.loc[wide["taxon_system"] == system, ["key", "taxon_class"]]
        wide_df = wide_df.merge(sub.rename(columns={"taxon_class": col}),
                                on="key", how="left")
    wide_df = wide_df.drop(columns=["key"])

    long_out = INTERIM / "taxonomy_long.parquet"
    wide_out = INTERIM / "taxonomy_wide.parquet"
    long_df.to_parquet(long_out, engine="pyarrow", compression="snappy", index=False)
    wide_df.to_parquet(wide_out, engine="pyarrow", compression="snappy", index=False)

    sys_counts = long_df["taxon_system"].value_counts().to_dict()
    n_objects = key.nunique()

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    PROVENANCE.joinpath("taxonomy.json").write_text(json.dumps({
        "source": "Asteroid Taxonomy V1.1, PDS4 bundle (PDS3 lineage EAR-A-5-DDR-TAXONOMY-V6.0)",
        "url": ZIP_URL,
        "urn": PDS_URN,
        "version": version,
        "doi": PDS_DOI,
        "date": retrieved,
        "sha256": sha256_file(bundle),
        "row_count": int(len(long_df)),
        "distinct_objects": int(n_objects),
        "per_system_counts": sys_counts,
        "wide_view_rule": ("one row per object in source, no dedup needed; "
                           "taxon_s3os2 in the wide view = Bus-like variant "
                           "(S3OS2_CLASS_BB); Tholen-like variant kept in the "
                           "long table as s3os2_th"),
        "citation": ("Neese, C. (Ed.), Asteroid Taxonomy V6.0, "
                     "EAR-A-5-DDR-TAXONOMY-V6.0, NASA Planetary Data System (2010); "
                     f"PDS4 migration V1.1 (2017), {PDS_URN}, doi:{PDS_DOI}, "
                     f"retrieved {retrieved}"),
        "license_note": "NASA PDS open data; cite the bundle DOI when redistributing derived tables.",
    }, indent=2))
    print("[provenance] wrote data/provenance/taxonomy.json")

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    our_numbers = set(long_df.loc[long_df["number_mp"].notna(), "number_mp"])
    g_any = int(gasp["number_mp"].isin(our_numbers).sum())

    mahlke_line = "GASP catalog not readable - cross-check skipped"
    try:
        gcat = pd.read_parquet(GASP_CATALOG, columns=["number_mp", "taxonomy"])
        gcat = gcat[gcat["taxonomy"].notna()]
        prio = {"bus_demeo": 0, "bus": 1, "tholen": 2}
        ours = (long_df[long_df["taxon_system"].isin(prio)]
                .assign(p=lambda d: d["taxon_system"].map(prio))
                .sort_values("p")
                .drop_duplicates("number_mp"))
        both = gcat.merge(ours[["number_mp", "taxon_system", "taxon_class"]],
                          on="number_mp", how="inner")
        agree = (both["taxonomy"].str[0].str.upper()
                 == both["taxon_class"].str[0].str.upper())
        mahlke_line = (f"objects in both {len(both):,}, simple-letter agreement "
                       f"{100*agree.mean():.1f}% (vs GASP 'taxonomy' column, "
                       f"n={len(gcat):,} external labels; report only)")
    except Exception as e:
        mahlke_line += f" ({e})"

    def wide_val(num, col):
        row = wide_df[wide_df["number_mp"] == num]
        return row.iloc[0][col] if len(row) and pd.notna(row.iloc[0][col]) else "-"

    apo = long_df[long_df["number_mp"] == 99942]
    if len(apo):
        apo_line = ", ".join(f"{r.taxon_system}={r.taxon_class}"
                             for r in apo.itertuples())
    else:
        apo_line = ("absent (compilation predates wide NEO coverage; "
                    "Sq label per Binzel et al. 2009 remains literature-only)")

    print()
    print("=== Coverage report (AP_A5b) ===")
    print("Systems found:")
    for system, _, _ in SYSTEMS:
        print(f"    {system:<10} {sys_counts.get(system, 0):,}")
    print(f"Distinct objects any system:   {n_objects:,}")
    print(f"GASP core matched (any):       {g_any:,} / {len(gasp):,} "
          f"({100*g_any/len(gasp):.2f}%)")
    print(f"Cross-check vs GASP Mahlke layer: {mahlke_line}")
    print(f"Vesta (4):    tholen {wide_val(4, 'taxon_tholen')}, "
          f"bus {wide_val(4, 'taxon_bus')}     (expect V / V)")
    print(f"Ceres (1):    tholen {wide_val(1, 'taxon_tholen')}, "
          f"bus {wide_val(1, 'taxon_bus')}     (expect G / C)")
    print(f"Apophis (99942): {apo_line}")
    print(f"Saved: data/interim/taxonomy_long.parquet  "
          f"{long_out.stat().st_size/1e6:.2f} MB")
    print(f"Saved: data/interim/taxonomy_wide.parquet  "
          f"{wide_out.stat().st_size/1e6:.2f} MB")

    status = DOCS / "apophis_status.md"
    if status.exists() and "| A5b |" not in status.read_text():
        lines = status.read_text().splitlines()
        last_row = max(i for i, l in enumerate(lines) if l.startswith("|"))
        lines.insert(last_row + 1,
                     f"| A5b | PDS Asteroid Taxonomy V1.1 | {apo_line} |")
        status.write_text("\n".join(lines) + "\n")
        print("[docs] appended Apophis line to docs/apophis_status.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
