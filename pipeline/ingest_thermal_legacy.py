#!/usr/bin/env python3
"""AP_A4: AKARI + IRAS/SIMPS thermal diameters and albedos.

Two independent legacy thermal surveys, kept as separate column groups
(akari_*, iras_*). Deliberately NO fusion into a "best" albedo here -
precedence rules are a Phase-B decision.

Source 1: AKARI Asteroid Catalog (AcuA) V1, Usui et al. 2011, PASJ 63,
1117. Retrieved from JAXA DARTS (primary distribution). The fixed-width
layout is parsed from the byte-by-byte block in ReadMe.AcuA.txt, never
guessed.

Source 2: IRAS Minor Planet Survey (SIMPS) V6.0, Tedesco et al. 2002,
AJ 123, 1056. PDS3 dataset IRAS-A-FPA-3-RDR-IMPS-V6.0 at the SBN;
product DIAMALB. Column positions are parsed from the PDS3 label.

Values are reported as archived (e.g. the famous IRAS Ceres ~848 km vs
modern ~940 km) and never "corrected".
"""

import gzip
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"

AKARI_BASE = "https://data.darts.isas.jaxa.jp/pub/akari/AKARI-IRC_Catalogue_AllSky_AcuA_1.0"
AKARI_DATA_URL = f"{AKARI_BASE}/AcuA_V1.txt.gz"
AKARI_README_URL = f"{AKARI_BASE}/ReadMe.AcuA.txt"

SIMPS_BASE = "https://sbnarchive.psi.edu/pds3/iras/IRAS_A_FPA_3_RDR_IMPS_V6_0/data"
SIMPS_TAB_URL = f"{SIMPS_BASE}/diamalb.tab"
SIMPS_LBL_URL = f"{SIMPS_BASE}/diamalb.lbl"
SIMPS_DATASET_ID = "IRAS-A-FPA-3-RDR-IMPS-V6.0"
SIMPS_DOI = "10.26033/pf3k-m168"

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
        print(f"[thermal] using cached {dest.name}")
        return dest
    print(f"[thermal] downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, headers={"User-Agent": "esox-master-catalog/A4"},
                      stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[thermal] saved {dest.name} ({dest.stat().st_size/1e3:.0f} kB)")
    return dest


def to_num(s):
    s = s.strip()
    if not s:
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def parse_akari_readme(readme_path):
    """Parse the CDS-style byte-by-byte block: [(label, start, end)] 1-based."""
    fields = []
    in_block = False
    for line in readme_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Byte-by-byte Description" in line:
            in_block = True
            continue
        if in_block:
            m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s+\S+\s+\S+\s+(\S+)", line)
            if m:
                fields.append((m.group(3), int(m.group(1)), int(m.group(2))))
            elif fields and re.match(r"^-{20,}", line):
                break
    if not fields:
        raise ValueError("could not parse byte-by-byte block from AKARI ReadMe")
    return fields


def ingest_akari():
    raw_dir = RAW / "akari"
    raw_dir.mkdir(parents=True, exist_ok=True)
    data_file = download(AKARI_DATA_URL, raw_dir / f"AcuA_V1_{TODAY}.txt.gz")
    readme_file = download(AKARI_README_URL, raw_dir / f"ReadMe.AcuA_{TODAY}.txt")

    fields = parse_akari_readme(readme_file)
    print(f"[akari] format spec: {[(n, a, b) for n, a, b in fields]}")
    need = {"NUMBER", "PROV_DES", "NID", "DIAMETER", "D_ERR", "ALBEDO", "A_ERR"}
    missing = need - {n for n, _, _ in fields}
    if missing:
        raise ValueError(f"AKARI ReadMe lacks expected labels: {missing}")

    rows = []
    with gzip.open(data_file, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            rec = {n: line[a - 1:b].strip() for n, a, b in fields}
            num = rec["NUMBER"]
            rows.append((
                int(num) if num.isdigit() else None,
                rec["PROV_DES"] or None,
                rec["NAME"] or None,
                to_num(rec["DIAMETER"]),
                to_num(rec["D_ERR"]),
                to_num(rec["ALBEDO"]),
                to_num(rec["A_ERR"]),
                int(rec["NID"]) if rec["NID"].isdigit() else None,
            ))
    df = pd.DataFrame(rows, columns=[
        "number_mp", "designation", "akari_name",
        "akari_diameter_km", "akari_diameter_err",
        "akari_albedo", "akari_albedo_err", "akari_n_obs"])
    df["number_mp"] = df["number_mp"].astype("Int64")
    df["akari_n_obs"] = df["akari_n_obs"].astype("Int64")
    print(f"[akari] parsed {len(df):,} rows")
    return df, data_file, readme_file


def parse_pds3_columns(lbl_path):
    """Parse OBJECT=COLUMN blocks from a PDS3 label: [(name, start, bytes)]."""
    cols, cur = [], None
    for raw in lbl_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if re.match(r"^OBJECT\s*=\s*COLUMN$", line):
            cur = {}
        elif cur is not None:
            m = re.match(r"^(\w+)\s*=\s*\"?([^\"]+?)\"?\s*$", line)
            if m:
                cur[m.group(1)] = m.group(2)
            if re.match(r"^END_OBJECT\s*=\s*COLUMN$", line):
                cols.append((cur["NAME"], int(cur["START_BYTE"]), int(cur["BYTES"])))
                cur = None
    if not cols:
        raise ValueError("no COLUMN objects found in PDS3 label")
    return cols


def ingest_simps():
    raw_dir = RAW / "iras_simps"
    raw_dir.mkdir(parents=True, exist_ok=True)
    tab_file = download(SIMPS_TAB_URL, raw_dir / f"diamalb_{TODAY}.tab")
    lbl_file = download(SIMPS_LBL_URL, raw_dir / f"diamalb_{TODAY}.lbl")

    cols = parse_pds3_columns(lbl_file)
    print(f"[iras] label: {len(cols)} columns, using AST_NUMBER, AST_DESIG, "
          "MEAN_ALBEDO, ALBEDO_UNC, DIAM, DIAM_UNC")
    need = {"AST_NUMBER", "AST_DESIG", "MEAN_ALBEDO", "ALBEDO_UNC", "DIAM", "DIAM_UNC"}
    missing = need - {n for n, _, _ in cols}
    if missing:
        raise ValueError(f"SIMPS label lacks expected columns: {missing}")
    pos = {n: (s, b) for n, s, b in cols}

    def cut(line, name):
        s, b = pos[name]
        return line[s - 1:s - 1 + b].strip()

    rows = []
    for line in tab_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        num = cut(line, "AST_NUMBER")
        rows.append((
            int(num) if num.isdigit() else None,
            cut(line, "AST_DESIG") or None,
            to_num(cut(line, "DIAM")),
            to_num(cut(line, "DIAM_UNC")),
            to_num(cut(line, "MEAN_ALBEDO")),
            to_num(cut(line, "ALBEDO_UNC")),
        ))
    df = pd.DataFrame(rows, columns=[
        "number_mp", "designation",
        "iras_diameter_km", "iras_diameter_err",
        "iras_albedo", "iras_albedo_err"])
    df["number_mp"] = df["number_mp"].astype("Int64")
    print(f"[iras] parsed {len(df):,} rows")
    return df, tab_file, lbl_file


def main():
    for d in (INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)

    akari, akari_data, akari_readme = ingest_akari()
    iras, iras_tab, iras_lbl = ingest_simps()

    for name, df in (("akari", akari), ("iras", iras)):
        numbered = df.loc[df["number_mp"].notna(), "number_mp"]
        n_dup = int(numbered.duplicated().sum())
        if n_dup:
            print(f"[{name}] WARNING {n_dup} duplicate number_mp rows")

    akari_out = INTERIM / "akari.parquet"
    iras_out = INTERIM / "iras_simps.parquet"
    akari.to_parquet(akari_out, engine="pyarrow", compression="snappy", index=False)
    iras.to_parquet(iras_out, engine="pyarrow", compression="snappy", index=False)

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    PROVENANCE.joinpath("akari.json").write_text(json.dumps({
        "source": "AKARI Asteroid Catalog (AcuA) V1.0, AKARI/IRC Mid-Infrared All-Sky Survey",
        "source_used": "JAXA DARTS (primary distribution)",
        "url": AKARI_DATA_URL,
        "readme_url": AKARI_README_URL,
        "version": "AcuA 1.0",
        "bibcode": "2011PASJ...63.1117U",
        "date": retrieved,
        "sha256": sha256_file(akari_data),
        "readme_sha256": sha256_file(akari_readme),
        "row_count": int(len(akari)),
        "citation": ("Usui, F., et al. (2011), Asteroid Catalog Using AKARI: AKARI/IRC "
                     "Mid-Infrared Asteroid Survey, PASJ 63, 1117-1138; "
                     f"AcuA V1.0 via JAXA DARTS, retrieved {retrieved}"),
        "license_note": "JAXA/ISAS DARTS public archive data.",
    }, indent=2))
    PROVENANCE.joinpath("iras_simps.json").write_text(json.dumps({
        "source": "Supplemental IRAS Minor Planet Survey (SIMPS), DIAMALB product",
        "url": SIMPS_TAB_URL,
        "label_url": SIMPS_LBL_URL,
        "dataset_id": SIMPS_DATASET_ID,
        "version": "V6.0",
        "doi": SIMPS_DOI,
        "date": retrieved,
        "sha256": sha256_file(iras_tab),
        "label_sha256": sha256_file(iras_lbl),
        "row_count": int(len(iras)),
        "citation": ("Tedesco, E.F., Noah, P.V., Noah, M., Price, S.D. (2002), The "
                     "Supplemental IRAS Minor Planet Survey, AJ 123, 1056-1085; "
                     f"{SIMPS_DATASET_ID}, NASA Planetary Data System, 2004, "
                     f"doi:{SIMPS_DOI}, retrieved {retrieved}"),
        "license_note": "NASA PDS open data; cite the dataset DOI when redistributing derived tables.",
    }, indent=2))
    print("[provenance] wrote data/provenance/akari.json and iras_simps.json")

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    akari_nums = set(akari.loc[akari["number_mp"].notna(), "number_mp"])
    iras_nums = set(iras.loc[iras["number_mp"].notna(), "number_mp"])
    g_akari = int(gasp["number_mp"].isin(akari_nums).sum())
    g_iras = int(gasp["number_mp"].isin(iras_nums).sum())
    g_union = int(gasp["number_mp"].isin(akari_nums | iras_nums).sum())
    overlap = len(akari_nums & iras_nums)

    def one(df, col, n=1):
        row = df[df["number_mp"] == n]
        return row.iloc[0][col] if len(row) else None

    apo_a = (akari["number_mp"] == 99942).any()
    apo_i = (iras["number_mp"] == 99942).any()
    apo_line = ("absent in both as expected (too small/faint for both surveys)"
                if not (apo_a or apo_i)
                else f"UNEXPECTED: AKARI={apo_a}, IRAS={apo_i}")

    a_mb = akari_out.stat().st_size / 1e6
    i_mb = iras_out.stat().st_size / 1e6
    print()
    print("=== Coverage report (AP_A4) ===")
    print(f"AKARI rows:                    {len(akari):,}")
    print(f"IRAS/SIMPS rows:               {len(iras):,}")
    print(f"Overlap AKARI x IRAS:          {overlap:,}")
    print(f"GASP core matched (AKARI):     {g_akari:,} / {len(gasp):,} ({100*g_akari/len(gasp):.2f}%)")
    print(f"GASP core matched (IRAS):      {g_iras:,} / {len(gasp):,} ({100*g_iras/len(gasp):.2f}%)")
    print(f"Union with any thermal D:      AKARI or IRAS: {g_union:,} ({100*g_union/len(gasp):.2f}%)")
    print(f"Apophis (99942) check:         {apo_line}")
    print(f"Ceres (1) cross-check:         akari_D {one(akari, 'akari_diameter_km')} km, "
          f"iras_D {one(iras, 'iras_diameter_km')} km")
    print("    NOTE: IRAS/SIMPS lists Ceres near ~848 km, well below the modern")
    print("    ~940 km value - known literature feature, reported as archived.")
    print(f"Saved: data/interim/akari.parquet  {a_mb:.2f} MB")
    print(f"Saved: data/interim/iras_simps.parquet  {i_mb:.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
