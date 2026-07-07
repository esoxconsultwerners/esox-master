#!/usr/bin/env python3
"""AP_A2: LCDB rotation periods.

Source: Asteroid Lightcurve Database (LCDB), Warner/Harris/Pravec,
archived as versioned PDS4 bundle at the NASA PDS Small Bodies Node.
Product: lc_summary.csv (one row per asteroid) + its PDS4 XML label,
which is the authoritative format spec (field names/order are parsed
from the label, never guessed).

Writes data/interim/lcdb.parquet plus provenance JSON and prints a
coverage report against the GASP core, including the mandatory
Apophis (99942) showcase line.
"""

import csv
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "lcdb"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"

BUNDLE_BASE = "https://sbnarchive.psi.edu/pds4/non_mission/ast-lightcurve-database_V4_0"
CSV_URL = f"{BUNDLE_BASE}/data/lc_summary.csv"
LABEL_URL = f"{BUNDLE_BASE}/data/lc_summary.xml"
RESOURCE_PAGE = "https://sbn.psi.edu/pds/resource/lc.html"

PDS_URN = "urn:nasa:pds:ast-lightcurve-database::4.0"
PDS_DOI = "10.26033/j3xc-3359"

NUMERIC_SENTINEL = -9.99

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
        print(f"[lcdb] using cached {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    print(f"[lcdb] downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, headers={"User-Agent": "esox-master-catalog/A2"},
                      stream=True, timeout=(30, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[lcdb] saved {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


def strip_ns(tag):
    return tag.split("}")[-1]


def parse_label(label_path):
    """Parse the PDS4 XML label: field names in order, table byte offset,
    declared record count and product version."""
    root = ET.parse(label_path).getroot()
    fields, offset, records, version = [], None, None, None
    for el in root.iter():
        tag = strip_ns(el.tag)
        if tag == "version_id" and version is None:
            version = el.text.strip()
        elif tag == "Table_Delimited":
            for sub in el.iter():
                st = strip_ns(sub.tag)
                if st == "offset":
                    offset = int(sub.text)
                elif st == "records":
                    records = int(sub.text)
                elif st == "Field_Delimited":
                    for f in sub:
                        if strip_ns(f.tag) == "name":
                            fields.append(f.text.strip())
    if not fields or offset is None or records is None:
        raise ValueError(f"could not parse table spec from {label_path.name}")
    print(f"[lcdb] label: version {version}, {len(fields)} fields, "
          f"table at byte {offset}, {records} records declared")
    return fields, offset, records, version


def clean_str(s):
    s = (s or "").strip()
    return None if s in ("", "-") else s


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return np.nan
    return np.nan if x == NUMERIC_SENTINEL else x


def quality_numeric(u):
    """U code to float: base digit, '+' adds 0.3, '-' subtracts 0.3."""
    if not u:
        return np.nan
    base = float(u[0])
    if u.endswith("+"):
        return base + 0.3
    if u.endswith("-"):
        return base - 0.3
    return base


def parse_summary(csv_path, fields, offset, records):
    with open(csv_path, encoding="utf-8", errors="replace", newline="") as f:
        f.seek(offset)
        rows = list(csv.reader(f))
    if len(rows) != records:
        print(f"[lcdb] WARNING parsed {len(rows)} rows, label declares {records}")
    bad = [r for r in rows if len(r) != len(fields)]
    if bad:
        raise ValueError(f"{len(bad)} rows with field count != {len(fields)}")
    raw = pd.DataFrame(rows, columns=fields)

    df = pd.DataFrame()
    number = pd.to_numeric(raw["Number"], errors="coerce")
    df["number_mp"] = number.where(number > 0).astype("Int64")
    df["designation"] = raw["Desig"].map(clean_str)
    df["lcdb_name"] = raw["Name"].map(clean_str)
    df["lcdb_period_h"] = raw["Period"].map(num)
    df["lcdb_period_flag"] = raw["PFlag"].map(clean_str)
    df["lcdb_quality_u"] = raw["U"].map(clean_str)
    df["lcdb_quality_numeric"] = df["lcdb_quality_u"].map(quality_numeric)
    df["lcdb_amp_min"] = raw["AmpMin"].map(num)
    df["lcdb_amp_max"] = raw["AmpMax"].map(num)
    df["lcdb_h"] = raw["H"].map(num)
    df["lcdb_albedo"] = raw["Albedo"].map(num)
    df["lcdb_diameter_km"] = raw["Diam"].map(num)
    df["lcdb_class"] = raw["Class"].map(clean_str)
    df["lcdb_binary_flag"] = raw["IsBinary"].map(clean_str)
    df["lcdb_notes"] = raw["Notes"].map(clean_str)

    neg_period = int((df["lcdb_period_h"] <= 0).sum())
    if neg_period:
        print(f"[lcdb] WARNING {neg_period} non-positive periods set to null")
        df.loc[df["lcdb_period_h"] <= 0, "lcdb_period_h"] = np.nan
    return df


def dedup(df):
    """LCDB summary is nominally one row per object; if duplicates exist
    among numbered objects, keep the row with the highest
    lcdb_quality_numeric (nulls last)."""
    numbered = df["number_mp"].notna()
    n_dup = int(df.loc[numbered, "number_mp"].duplicated().sum())
    if n_dup:
        print(f"[lcdb] deduplicating {n_dup} duplicate numbered rows "
              f"(keeping highest lcdb_quality_numeric)")
        df = df.sort_values("lcdb_quality_numeric", ascending=False, na_position="last")
        keep = ~(df["number_mp"].duplicated() & df["number_mp"].notna())
        df = df[keep].sort_index()
    return df


def main():
    for d in (RAW, INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)

    csv_file = download(CSV_URL, RAW / f"lc_summary_{TODAY}.csv")
    label_file = download(LABEL_URL, RAW / f"lc_summary_{TODAY}.xml")

    fields, offset, records, version = parse_label(label_file)
    df = parse_summary(csv_file, fields, offset, records)
    df = dedup(df)

    out = INTERIM / "lcdb.parquet"
    df.to_parquet(out, engine="pyarrow", compression="snappy", index=False)

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    prov = {
        "source": "Asteroid Lightcurve Database (LCDB), PDS4 bundle V4.0, NASA PDS Small Bodies Node",
        "url": CSV_URL,
        "label_url": LABEL_URL,
        "resource_page": RESOURCE_PAGE,
        "version": version,
        "pds_urn": PDS_URN,
        "doi": PDS_DOI,
        "date": retrieved,
        "sha256": sha256_file(csv_file),
        "label_sha256": sha256_file(label_file),
        "row_count": int(len(df)),
        "rows_with_period": int(df["lcdb_period_h"].notna().sum()),
        "citation": ("Warner, B.D., Harris, A.W., Pravec, P., Icarus 202, 134-146 (2009); "
                     "updated PDS4 bundle V4.0 (data through 2021-09-02), "
                     f"{PDS_URN}, doi:{PDS_DOI}, retrieved {retrieved}"),
        "license_note": "NASA PDS open data; cite the bundle DOI when redistributing derived tables.",
    }
    PROVENANCE.joinpath("lcdb.json").write_text(json.dumps(prov, indent=2))
    print(f"[provenance] wrote data/provenance/lcdb.json")

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    numbered_keys = df.loc[df["number_mp"].notna(), "number_mp"]
    gasp_hit = gasp["number_mp"].isin(numbered_keys)
    gasp_matched = int(gasp_hit.sum())
    reliable = df["lcdb_quality_numeric"] >= 2.0
    reliable_keys = df.loc[reliable & df["number_mp"].notna(), "number_mp"]
    gasp_reliable = int(gasp["number_mp"].isin(reliable_keys).sum())

    n_period = int(df["lcdb_period_h"].notna().sum())
    n_numbered = int(df["number_mp"].notna().sum())

    apo = df[df["number_mp"] == 99942]
    if len(apo):
        r = apo.iloc[0]
        apo_line = (f"period {r['lcdb_period_h']} h, U {r['lcdb_quality_u']}, "
                    f"amp {r['lcdb_amp_min']}-{r['lcdb_amp_max']}"
                    + (f", notes {r['lcdb_notes']} (T = tumbler/NPA)" if r["lcdb_notes"] else ""))
    else:
        apo_line = "NOT FOUND in LCDB summary"

    size_mb = out.stat().st_size / 1e6
    print()
    print("=== Coverage report (AP_A2) ===")
    print(f"Total LCDB rows:               {len(df):,}")
    print(f"With period:                   {n_period:,}")
    print(f"Reliable (U >= 2):             {int(reliable.sum()):,}")
    print(f"Numbered / unnumbered:         {n_numbered:,} / {len(df) - n_numbered:,}")
    print(f"GASP core matched:             {gasp_matched:,} / {len(gasp):,} "
          f"({100*gasp_matched/len(gasp):.2f}%)")
    print(f"GASP core with reliable P:     {gasp_reliable:,} ({100*gasp_reliable/len(gasp):.2f}%)")
    print(f"Apophis (99942) check:         {apo_line}")
    print(f"Saved: data/interim/lcdb.parquet  {size_mb:.1f} MB")

    print()
    print("Spot checks vs literature:")
    for n, lit in ((433, "5.27 h"), (216, "5.385 h")):
        row = df[df["number_mp"] == n]
        if len(row):
            r = row.iloc[0]
            print(f"  ({n}) {r['lcdb_name']}: LCDB {r['lcdb_period_h']} h "
                  f"(U {r['lcdb_quality_u']}) | literature {lit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
