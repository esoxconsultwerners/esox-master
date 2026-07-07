#!/usr/bin/env python3
"""AP_A3: Nesvorny HCM asteroid families.

Source: PDS4 bundle "Nesvorny HCM Asteroid Families" V2.0 (2024) at the
NASA PDS Small Bodies Node - newer than the 2015 release used in GASP v1.
The whole bundle zip is downloaded once (static, versioned, citable) and
all tables are read directly from the cached zip. PDS4 XML labels are
parsed for every table (fixed-width field positions for the 2015 .tab
files, delimited field names for the 2024 .csv files) - column layouts
are never guessed.

fam_id scheme (the 2024 products carry no family numbers):
  - families_2015: FAMILY_NUMBER from the product (3-digit HCM id)
  - families_2024: 1000000 + parent asteroid number; the two families
    with unnumbered parents get sequential ids from 1900000
    (deterministic: sorted by file name)
fam_membership: '2015' (established list) or '2024' (new families of
the V2.0 release); the product has no per-row core/halo flag.
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

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "nesvorny_families"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"

ZIP_URL = "https://sbnarchive.psi.edu/pds4/non_mission/ast.nesvorny.families_V2_0.zip"
RESOURCE_PAGE = "https://sbn.psi.edu/pds/resource/nesvornyfam.html"
PDS_URN = "urn:nasa:pds:ast.nesvorny.families::2.0"
PDS_DOI = "10.26033/5hyq-6k90"

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
    dest = RAW / f"ast.nesvorny.families_V2_0_{TODAY}.zip"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[families] using cached {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    print(f"[families] downloading {ZIP_URL}")
    tmp = dest.with_suffix(".part")
    with requests.get(ZIP_URL, headers={"User-Agent": "esox-master-catalog/A3"},
                      stream=True, timeout=(30, 1800)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[families] saved {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


def strip_ns(tag):
    return tag.split("}")[-1]


def parse_char_label(xml_text):
    """Fixed-width Table_Character label: [(name, loc_1based, length)], records."""
    root = ET.fromstring(xml_text)
    fields, records = [], None
    for el in root.iter():
        tag = strip_ns(el.tag)
        if tag == "Table_Character":
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
    if not fields or records is None:
        raise ValueError("no Table_Character spec found in label")
    return fields, records


def parse_delim_label(xml_text):
    """Table_Delimited label: field names in order, records."""
    root = ET.fromstring(xml_text)
    names, records = [], None
    for el in root.iter():
        tag = strip_ns(el.tag)
        if tag == "Table_Delimited":
            for sub in el.iter():
                st = strip_ns(sub.tag)
                if st == "records":
                    records = int(sub.text)
                elif st == "Field_Delimited":
                    for f in sub:
                        if strip_ns(f.tag) == "name":
                            names.append(f.text.strip())
    if not names or records is None:
        raise ValueError("no Table_Delimited spec found in label")
    return names, records


def read_char_table(text, fields):
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        rows.append({n: line[loc - 1: loc - 1 + ln].strip()
                     for n, loc, ln in fields})
    return rows


def bundle_version(zf, names):
    for n in names:
        if re.search(r"/bundle_[^/]*\.xml$", n):
            root = ET.fromstring(zf.read(n))
            for el in root.iter():
                if strip_ns(el.tag) == "version_id":
                    return el.text.strip()
    return None


def new_family_names(zf, names):
    """Map parent number (or compacted designation) -> display name, from
    document/list_of_new_families_2024.txt."""
    doc = next((n for n in names if n.endswith("list_of_new_families_2024.txt")), None)
    by_number, by_compact = {}, {}
    if not doc:
        return by_number, by_compact
    for line in zf.read(doc).decode("utf-8", "replace").splitlines():
        m = re.match(r"^(\S.*?)\s{2,}(\d+)\s+(\d+)\s*$", line.rstrip())
        if not m or m.group(1).startswith("-"):
            continue
        name = m.group(1).strip()
        num_m = re.match(r"^(\d+)\s+(.+)$", name)
        if num_m:
            by_number[int(num_m.group(1))] = name
        by_compact[re.sub(r"[^a-z0-9]", "", name.lower())] = name
    return by_number, by_compact


def ingest_2015(zf, names):
    rows = []
    labels = sorted(n for n in names
                    if "/data/families_2015/" in n and n.endswith(".xml"))
    for label_name in labels:
        tab_name = label_name[:-4] + ".tab"
        if tab_name not in names:
            continue
        fields, records = parse_char_label(zf.read(label_name))
        recs = read_char_table(zf.read(tab_name).decode("utf-8", "replace"), fields)
        if len(recs) != records:
            print(f"[families] WARNING {Path(tab_name).name}: "
                  f"{len(recs)} rows, label declares {records}")
        for r in recs:
            num = int(r["AST_NUMBER"])
            rows.append((num if num > 0 else None, None,
                         int(r["FAMILY_NUMBER"]), r["FAMILY_NAME"], "2015", None))
    print(f"[families] families_2015: {len(labels)} tables, {len(rows):,} members")
    return rows


def ingest_2024(zf, names, by_number, by_compact):
    rows = []
    labels = sorted(n for n in names
                    if "/data/families_2024/" in n and n.endswith(".xml"))
    zero_parent_seq = 1900000
    for label_name in labels:
        csv_name = label_name[:-4] + ".csv"
        if csv_name not in names:
            continue
        fields, records = parse_delim_label(zf.read(label_name))
        stem = Path(csv_name).stem
        m = re.match(r"^(highi|inner|middle|outer|pristine)_(\d+)_(.+?)_fam3$", stem)
        if not m:
            raise ValueError(f"unexpected 2024 family file name: {stem}")
        zone, parent, name_token = m.group(1), int(m.group(2)), m.group(3)
        if parent > 0:
            fam_id = 1000000 + parent
            fam_name = by_number.get(parent, f"{parent} {name_token.title()}")
        else:
            fam_id = zero_parent_seq
            zero_parent_seq += 1
            fam_name = by_compact.get(re.sub(r"[^a-z0-9]", "", name_token),
                                      name_token.upper())
        text = zf.read(csv_name).decode("utf-8", "replace")
        recs = list(csv.reader(io.StringIO(text)))
        recs = [r for r in recs if r and any(c.strip() for c in r)]
        if len(recs) != records:
            print(f"[families] WARNING {Path(csv_name).name}: "
                  f"{len(recs)} rows, label declares {records}")
        idx = {f: i for i, f in enumerate(fields)}
        i_name = idx["UNPACKED_NAME"]
        for r in recs:
            unpacked = r[i_name].strip()
            if unpacked.isdigit():
                rows.append((int(unpacked), None, fam_id, fam_name, "2024", zone))
            else:
                rows.append((None, unpacked, fam_id, fam_name, "2024", zone))
    print(f"[families] families_2024: {len(labels)} tables, "
          f"{len(rows):,} members")
    return rows


def main():
    for d in (RAW, INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)

    bundle = download_bundle()
    zf = zipfile.ZipFile(bundle)
    names = set(zf.namelist())

    version = bundle_version(zf, names)
    print(f"[families] bundle version_id: {version}")

    by_number, by_compact = new_family_names(zf, names)
    rows = ingest_2015(zf, names) + ingest_2024(zf, names, by_number, by_compact)

    df = pd.DataFrame(rows, columns=[
        "number_mp", "designation", "fam_id", "fam_name",
        "fam_membership", "fam_zone"])
    df["number_mp"] = df["number_mp"].astype("Int64")
    df["fam_id"] = df["fam_id"].astype("int64")

    before = len(df)
    df = df.drop_duplicates(subset=["number_mp", "designation", "fam_id"])
    if len(df) != before:
        print(f"[families] dropped {before - len(df)} exact duplicate memberships")

    out = INTERIM / "families.parquet"
    df.to_parquet(out, engine="pyarrow", compression="snappy", index=False)

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    prov = {
        "source": "Nesvorny HCM Asteroid Families, PDS4 bundle V2.0, NASA PDS Small Bodies Node",
        "url": ZIP_URL,
        "resource_page": RESOURCE_PAGE,
        "urn": PDS_URN,
        "version": version,
        "doi": PDS_DOI,
        "date": retrieved,
        "sha256": sha256_file(bundle),
        "row_count": int(len(df)),
        "families_2015_members": int((df["fam_membership"] == "2015").sum()),
        "families_2024_members": int((df["fam_membership"] == "2024").sum()),
        "distinct_families": int(df["fam_id"].nunique()),
        "fam_id_note": ("2015 families keep their HCM FAMILY_NUMBER; 2024 families "
                        "(no id in product) use 1000000 + parent asteroid number, "
                        "unnumbered parents get sequential ids from 1900000"),
        "citation": ("Nesvorny, D., Broz, M., Carruba, V. (2015), Identification and "
                     "Dynamical Properties of Asteroid Families, Asteroids IV, "
                     "Univ. of Arizona Press; via Nesvorny, D. (2024), Nesvorny HCM "
                     f"Asteroid Families V2.0, {PDS_URN}, NASA Planetary Data System, "
                     f"doi:{PDS_DOI}, retrieved {retrieved}"),
        "license_note": "NASA PDS open data; cite the bundle DOI when redistributing derived tables.",
    }
    PROVENANCE.joinpath("families.json").write_text(json.dumps(prov, indent=2))
    print("[provenance] wrote data/provenance/families.json")

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    member_numbers = df.loc[df["number_mp"].notna(), "number_mp"].unique()
    gasp_matched = int(gasp["number_mp"].isin(member_numbers).sum())
    gasp_rate = 100 * gasp_matched / len(gasp)

    n_numbered = int(df["number_mp"].notna().sum())
    top5 = (df.groupby(["fam_id", "fam_name"]).size()
              .sort_values(ascending=False).head(5))

    apo = df[df["number_mp"] == 99942]
    if len(apo):
        r = apo.iloc[0]
        apo_line = f"PRESENT: fam_id {r['fam_id']} ({r['fam_name']}) - unexpected, see notes"
    else:
        apo_line = "absent as expected (NEO; families are main-belt)"

    size_mb = out.stat().st_size / 1e6
    print()
    print("=== Coverage report (AP_A3) ===")
    print(f"Total family members:          {len(df):,}")
    print(f"Distinct families:             {df['fam_id'].nunique():,}")
    print("Largest families:")
    for (fid, fname), cnt in top5.items():
        print(f"    {fname:<22} (fam_id {fid}): {cnt:,}")
    print(f"Numbered / unnumbered:         {n_numbered:,} / {len(df) - n_numbered:,}")
    print(f"GASP core matched:             {gasp_matched:,} / {len(gasp):,} ({gasp_rate:.2f}%)")
    print(f"    (GASP v1 with 2015 lists: 39.40% -> delta {gasp_rate - 39.40:+.2f} pp)")
    print(f"Apophis (99942) check:         {apo_line}")
    print(f"Saved: data/interim/families.parquet  {size_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
