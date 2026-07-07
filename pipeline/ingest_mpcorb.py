#!/usr/bin/env python3
"""AP_A1: orbit backbone.

Spine:      MPC Extended Orbit File (mpcorb_extended.json.gz, ~1.4M objects)
Enrichment: JPL SBDB Query API (numbered asteroids, one bulk call)

Writes data/interim/orbits_backbone.parquet plus provenance JSON and prints
a coverage report against the GASP core (19,190 objects, key number_mp).
"""

import gzip
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROVENANCE = ROOT / "data" / "provenance"

MPCORB_URL = "https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz"
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"

SBDB_FIELDS = (
    "spkid,full_name,pdes,neo,pha,H,diameter,albedo,class,"
    "e,a,q,i,om,w,ma,ad,per_y,moid,moid_jup,t_jup,"
    "condition_code,n_obs_used,first_obs,last_obs"
)

TODAY = datetime.now(timezone.utc).strftime("%Y%m%d")

MPC_EXPECTED_KEYS = [
    "Number", "Principal_desig", "Name", "H", "G", "Epoch",
    "a", "e", "i", "Node", "Peri", "M", "n", "Orbit_type",
    "NEO_flag", "PHA_flag", "Num_obs", "Arc_years", "Arc_length", "Last_obs",
]


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download_mpcorb():
    dest = RAW / f"mpcorb_extended_{TODAY}.json.gz"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[mpcorb] using cached {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    print(f"[mpcorb] downloading {MPCORB_URL}")
    tmp = dest.with_suffix(".part")
    with requests.get(
        MPCORB_URL,
        headers={"Accept-Encoding": "identity"},
        stream=True,
        timeout=(30, 1800),
    ) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[mpcorb] saved {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


def iter_json_records(gz_path):
    """Yield records from a gzipped JSON file that is either a single JSON
    array or JSON lines. The array case is parsed incrementally to keep
    memory bounded."""
    with gzip.open(gz_path, "rt", encoding="utf-8") as f:
        c = f.read(1)
        while c and c.isspace():
            c = f.read(1)
        if not c:
            return
        if c == "[":
            dec = json.JSONDecoder()
            buf = ""
            while True:
                chunk = f.read(1 << 22)
                buf += chunk
                pos = 0
                while True:
                    while pos < len(buf) and buf[pos] in " \t\r\n,":
                        pos += 1
                    if pos < len(buf) and buf[pos] == "]":
                        return
                    try:
                        obj, pos = dec.raw_decode(buf, pos)
                    except json.JSONDecodeError:
                        break
                    yield obj
                buf = buf[pos:]
                if not chunk:
                    if buf.strip(" \t\r\n,]"):
                        raise ValueError(
                            f"unparsed tail of {len(buf)} chars in {gz_path.name}"
                        )
                    return
        else:
            line = c + f.readline()
            while line:
                s = line.strip().rstrip(",")
                if s and s not in ("[", "]"):
                    yield json.loads(s)
                line = f.readline()


def to_int(v):
    if v is None:
        return None
    try:
        return int(str(v).strip("() "))
    except (ValueError, TypeError):
        return None


def to_float(v):
    if v is None or v == "":
        return np.nan
    try:
        return float(v)
    except (ValueError, TypeError):
        return np.nan


def arc_years(rec):
    ay = rec.get("Arc_years")
    if ay:
        try:
            y0, y1 = str(ay).split("-")
            return float(int(y1) - int(y0))
        except (ValueError, TypeError):
            pass
    al = rec.get("Arc_length")
    if al is not None:
        try:
            return float(al) / 365.25
        except (ValueError, TypeError):
            pass
    return np.nan


def parse_mpcorb(gz_path):
    print("[mpcorb] parsing (streaming) ...")
    rows = []
    seen_keys = set()
    first_rec = None
    n_checked = 0
    for rec in iter_json_records(gz_path):
        if n_checked < 1000:
            if first_rec is None:
                first_rec = rec
            seen_keys.update(rec.keys())
            n_checked += 1
            if n_checked == 1000:
                missing = [k for k in MPC_EXPECTED_KEYS if k not in seen_keys]
                if missing:
                    print(f"[mpcorb] WARNING expected fields missing from schema (first 1000 records): {missing}")
                    print(f"[mpcorb] actual top-level keys of first record: {sorted(first_rec.keys())}")
        rows.append((
            to_int(rec.get("Number")),
            rec.get("Principal_desig"),
            rec.get("Name"),
            to_float(rec.get("H")),
            to_float(rec.get("G")),
            to_float(rec.get("Epoch")),
            to_float(rec.get("a")),
            to_float(rec.get("e")),
            to_float(rec.get("i")),
            to_float(rec.get("Node")),
            to_float(rec.get("Peri")),
            to_float(rec.get("M")),
            to_float(rec.get("n")),
            rec.get("Orbit_type"),
            bool(rec.get("NEO_flag", 0)),
            bool(rec.get("PHA_flag", 0)),
            to_int(rec.get("Num_obs")),
            arc_years(rec),
            rec.get("Last_obs"),
        ))
        if len(rows) % 200000 == 0:
            print(f"[mpcorb]   {len(rows):,} records parsed")
    cols = [
        "number_mp", "designation", "mpc_name", "mpc_h", "mpc_g", "mpc_epoch",
        "mpc_a", "mpc_e", "mpc_i", "mpc_node", "mpc_peri", "mpc_m", "mpc_n",
        "mpc_orbit_type", "mpc_neo_flag", "mpc_pha_flag", "mpc_num_obs",
        "mpc_arc_years", "mpc_last_obs",
    ]
    df = pd.DataFrame(rows, columns=cols)
    df["number_mp"] = df["number_mp"].astype("Int64")
    df["mpc_num_obs"] = df["mpc_num_obs"].astype("Int64")
    print(f"[mpcorb] parsed {len(df):,} records")
    return df


def fetch_sbdb():
    """One bulk call for all numbered asteroids (stable NASA bulk API,
    accepted exception to the static-file rule). Cached under data/raw/sbdb/."""
    sbdb_dir = RAW / "sbdb"
    sbdb_dir.mkdir(parents=True, exist_ok=True)
    cache = sbdb_dir / f"sbdb_numbered_{TODAY}.json"
    params = {"fields": SBDB_FIELDS, "sb-ns": "n", "sb-kind": "a"}
    if not (cache.exists() and cache.stat().st_size > 0):
        print("[sbdb] querying JPL SBDB (all numbered asteroids, one bulk call)")
        tmp = cache.with_suffix(".part")
        for attempt in range(1, 4):
            try:
                with requests.get(SBDB_URL, params=params, stream=True,
                                  timeout=(30, 1800)) as r:
                    r.raise_for_status()
                    with open(tmp, "wb") as f:
                        for chunk in r.iter_content(chunk_size=1 << 20):
                            f.write(chunk)
                tmp.rename(cache)
                break
            except (requests.RequestException, OSError) as e:
                print(f"[sbdb] attempt {attempt} failed: {e}")
                if attempt == 3:
                    raise
                time.sleep(15)
        print(f"[sbdb] saved {cache.name} ({cache.stat().st_size/1e6:.1f} MB)")
    else:
        print(f"[sbdb] using cached {cache.name} ({cache.stat().st_size/1e6:.1f} MB)")

    with open(cache) as f:
        payload = json.load(f)
    fields = payload["fields"]
    data = payload["data"]
    df = pd.DataFrame(data, columns=fields)
    del payload, data
    df.columns = ["sbdb_" + c for c in df.columns]

    df["number_mp"] = pd.to_numeric(df["sbdb_pdes"], errors="coerce").astype("Int64")
    n_bad = int(df["number_mp"].isna().sum())
    if n_bad:
        print(f"[sbdb] WARNING {n_bad} rows without numeric pdes dropped")
        df = df[df["number_mp"].notna()]

    for c in ["sbdb_H", "sbdb_diameter", "sbdb_albedo", "sbdb_e", "sbdb_a",
              "sbdb_q", "sbdb_i", "sbdb_om", "sbdb_w", "sbdb_ma", "sbdb_ad",
              "sbdb_per_y", "sbdb_moid", "sbdb_moid_jup", "sbdb_t_jup"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["sbdb_n_obs_used"] = pd.to_numeric(df["sbdb_n_obs_used"], errors="coerce").astype("Int64")
    df["sbdb_neo"] = df["sbdb_neo"].eq("Y")
    df["sbdb_pha"] = df["sbdb_pha"].eq("Y")
    df["sbdb_full_name"] = df["sbdb_full_name"].str.strip()
    print(f"[sbdb] parsed {len(df):,} numbered asteroids")
    return df, params, cache


def write_provenance(name, doc):
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    path = PROVENANCE / f"{name}.json"
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
    print(f"[provenance] wrote {path.relative_to(ROOT)}")


def main():
    for d in (RAW, INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)

    mpcorb_file = download_mpcorb()
    mpcorb_sha = sha256_file(mpcorb_file)
    mpc = parse_mpcorb(mpcorb_file)

    sbdb, sbdb_params, sbdb_cache = fetch_sbdb()

    merged = mpc.merge(sbdb, on="number_mp", how="left")
    n_matched = int(merged["sbdb_spkid"].notna().sum())
    unmatched = merged[merged["number_mp"].notna() & merged["sbdb_spkid"].isna()]
    print(f"[merge] MPCORB rows: {len(mpc):,}  SBDB rows: {len(sbdb):,}  "
          f"matched: {n_matched:,}")
    if len(unmatched):
        print(f"[merge] numbered-but-unmatched: {len(unmatched):,}, examples:")
        for _, r in unmatched.head(5).iterrows():
            print(f"[merge]   number_mp={r['number_mp']} designation={r['designation']}")

    out = INTERIM / "orbits_backbone.parquet"
    merged.to_parquet(out, engine="pyarrow", compression="snappy", index=False)

    retrieved = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    write_provenance("mpcorb", {
        "source": "IAU Minor Planet Center, Extended Orbit File (MPCORB)",
        "url": MPCORB_URL,
        "date": retrieved,
        "file": str(mpcorb_file.relative_to(ROOT)),
        "sha256": mpcorb_sha,
        "size_bytes": mpcorb_file.stat().st_size,
        "row_count": int(len(mpc)),
        "citation": f"IAU Minor Planet Center, MPCORB extended file, retrieved {retrieved}",
        "license_note": ("MPC data policy applies to raw redistribution; raw file is "
                         "cached locally only and the project publishes only derived tables."),
    })
    write_provenance("sbdb", {
        "source": "JPL Small-Body Database Query API",
        "url": SBDB_URL,
        "date": retrieved,
        "params": sbdb_params,
        "file": str(sbdb_cache.relative_to(ROOT)),
        "sha256": sha256_file(sbdb_cache),
        "row_count": int(len(sbdb)),
        "citation": f"JPL Small-Body Database, NASA/JPL-Caltech Solar System Dynamics, retrieved {retrieved}",
        "license_note": "NASA/JPL data, publicly available; cached raw responses kept locally.",
    })

    gasp_keys = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    gasp_n = len(gasp_keys)
    gasp_matched = int(gasp_keys["number_mp"].isin(
        merged.loc[merged["number_mp"].notna(), "number_mp"]).sum())

    n_numbered = int(merged["number_mp"].notna().sum())
    n_neo = int(merged["mpc_neo_flag"].sum())
    n_pha = int(merged["mpc_pha_flag"].sum())
    n_diam = int(merged["sbdb_diameter"].notna().sum())
    n_alb = int(merged["sbdb_albedo"].notna().sum())
    size_mb = out.stat().st_size / 1e6

    print()
    print("=== Coverage report (AP_A1) ===")
    print(f"Total objects:                {len(merged):,}")
    print(f"Numbered:                     {n_numbered:,}")
    print(f"NEOs:                         {n_neo:,}   PHAs: {n_pha:,}")
    print(f"SBDB diameter available:      {n_diam:,} ({100*n_diam/len(merged):.1f}%)")
    print(f"SBDB albedo available:        {n_alb:,} ({100*n_alb/len(merged):.1f}%)")
    print(f"GASP core matched:            {gasp_matched:,} / {gasp_n:,} "
          f"({100*gasp_matched/gasp_n:.2f}%)")
    print(f"Saved: data/interim/orbits_backbone.parquet  {size_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
