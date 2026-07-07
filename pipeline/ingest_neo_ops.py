#!/usr/bin/env python3
"""AP_A6: NHATS accessibility + NEOCC/Sentry risk data.

Three living lists; the retrieval datetime IS the version (documented
NASA/ESA endpoints are the citable primary source here - accepted
exception to the static-file rule, as with SBDB in AP_A1).

  nhats_*:  JPL NHATS accessibility table (summary mode, API defaults:
            dv <= 12 km/s, duration <= 450 d, stay >= 8 d, launch
            window 2020-2045, no H / OCC filter). The API provides the
            next OBSERVATION window (obs_start/obs_end), not a launch
            window - kept as nhats_obs_start / nhats_obs_end.
  neocc_*:  ESA NEOCC risk list (pipe-delimited text download).
  sentry_*: JPL Sentry summary (cross-reference only - we ingest
            monitoring results, never compute impact probabilities).

Key resolution (important for unnumbered objects / later LSST work):
the canonical join key is str(number_mp) for numbered objects, else
the provisional designation normalized to MPC spacing ("2023VD3" ->
"2023 VD3"). NEOCC prints compact designations; JPL prints spaced
ones - both normalize to the same key.
"""

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
DOCS = ROOT / "docs"

NHATS_URL = "https://ssd-api.jpl.nasa.gov/nhats.api"
SENTRY_URL = "https://ssd-api.jpl.nasa.gov/sentry.api"
NEOCC_URL = "https://neo.ssa.esa.int/PSDB-portlet/download?file=esa_risk_list"

TODAY = datetime.now(timezone.utc).strftime("%Y%m%d")
NOW_ISO = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download(url, dest, ua="esox-master-catalog/A6"):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[neo-ops] using cached {dest.name}")
        return dest
    print(f"[neo-ops] downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, headers={"User-Agent": ua}, stream=True,
                      timeout=(30, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    tmp.rename(dest)
    print(f"[neo-ops] saved {dest.name} ({dest.stat().st_size/1e3:.0f} kB)")
    return dest


def to_num(v):
    if v is None or v == "":
        return np.nan
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def normalize_desig(s):
    """'2023VD3' -> '2023 VD3'; spaced forms pass through."""
    s = (s or "").strip().upper()
    m = re.match(r"^(\d{4})\s*([A-Z]{1,2}\d*)$", s)
    if m:
        return f"{m.group(1)} {m.group(2)}"
    return s or None


def resolve_key(des):
    """-> (key, number_mp, designation)."""
    des = (des or "").strip()
    if des.isdigit():
        return des, int(des), None
    d = normalize_desig(des)
    return d, None, d


def ingest_nhats():
    raw = download(NHATS_URL, RAW / "nhats" / f"nhats_{TODAY}.json")
    payload = json.loads(raw.read_text())
    rows = []
    for r in payload["data"]:
        key, num, desig = resolve_key(r["des"])
        rows.append({
            "key": key, "number_mp": num, "designation": desig,
            "nhats_min_dv_kms": to_num(r["min_dv"]["dv"]),
            "nhats_min_dv_duration_d": to_num(r["min_dv"]["dur"]),
            "nhats_min_dur_d": to_num(r["min_dur"]["dur"]),
            "nhats_n_trajectories": int(r["n_via_traj"]),
            "nhats_h": to_num(r.get("h")),
            "nhats_size_min_m": to_num(r.get("min_size")),
            "nhats_size_max_m": to_num(r.get("max_size")),
            "nhats_obs_start": r.get("obs_start"),
            "nhats_obs_end": r.get("obs_end"),
        })
    df = pd.DataFrame(rows)
    print(f"[nhats] {len(df):,} objects (API count: {payload['count']})")
    return df, raw, payload["signature"]


def ingest_sentry():
    raw = download(SENTRY_URL, RAW / "sentry" / f"sentry_{TODAY}.json")
    payload = json.loads(raw.read_text())
    rows = []
    for r in payload["data"]:
        key, num, desig = resolve_key(r["des"])
        rows.append({
            "key": key, "number_mp": num, "designation": desig,
            "sentry_ps_max": to_num(r.get("ps_max")),
            "sentry_ps_cum": to_num(r.get("ps_cum")),
            "sentry_ts_max": to_num(r.get("ts_max")),
            "sentry_ip": to_num(r.get("ip")),
            "sentry_n_imp": int(r["n_imp"]) if r.get("n_imp") is not None else None,
        })
    df = pd.DataFrame(rows)
    print(f"[sentry] {len(df):,} objects (API count: {payload['count']})")
    return df, raw, payload["signature"]


def ingest_neocc():
    raw = download(NEOCC_URL, RAW / "neocc" / f"esa_risk_list_{TODAY}.txt",
                   ua="Mozilla/5.0 (esox-master-catalog/A6)")
    text = raw.read_text(encoding="utf-8", errors="replace")
    last_update = None
    rows = []
    for line in text.splitlines():
        if line.startswith("Last Update"):
            last_update = line.split(":", 1)[1].strip()
            continue
        if "|" not in line or "Num/des" in line or line.lstrip().startswith(("Object", "AAAA")):
            continue
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 11 or not parts[0]:
            continue
        toks = parts[0].split()
        des = toks[0]
        key, num, desig = resolve_key(des)
        if num is not None and len(toks) == 1:
            pass
        rows.append({
            "key": key, "number_mp": num, "designation": desig,
            "neocc_name": " ".join(toks[1:]) or None,
            "neocc_diameter_m": to_num(parts[1]),
            "neocc_vi_max_datetime": parts[3] or None,
            "neocc_ip_max": to_num(parts[4]),
            "neocc_risk_ps_max": to_num(parts[5]),
            "neocc_risk_ts": to_num(parts[6]),
            "neocc_velocity_kms": to_num(parts[7]),
            "neocc_vi_years": parts[8] or None,
            "neocc_ip_cum": to_num(parts[9]),
            "neocc_risk_ps_cum": to_num(parts[10]),
        })
    df = pd.DataFrame(rows)
    print(f"[neocc] {len(df):,} objects (list last update: {last_update})")
    return df, raw, last_update


def dedup(df, name):
    n = int(df["key"].duplicated().sum())
    if n:
        print(f"[{name}] dropping {n} duplicate keys (keeping first)")
        df = df.drop_duplicates("key")
    return df


def main():
    for d in (INTERIM, PROVENANCE):
        d.mkdir(parents=True, exist_ok=True)
    for sub in ("nhats", "sentry", "neocc"):
        (RAW / sub).mkdir(parents=True, exist_ok=True)

    nhats, nhats_raw, nhats_sig = ingest_nhats()
    sentry, sentry_raw, sentry_sig = ingest_sentry()
    neocc, neocc_raw, neocc_update = ingest_neocc()

    nhats = dedup(nhats, "nhats")
    sentry = dedup(sentry, "sentry")
    neocc = dedup(neocc, "neocc")

    merged = nhats.merge(neocc, on="key", how="outer", suffixes=("", "_neocc"))
    merged = merged.merge(sentry, on="key", how="outer", suffixes=("", "_sentry"))
    for col in ("number_mp", "designation"):
        for suf in ("_neocc", "_sentry"):
            if col + suf in merged:
                merged[col] = merged[col].combine_first(merged[col + suf])
                merged = merged.drop(columns=[col + suf])
    merged["number_mp"] = merged["number_mp"].astype("Int64")
    merged = merged.drop(columns=["key"])

    out = INTERIM / "neo_ops.parquet"
    merged.to_parquet(out, engine="pyarrow", compression="snappy", index=False)

    prov_common = {
        "date": NOW_ISO[:10],
        "retrieved_utc": NOW_ISO,
        "note": "living list - the retrieval datetime is the version",
    }
    PROVENANCE.joinpath("nhats.json").write_text(json.dumps({
        "source": "NASA/JPL NHATS Data API (summary mode)",
        "url": NHATS_URL,
        "api_signature": nhats_sig,
        "default_filters": ("dv<=12 km/s, duration<=450 d, stay>=8 d, "
                            "launch window 2020-2045, no H filter, no OCC filter"),
        "next_window_note": ("API provides next observation window "
                             "(obs_start/obs_end), not a launch window; kept as "
                             "nhats_obs_start/nhats_obs_end"),
        **prov_common,
        "sha256": sha256_file(nhats_raw),
        "row_count": int(len(nhats)),
        "citation": ("NASA/JPL Near-Earth Object Human Space Flight Accessible "
                     f"Targets Study (NHATS), {NHATS_URL}, retrieved {NOW_ISO}"),
    }, indent=2))
    PROVENANCE.joinpath("neocc_risk.json").write_text(json.dumps({
        "source": "ESA NEO Coordination Centre Risk List",
        "url": NEOCC_URL,
        "list_last_update": neocc_update,
        **prov_common,
        "sha256": sha256_file(neocc_raw),
        "row_count": int(len(neocc)),
        "citation": ("ESA NEO Coordination Centre risk list, https://neo.ssa.esa.int, "
                     f"list update {neocc_update}, retrieved {NOW_ISO}"),
    }, indent=2))
    PROVENANCE.joinpath("sentry.json").write_text(json.dumps({
        "source": "NASA/JPL Sentry Impact Monitoring API (summary mode)",
        "url": SENTRY_URL,
        "api_signature": sentry_sig,
        "scope_note": ("cross-reference only: monitoring results ingested as "
                       "columns, impact probabilities never computed by us"),
        **prov_common,
        "sha256": sha256_file(sentry_raw),
        "row_count": int(len(sentry)),
        "citation": ("NASA/JPL Center for NEO Studies Sentry impact monitoring, "
                     f"{SENTRY_URL}, retrieved {NOW_ISO}"),
    }, indent=2))
    print("[provenance] wrote nhats.json, neocc_risk.json, sentry.json")

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    our_numbers = set(merged.loc[merged["number_mp"].notna(), "number_mp"])
    g_any = int(gasp["number_mp"].isin(our_numbers).sum())

    backbone = pd.read_parquet(INTERIM / "orbits_backbone.parquet",
                               columns=["number_mp", "designation", "mpc_neo_flag"])
    neos = backbone[backbone["mpc_neo_flag"]]
    our_desigs = set(merged.loc[merged["designation"].notna(), "designation"])
    bb_hit = (neos["number_mp"].isin(our_numbers)
              | neos["designation"].isin(our_desigs))
    n_backbone = int(bb_hit.sum())

    overlap = len(set(neocc["key"]) & set(sentry["key"]))
    n_numbered = int(merged["number_mp"].notna().sum())

    apo_nhats = nhats[nhats["number_mp"] == 99942]
    apo_neocc = (neocc["number_mp"] == 99942).any() or \
        neocc["designation"].eq("2004 MN4").any()
    apo_sentry = (sentry["number_mp"] == 99942).any() or \
        sentry["designation"].eq("2004 MN4").any()
    risk_line = ("ABSENT from NEOCC and Sentry as expected (removed 2021 after "
                 "radar ruled out impacts for 100+ years)"
                 if not (apo_neocc or apo_sentry)
                 else f"UNEXPECTED presence: NEOCC={apo_neocc}, Sentry={apo_sentry}")
    if len(apo_nhats):
        r = apo_nhats.iloc[0]
        nhats_line = (f"PRESENT in NHATS: min dv {r['nhats_min_dv_kms']} km/s "
                      f"({r['nhats_min_dv_duration_d']:.0f} d), "
                      f"{r['nhats_n_trajectories']:,} viable trajectories")
    else:
        nhats_line = "absent from NHATS"

    print()
    print("=== Coverage report (AP_A6) ===")
    print(f"NHATS objects:                 {len(nhats):,}")
    print(f"NEOCC risk list objects:       {len(neocc):,}")
    print(f"Sentry objects:                {len(sentry):,}")
    print(f"NEOCC x Sentry overlap:        {overlap:,}  (sanity: should be large)")
    print(f"Numbered / unnumbered total:   {n_numbered:,} / {len(merged) - n_numbered:,}")
    print(f"GASP core matched (any):       {g_any:,} / {len(gasp):,} "
          f"({100*g_any/len(gasp):.2f}%)")
    print(f"Backbone matched (any):        {n_backbone:,} of {len(neos):,} orbit backbone NEOs")
    print(f"Apophis (99942) check:")
    print(f"    risk lists: {risk_line}")
    print(f"    NHATS: {nhats_line}")
    print(f"Saved: data/interim/neo_ops.parquet  {out.stat().st_size/1e6:.2f} MB")

    status = DOCS / "apophis_status.md"
    if status.exists() and "| A6 |" not in status.read_text():
        lines = status.read_text().splitlines()
        last_row = max(i for i, l in enumerate(lines) if l.startswith("|"))
        lines.insert(last_row + 1,
                     f"| A6 | NEOCC + Sentry risk lists | Absent from both, as expected — "
                     f"removed 2021 after the radar campaign ruled out impacts for 100+ years "
                     f"(2004: Torino 4, highest ever measured) |")
        lines.insert(last_row + 2,
                     f"| A6 | JPL NHATS | {nhats_line}; radar windows 2028-09 (Arecibo-class) "
                     f"/ 2029-04 (Goldstone) archived in raw response |")
        status.write_text("\n".join(lines) + "\n")
        print("[docs] appended Apophis lines to docs/apophis_status.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
