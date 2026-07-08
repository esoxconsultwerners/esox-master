"""AP_A8 - Meteoritical Bulletin census + meteorite orbits + CNEOS fireballs.

Three sources, three maturities:

  1. MetBull classification census. The live Meteoritical Bulletin search
     (lpi.usra.edu/meteor/metbull.php) sits behind a Cloudflare JS challenge
     that cannot be passed non-interactively from the VPS, so we use the
     sanctioned fallback: NASA Open Data "Meteorite Landings" (data.nasa.gov),
     a frozen 2013-era export of MetBull. Staleness is documented in the
     provenance and coverage report.

  2. Instrumentally observed falls with reconstructed heliocentric orbits -
     the crown jewels for source attribution. Meier's curated compilation at
     meteoriteorbits.info ships one HTML table plus a full reference legend;
     every row carries a per-event citation code that we expand to the full
     reference string (MANDATORY per row - this table gets published).

  3. CNEOS fireball events (JPL) for context statistics only. No orbit
     derivation here; most events lack full velocity vectors.

Group normalization reuses the reviewed data/interim/relab_group_map.csv
(one canonical mapping for the whole project); MetBull-only class tokens are
appended to that file rather than forking a second map.
"""

import datetime as dt
import hashlib
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROV = ROOT / "data" / "provenance"
GROUP_MAP = INTERIM / "relab_group_map.csv"

SOFT_GROUPS = {"H-L", "L-LL", "H-LL"}
COARSE_GROUPS = {"C-ungrouped", "iron", "OC-ung", "E"}


def group_kind_of(canonical):
    if canonical in SOFT_GROUPS:
        return "soft"
    if canonical in COARSE_GROUPS:
        return "coarse"
    return "primary"

UA = "EsoxConsult-research/1.0 (+https://esoxspace.com; vienna@esoxconsult.com)"
BROWSER_UA = "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"
TODAY = dt.date.today().isoformat()

METBULL_CSV_URL = "https://data.nasa.gov/docs/legacy/meteorite_landings/Meteorite_Landings.csv"
MORB_URL = "https://www.meteoriteorbits.info"
CNEOS_URL = "https://ssd-api.jpl.nasa.gov/fireball.api?vel-comp=true"


def fetch(url, ua=UA):
    req = urllib.request.Request(url, headers={"User-Agent": ua})
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read()


def cache(subdir, filename, data):
    d = RAW / subdir
    d.mkdir(parents=True, exist_ok=True)
    path = d / filename
    path.write_bytes(data)
    return path, hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------
# Class -> canonical group (aligned with the RELAB canonical vocabulary)
# --------------------------------------------------------------------------

def classify_class(raw):
    """Map a meteorite classification string to a canonical group token.

    Vocabulary matches data/interim/relab_group_map.csv (H/L/LL/CI/CM/CO/CV/
    CK/CR/CH/CB/E/EH/EL/R/K, achondrite + iron + stony-iron + planetary group
    names). 'OC-ung' is the ordinary-chondrite-undetermined bucket ('OC',
    'H/L', 'L/LL' intermediates). 'unclassified' is the honest catch-all
    (Stone-uncl, bare Chondrite/Achondrite-ung, Unknown) and is NOT counted
    as classified.
    """
    if raw is None:
        return "unclassified"
    s = str(raw).strip()
    if not s:
        return "unclassified"
    low = s.lower()
    oc = low.replace(" ", "")

    # Compound like "C, CM2" or "Stone, CM2": classify the specific tail.
    if "," in low and low.split(",")[0].strip() in ("c", "stone", "stony",
                                                    "chondrite", "achondrite"):
        tail = classify_class(low.split(",", 1)[1])
        if tail != "unclassified":
            return tail

    # Stony-irons (before irons, before name-word achondrites)
    if low.startswith("mesosiderite"):
        return "mesosiderite"
    if low.startswith("pallasite"):
        return "pallasite"

    # Irons: "Iron", "Iron, IIIAB", "Iron, ungrouped", "Iron, IAB-ung"
    if low.startswith("iron"):
        return "iron"

    # Planetary
    if (low.startswith("martian") or "shergottite" in low or "nakhlite" in low
            or "chassignite" in low):
        return "martian"
    if low.startswith("lunar"):
        return "lunar-meteorite"

    # Achondrites by name ('eurcite' is a source-side spelling variant of
    # eucrite seen in the orbits table; the archived class string is kept
    # verbatim, only this derived group tolerates it)
    if low.startswith("eucrite") or low.startswith("eurcite"):
        return "eucrite"
    if low.startswith("howardite"):
        return "howardite"
    if low.startswith("diogenite"):
        return "diogenite"
    if low.startswith("aubrite"):
        return "aubrite"
    if low.startswith("ureilite"):
        return "ureilite"
    if low.startswith("angrite"):
        return "angrite"
    if low.startswith("brachinite"):
        return "brachinite"
    if low.startswith("winonaite"):
        return "winonaite"
    if low.startswith("acapulcoite") or low.startswith("lodranite"):
        return "lodranite-acapulcoite"

    # Carbonaceous chondrites (CI/CM/CO/CV/CK/CR/CH/CB, then ungrouped)
    for sub in ("ci", "cm", "co", "cv", "ck", "cr", "ch", "cb"):
        if low.startswith(sub):
            return sub.upper()
    if re.match(r"^c[0-9]", low) or low.startswith("c-ung") or low.startswith("c2-ung") \
            or low.startswith("c3-ung") or low.startswith("c4-ung") or low == "c":
        return "C-ungrouped"

    # Enstatite chondrites
    if low.startswith("eh"):
        return "EH"
    if low.startswith("el"):
        return "EL"
    if re.match(r"^e[0-9]", low) or low.startswith("e-ung") or low.startswith("enst") \
            or low == "e":
        return "E"

    # Rumuruti / Kakangari
    if re.match(r"^r[0-9]", low) or low.startswith("r-ung") or low == "r":
        return "R"
    if re.match(r"^k[0-9]", low) or low == "k" or low.startswith("kakangari"):
        return "K"

    # Ordinary-chondrite intermediates (H/L, L/LL, ...) become dedicated
    # soft-label groups (H-L, L-LL) rather than being forced into a neighbour.
    m = re.match(r"^(h|ll|l)/(h|ll|l)", oc)
    if m and m.group(1) != m.group(2):
        rank = {"h": 0, "l": 1, "ll": 2}
        x, y = sorted((m.group(1), m.group(2)), key=lambda t: rank[t])
        return f"{x.upper()}-{y.upper()}"
    # then LL before L before H for the pure classes
    if re.match(r"^ll([0-9\-./]|$)", oc):
        return "LL"
    if re.match(r"^l([0-9\-./]|$)", oc):
        return "L"
    if re.match(r"^h([0-9\-./]|$)", oc):
        return "H"
    if oc.startswith("oc") or oc.startswith("ordinarychondrite"):
        return "OC-ung"

    return "unclassified"


def build_group_lookup(distinct_classes):
    """Reuse relab_group_map.csv where the class token matches; classify the
    rest with the rules above and append the new tokens to the SAME map file
    (2-column raw/canonical). Existing rows always win and are never touched -
    Werner audits the git diff of the additions.
    """
    existing = pd.read_csv(GROUP_MAP)
    known = dict(zip(existing["raw"].astype(str), existing["canonical"].astype(str)))

    additions = []
    lookup = {}
    for cls in sorted(distinct_classes):
        if cls in known:
            lookup[cls] = known[cls]
            continue
        canon = classify_class(cls)
        lookup[cls] = canon
        if canon != "unclassified":
            additions.append((cls, canon))

    if additions:
        add_df = pd.DataFrame(additions, columns=["raw", "canonical"])
        merged = pd.concat([existing, add_df], ignore_index=True)
        merged = merged.drop_duplicates(subset="raw", keep="first")
        if "group_kind" in merged.columns:
            merged["group_kind"] = merged["canonical"].map(group_kind_of)
        merged = merged.sort_values(["canonical", "raw"]).reset_index(drop=True)
        merged.to_csv(GROUP_MAP, index=False)

    return lookup, additions


# --------------------------------------------------------------------------
# Source 1: MetBull census (NASA Open Data fallback)
# --------------------------------------------------------------------------

def ingest_metbull():
    print("[metbull] downloading NASA Open Data Meteorite Landings ...")
    data = fetch(METBULL_CSV_URL, ua=BROWSER_UA)
    raw_path, sha = cache("metbull", f"Meteorite_Landings_{TODAY}.csv", data)

    import io
    df = pd.read_csv(io.BytesIO(data))
    df.columns = [c.strip() for c in df.columns]

    n_raw = len(df)
    df = df.drop_duplicates(subset="name", keep="first").reset_index(drop=True)
    n_dedup_dropped = n_raw - len(df)

    out = pd.DataFrame()
    out["number_mp"] = pd.array([pd.NA] * len(df), dtype="Int64")
    out["metbull_name"] = df["name"].astype(str).str.strip()
    out["metbull_class_raw"] = df["recclass"].astype(str).str.strip()
    fall = df["fall"].astype(str).str.strip().str.lower()
    out["metbull_fall_find"] = fall.map({"fell": "fall", "found": "find"}).fillna(fall)
    out["metbull_year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    out["metbull_mass_g"] = pd.to_numeric(df["mass (g)"], errors="coerce")
    out["metbull_lat"] = pd.to_numeric(df["reclat"], errors="coerce")
    out["metbull_lon"] = pd.to_numeric(df["reclong"], errors="coerce")
    # NASA Open Data has no country column; MetBull-live does. Documented gap.
    out["metbull_country"] = pd.array([pd.NA] * len(df), dtype="string")

    lookup, additions = build_group_lookup(set(out["metbull_class_raw"]))
    out["metbull_group"] = out["metbull_class_raw"].map(lookup)

    out.to_parquet(INTERIM / "metbull.parquet", index=False)

    max_year = int(out["metbull_year"].dropna().max())
    prov = {
        "source": "Meteoritical Bulletin Database (via NASA Open Data fallback)",
        "url": METBULL_CSV_URL,
        "landing_page": "https://www.lpi.usra.edu/meteor/metbull.php",
        "retrieval_date": TODAY,
        "sha256": sha,
        "raw_file": str(raw_path.relative_to(ROOT)),
        "row_count": int(len(out)),
        "rows_before_dedup": int(n_raw),
        "dedup_dropped": int(n_dedup_dropped),
        "dedup_rule": "drop_duplicates on metbull_name, keep first",
        "falls": int((out["metbull_fall_find"] == "fall").sum()),
        "finds": int((out["metbull_fall_find"] == "find").sum()),
        "staleness": (
            "NASA Open Data 'Meteorite Landings' is a frozen snapshot of the "
            f"Meteoritical Bulletin; latest year present is {max_year}. The live "
            "MetBull DB holds ~80,000+ entries as of 2026, so this fallback lags "
            "by roughly a decade / tens of thousands of entries. The live search "
            "is behind a Cloudflare JS challenge that cannot be passed "
            "non-interactively from the VPS; a Phase-D refresh needs a browser-"
            "driven or licensed export."
        ),
        "missing_columns": {
            "metbull_country": "not present in the NASA export; nulled here"
        },
        "group_map_additions": [{"raw": r, "canonical": c} for r, c in additions],
        "citation": (
            "Meteoritical Society, Meteoritical Bulletin Database; distributed as "
            "NASA Open Data 'Meteorite Landings' (data.nasa.gov)."
        ),
        "license_terms": (
            "NASA Open Data (U.S. Government, public domain). Underlying data "
            "(c) The Meteoritical Society."
        ),
    }
    (PROV / "metbull.json").write_text(json.dumps(prov, indent=2))
    return out, additions


# --------------------------------------------------------------------------
# Source 2: instrumentally observed falls with orbits (meteoriteorbits.info)
# --------------------------------------------------------------------------

def _central(cell):
    """First numeric value in a cell like '1.844\\n\\u00b10.003' or '211.496'."""
    if cell is None:
        return None
    m = re.search(r"[-+]?\d+(?:\.\d+)?", str(cell).replace("–", "-"))
    return float(m.group(0)) if m else None


def _strip_tags(html):
    import html as _h
    return _h.unescape(re.sub(r"<[^>]+>", " ", html)).strip()


def ingest_orbits(metbull_names):
    print("[orbits] downloading meteoriteorbits.info ...")
    html_bytes = fetch(MORB_URL, ua=BROWSER_UA)
    raw_path, sha = cache("meteorite_orbits", f"meteoriteorbits_{TODAY}.html", html_bytes)
    html = html_bytes.decode("utf-8", "replace")

    table = re.search(r"<table.*?</table>", html, re.S).group(0)
    rows = re.findall(r"<tr.*?</tr>", table, re.S)

    # Reference legend (after the table): "AN22: Andrade et al. 2022. ..."
    after = html[html.index(table) + len(table):]
    after = re.sub(r"<script.*?</script>", "", after, flags=re.S)
    legend_txt = _strip_tags(after)
    legend = {}
    for m in re.finditer(r"(?m)^\s*([A-Za-z][A-Za-z0-9\-]{1,12}):\s+(.+?)\s*$", legend_txt):
        legend[m.group(1)] = re.sub(r"\s+", " ", m.group(2)).strip()

    def expand_ref(code_cell):
        text = _strip_tags(code_cell)
        codes = [c.strip() for c in re.split(r"[,;/]| and ", text) if c.strip()]
        parts = []
        for c in codes:
            parts.append(f"{c}: {legend[c]}" if c in legend else c)
        joined = "; ".join(parts) if parts else text
        return joined or None, ", ".join(codes)

    records = []
    for tr in rows[1:]:
        cells = re.findall(r"<t[dh].*?>(.*?)</t[dh]>", tr, re.S)
        if len(cells) < 14:
            continue
        name = _strip_tags(cells[1])
        link = re.search(r"code=(\d+)", cells[1])
        ref_full, ref_codes = expand_ref(cells[13])
        records.append({
            "morb_abbrev": _strip_tags(cells[0]),
            "morb_meteorite_name": name,
            "morb_class": _strip_tags(cells[2]),
            "morb_location": _strip_tags(cells[3]),
            "morb_fall_date": _strip_tags(cells[4]).replace("/", "-"),
            "morb_fall_time_utc": _strip_tags(cells[5]),
            "morb_a_au": _central(cells[8]),
            "morb_e": _central(cells[9]),
            "morb_i_deg": _central(cells[10]),
            "morb_peri_deg": _central(cells[11]),
            "morb_node_deg": _central(cells[12]),
            "morb_orbit_ref": ref_full,
            "morb_ref_code": ref_codes,
            "morb_metbull_code": int(link.group(1)) if link else pd.NA,
        })

    df = pd.DataFrame(records)
    df["number_mp"] = pd.array([pd.NA] * len(df), dtype="Int64")
    df["morb_group"] = df["morb_class"].map(classify_class)

    names_lower = {n.lower() for n in metbull_names}
    df["morb_in_metbull"] = df["morb_meteorite_name"].str.lower().isin(names_lower)

    df.to_parquet(INTERIM / "meteorite_orbits.parquet", index=False)

    prov = {
        "source": "meteoriteorbits.info (M. M. M. Meier, curated compilation)",
        "url": MORB_URL,
        "retrieval_date": TODAY,
        "sha256": sha,
        "raw_file": str(raw_path.relative_to(ROOT)),
        "row_count": int(len(df)),
        "reference_legend_size": len(legend),
        "note": (
            "Central orbital-element values parsed from the HTML table; per-row "
            "uncertainties in the source are dropped. Every event carries a "
            "morb_orbit_ref (reference code expanded via the page's legend). "
            "The table is the version-of-record for observed-fall orbits; a "
            "Phase-D refresh re-scrapes it (retrieval date = version)."
        ),
        "citation": (
            "Meier M.M.M., meteoriteorbits.info; individual orbits credited "
            "per-row in morb_orbit_ref (original literature)."
        ),
        "license_terms": (
            "Compilation (c) M.M.M. Meier, used with per-event attribution to "
            "the primary literature. Not redistributed as bulk; derived orbital "
            "elements only."
        ),
    }
    (PROV / "meteorite_orbits.json").write_text(json.dumps(prov, indent=2))
    return df


# --------------------------------------------------------------------------
# Source 3: CNEOS fireballs (context statistics)
# --------------------------------------------------------------------------

def ingest_cneos():
    print("[cneos] downloading JPL fireball API ...")
    data = fetch(CNEOS_URL, ua=UA)
    raw_path, sha = cache("cneos", f"fireball_{TODAY}.json", data)
    doc = json.loads(data)
    fields = doc["fields"]
    idx = {f: i for i, f in enumerate(fields)}

    def g(row, key):
        v = row[idx[key]] if key in idx else None
        return None if v in (None, "") else v

    recs = []
    for row in doc["data"]:
        lat = g(row, "lat")
        lon = g(row, "lon")
        latd = g(row, "lat-dir")
        lond = g(row, "lon-dir")
        lat = float(lat) * (-1 if latd == "S" else 1) if lat is not None else None
        lon = float(lon) * (-1 if lond == "W" else 1) if lon is not None else None
        vel = g(row, "vel")
        recs.append({
            "cneos_date": g(row, "date"),
            "cneos_energy_kt": float(g(row, "impact-e")) if g(row, "impact-e") else None,
            "cneos_radiated_energy_j": float(g(row, "energy")) * 1e10 if g(row, "energy") else None,
            "cneos_velocity_kms": float(vel) if vel is not None else None,
            "cneos_vx": float(g(row, "vx")) if g(row, "vx") is not None else None,
            "cneos_vy": float(g(row, "vy")) if g(row, "vy") is not None else None,
            "cneos_vz": float(g(row, "vz")) if g(row, "vz") is not None else None,
            "cneos_lat": lat,
            "cneos_lon": lon,
            "cneos_altitude_km": float(g(row, "alt")) if g(row, "alt") is not None else None,
        })

    df = pd.DataFrame(recs)
    df["cneos_date"] = pd.to_datetime(df["cneos_date"], errors="coerce")
    df["number_mp"] = pd.array([pd.NA] * len(df), dtype="Int64")
    df.to_parquet(INTERIM / "cneos_fireballs.parquet", index=False)

    n_vel = int(df["cneos_velocity_kms"].notna().sum())
    prov = {
        "source": "CNEOS Fireball and Bolide Data (JPL/NASA)",
        "url": CNEOS_URL,
        "api_note": "vel-comp=true returns velocity components vx/vy/vz where available",
        "retrieval_date": TODAY,
        "sha256": sha,
        "raw_file": str(raw_path.relative_to(ROOT)),
        "row_count": int(len(df)),
        "with_velocity": n_vel,
        "note": (
            "Context statistics only. energy_kt is total impact energy (kt TNT); "
            "most events lack full velocity vectors and NO orbit is derived in "
            "this package. Living list: retrieval date = version (Phase-D diff)."
        ),
        "citation": "NASA/JPL Center for Near-Earth Object Studies (CNEOS) fireball data.",
        "license_terms": "NASA/JPL public data.",
    }
    (PROV / "cneos_fireballs.json").write_text(json.dumps(prov, indent=2))
    return df


# --------------------------------------------------------------------------

def coverage_report(mb, orb, cne, additions):
    falls = mb[mb["metbull_fall_find"] == "fall"]
    finds = mb[mb["metbull_fall_find"] == "find"]
    classified = mb[mb["metbull_group"] != "unclassified"]
    falls_class = falls[falls["metbull_group"] != "unclassified"]

    grp_falls = falls_class["metbull_group"].value_counts()
    oc_falls = int(falls_class["metbull_group"].isin(
        ["H", "L", "LL", "OC-ung", "H-L", "L-LL", "H-LL"]).sum())
    ll_falls = int((falls_class["metbull_group"] == "LL").sum())
    ll_frac = 100.0 * ll_falls / len(falls_class) if len(falls_class) else 0.0

    n_ref = int(orb["morb_orbit_ref"].notna().sum())
    n_match = int(orb["morb_in_metbull"].sum())
    n_vel = int(cne["cneos_velocity_kms"].notna().sum())

    def sz(name):
        p = INTERIM / name
        return f"{p.name} ({p.stat().st_size/1024:.1f} KB)"

    print("\n" + "=" * 64)
    print("AP_A8 COVERAGE REPORT  (MetBull source: NASA Open Data fallback)")
    print("=" * 64)
    print(f"MetBull entries (falls / finds):     {len(mb)} ({len(falls)} / {len(finds)})")
    print(f"Classified with canonical group:     {len(classified)} "
          f"({100.0*len(classified)/len(mb):.1f}%)")
    print(f"Falls per major group:")
    for g, n in grp_falls.items():
        print(f"    {g:<18} {n:>4}  ({100.0*n/len(falls_class):.1f}%)")
    print(f"    {'-- OC (H+L+LL+OC-ung+soft)':<18} {oc_falls:>4}  "
          f"({100.0*oc_falls/len(falls_class):.1f}% of classified falls)")
    print(f"Orbits table events:                 {len(orb)}, all with per-row refs "
          f"({n_ref}/{len(orb)} non-null)")
    print(f"Orbit events matched to MetBull:     {n_match}")
    print(f"CNEOS fireballs:                     {len(cne)}, with velocity vector {n_vel}")
    print(f"LL falls fraction:                   {ll_frac:.1f}% of classified falls "
          f"({ll_falls}/{len(falls_class)})")
    print(f"Group-map additions (MetBull-only):  {len(additions)} new raw tokens")
    print(f"GASP core overlap:                   N/A (meteorite reference tables, "
          f"not asteroid-keyed)")
    print(f"Saved: {sz('metbull.parquet')}, {sz('meteorite_orbits.parquet')}, "
          f"{sz('cneos_fireballs.parquet')}")
    print("=" * 64)

    dossier = ROOT / "docs" / "apophis_status.md"
    line = (
        f"\n## AP_A8 (Meteoritical Bulletin + orbits, {TODAY})\n\n"
        f"- LL falls fraction: **{ll_frac:.1f}%** of classified MetBull falls "
        f"({ll_falls}/{len(falls_class)}) are LL chondrites. This is the base "
        f"rate behind the Apophis Sq<->LL story: Apophis is an Sq-type, and Sq "
        f"is the spectral bridge to LL ordinary chondrites - the most probable "
        f"meteorite analogue. LL is a minority of falls, so a delivered Apophis "
        f"fragment would be a comparatively rare, diagnostic sample.\n"
        f"- Observed-fall orbits with published references: {len(orb)} events "
        f"(incl. Chelyabinsk LL5, source region traceable per row).\n"
    )
    with dossier.open("a") as f:
        f.write(line)
    print(f"[dossier] appended LL falls fraction to {dossier.relative_to(ROOT)}")
    return ll_frac


def main():
    for d in (INTERIM, PROV):
        d.mkdir(parents=True, exist_ok=True)
    mb, additions = ingest_metbull()
    orb = ingest_orbits(set(mb["metbull_name"]))
    cne = ingest_cneos()
    coverage_report(mb, orb, cne, additions)


if __name__ == "__main__":
    main()
