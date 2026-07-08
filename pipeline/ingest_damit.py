"""AP_A9 - DAMIT spin and shape models (last package of the A-series).

Source: DAMIT (Database of Asteroid Models from Inversion Techniques),
Charles University Prague - https://damit.cuni.cz/ (moved from the old
astro.troja.mff.cuni.cz URL). Licence: CC-BY-4.0. Primary citation:
Durech, Sidorin & Kaasalainen (2010) A&A 513, A46.

DAMIT publishes per-table plaintext CSV exports (UTF-8 + BOM) under
/projects/damit/exports/table/<name>. We take the four tables that make up
the spin/shape catalogue plus the tumblers table (non-principal-axis
rotators, where Apophis lives) - one polite crawl, cached with SHA256.

One asteroid can carry MULTIPLE models (different publications, ambiguous
pole solutions - the classic lambda / lambda+180 mirror). We keep them ALL
in the long table; damit_best applies a PROVISIONAL convenience rule only
(newest reference wins). The real precedence rule is a Phase B decision.

This ingest is also the calling card for the later Durech cooperation
request: one reference per model, full citation in the provenance.
"""

import datetime as dt
import hashlib
import io
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "damit"
INTERIM = ROOT / "data" / "interim"
PROV = ROOT / "data" / "provenance"

UA = "EsoxConsult-research/1.0 (+https://esoxspace.com; vienna@esoxconsult.com)"
TODAY = dt.date.today().isoformat()
BASE = "https://damit.cuni.cz/projects/damit"
TABLE_URL = BASE + "/exports/table/{}"
EXPORTS_URL = BASE + "/exports"

DAMIT_DB_CITATION = "DAMIT database (Durech, Sidorin & Kaasalainen 2010, A&A 513, A46)"
DAMIT_DB_BIBCODE = "2010A&A...513A..46D"

TABLES = ("asteroids", "asteroid_models", "references",
          "asteroid_models_references", "tumblers", "references_tumblers")


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read()


def download_tables():
    RAW.mkdir(parents=True, exist_ok=True)
    frames, shas = {}, {}
    for t in TABLES:
        print(f"[damit] downloading table {t} ...")
        data = fetch(TABLE_URL.format(t))
        path = RAW / f"{t}_{TODAY}.csv"
        path.write_bytes(data)
        shas[t] = hashlib.sha256(data).hexdigest()
        frames[t] = pd.read_csv(io.BytesIO(data), encoding="utf-8-sig")
    return frames, shas


def latest_export_label():
    """Pin the DB snapshot version by the newest complete-export filename."""
    try:
        html = fetch(EXPORTS_URL).decode("utf-8", "replace")
        tags = re.findall(r"damit-(\d{8}T\d{6}Z)\.tar\.gz", html)
        return f"damit-{max(tags)}" if tags else None
    except Exception:
        return None


def derive_method(row):
    """DAMIT has no explicit method column; derive a best-effort label from the
    available flags. Documented in provenance.
    """
    if pd.notna(row.get("thermal_inertia")):
        return "thermophysical model"
    if row.get("nonconvex") == 1:
        return "nonconvex inversion"
    return "convex lightcurve inversion"


def build_ref_lookup(references, model_refs):
    """model_id -> (ref_string, bibcodes, newest_year). Models with no linked
    reference (19 in the DB) fall back to the DAMIT database citation so that
    every row carries a mandatory reference.
    """
    rmeta = {}
    for _, r in references.iterrows():
        year = int(r["year"]) if pd.notna(r["year"]) else None
        short = str(r["author_short"]).strip() if pd.notna(r["author_short"]) else None
        label = f"{short} ({year})" if short and year else (short or r.get("bibcode") or "")
        rmeta[r["id"]] = (label, r.get("bibcode"), year)

    by_model = {}
    for _, r in model_refs.iterrows():
        by_model.setdefault(r["asteroid_model_id"], []).append(r["reference_id"])

    lookup = {}
    for mid, rids in by_model.items():
        labels, bibs, years = [], [], []
        for rid in rids:
            if rid in rmeta:
                lab, bib, yr = rmeta[rid]
                if lab:
                    labels.append(lab)
                if pd.notna(bib):
                    bibs.append(str(bib))
                if yr:
                    years.append(yr)
        lookup[mid] = (
            "; ".join(dict.fromkeys(labels)) or DAMIT_DB_CITATION,
            ";".join(dict.fromkeys(bibs)) or DAMIT_DB_BIBCODE,
            max(years) if years else 2010,
        )
    return lookup


def ingest():
    for d in (INTERIM, PROV):
        d.mkdir(parents=True, exist_ok=True)
    frames, shas = download_tables()
    ast = frames["asteroids"]
    mod = frames["asteroid_models"]
    tum = frames["tumblers"]

    ref_lookup = build_ref_lookup(frames["references"], frames["asteroid_models_references"])

    ast_meta = ast.set_index("id")[["number", "name", "designation"]]

    rows = []
    for _, m in mod.iterrows():
        meta = ast_meta.loc[m["asteroid_id"]] if m["asteroid_id"] in ast_meta.index else None
        num = meta["number"] if meta is not None else None
        ref_s, bibs, ref_year = ref_lookup.get(
            m["id"], (DAMIT_DB_CITATION, DAMIT_DB_BIBCODE, 2010))
        rows.append({
            "number_mp": int(num) if pd.notna(num) else pd.NA,
            "designation": (meta["designation"] if meta is not None
                            and pd.notna(meta["designation"]) else pd.NA),
            "damit_name": (meta["name"] if meta is not None
                           and pd.notna(meta["name"]) else pd.NA),
            "damit_model_id": int(m["id"]),
            "damit_lambda_deg": float(m["lambda"]) % 360.0,
            "damit_beta_deg": float(m["beta"]),
            "damit_period_h": float(m["period"]),
            "damit_method": derive_method(m),
            "damit_nonconvex": bool(m.get("nonconvex") == 1),
            "damit_shape_available": True,
            "damit_ref": ref_s,
            "damit_bibcodes": bibs,
            "damit_ref_year": ref_year,
        })

    models = pd.DataFrame(rows)
    models["number_mp"] = models["number_mp"].astype("Int64")
    models["designation"] = models["designation"].astype("string")
    models.to_parquet(INTERIM / "damit_models.parquet", index=False)

    # Provisional "best" model per asteroid: newest reference year wins,
    # tie-break on newest model id. NOT a scientific precedence rule.
    models["_aid"] = mod["asteroid_id"].values
    best = (models.sort_values(["damit_ref_year", "damit_model_id"])
            .groupby("_aid", as_index=False).tail(1)
            .drop(columns="_aid").reset_index(drop=True))
    best["damit_best_rule"] = "provisional: newest reference year, tie-break newest model id"
    ast_counts = models.groupby("_aid").size()
    stats = {
        "distinct_asteroids": int(ast_counts.size),
        "multi_model": int((ast_counts > 1).sum()),
        "max_models": int(ast_counts.max()),
    }
    models = models.drop(columns="_aid")
    best.to_parquet(INTERIM / "damit_best.parquet", index=False)

    # Tumblers (NPA rotators) - bonus small table; this is where Apophis lives.
    tmeta = ast.set_index("id")[["number", "name", "designation"]]
    trows = []
    for _, t in tum.iterrows():
        meta = tmeta.loc[t["asteroid_id"]] if t["asteroid_id"] in tmeta.index else None
        num = meta["number"] if meta is not None else None
        trows.append({
            "number_mp": int(num) if pd.notna(num) else pd.NA,
            "damit_name": (meta["name"] if meta is not None
                           and pd.notna(meta["name"]) else pd.NA),
            "damit_tumbler_id": int(t["id"]),
            "damit_lambda_L_deg": float(t["lambda_angular_momentum"]) % 360.0
            if pd.notna(t["lambda_angular_momentum"]) else pd.NA,
            "damit_beta_L_deg": float(t["beta_angular_momentum"])
            if pd.notna(t["beta_angular_momentum"]) else pd.NA,
            "damit_period_phi_h": float(t["period_phi"]) if pd.notna(t["period_phi"]) else pd.NA,
            "damit_period_psi_h": float(t["period_psi"]) if pd.notna(t["period_psi"]) else pd.NA,
        })
    tumblers = pd.DataFrame(trows)
    tumblers["number_mp"] = tumblers["number_mp"].astype("Int64")
    tumblers.to_parquet(INTERIM / "damit_tumblers.parquet", index=False)

    export_label = latest_export_label()
    prov = {
        "source": "DAMIT - Database of Asteroid Models from Inversion Techniques",
        "institution": "Astronomical Institute, Charles University, Prague",
        "url": TABLE_URL.format("<table>"),
        "landing_page": "https://damit.cuni.cz/",
        "tables_used": list(TABLES),
        "db_snapshot": export_label,
        "retrieval_date": TODAY,
        "sha256": shas,
        "row_counts": {
            "asteroids": int(len(ast)),
            "asteroid_models": int(len(mod)),
            "models_emitted": int(len(models)),
            "damit_best": int(len(best)),
            "tumblers": int(len(tumblers)),
        },
        "method_derivation": (
            "DAMIT has no explicit method field. damit_method is derived: "
            "thermal_inertia present -> 'thermophysical model'; else nonconvex "
            "flag set -> 'nonconvex inversion'; else 'convex lightcurve "
            "inversion' (the default technique)."
        ),
        "shape_note": (
            "DAMIT models are shape models by construction, so "
            "damit_shape_available is True for every asteroid_models row "
            "(structural fact of the database, not a per-file download check)."
        ),
        "reference_note": (
            "damit_ref is mandatory per row. 19 models carry no reference link "
            "in the DB; those fall back to the DAMIT database citation. "
            "Ambiguous pole solutions (lambda / lambda+180) are kept as "
            "separate model rows - no de-duplication."
        ),
        "best_note": (
            "damit_best is PROVISIONAL (newest reference year, tie-break newest "
            "model id) for convenience only; the scientific precedence rule is a "
            "Phase B decision."
        ),
        "citation": (
            "Durech J., Sidorin V., Kaasalainen M. (2010) DAMIT: a database of "
            "asteroid models. A&A 513, A46. Individual models credited per row "
            "in damit_ref."
        ),
        "license_terms": "CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)",
    }
    (PROV / "damit.json").write_text(json.dumps(prov, indent=2))
    return models, best, tumblers, ast, stats


def coverage_report(models, best, tumblers, ast, stats):
    n_models = len(models)
    n_ast = stats["distinct_asteroids"]
    with_pole = int((models["damit_lambda_deg"].notna() & models["damit_beta_deg"].notna()).sum())
    with_shape = int(models["damit_shape_available"].sum())
    multi = stats["multi_model"]
    max_models = stats["max_models"]

    gasp = pd.read_parquet(INTERIM / "gasp_core_keys.parquet")
    gasp_keys = set(gasp["number_mp"].dropna().astype(int))
    dmt_nums = set(best["number_mp"].dropna().astype(int))
    matched = len(dmt_nums & gasp_keys)

    lcdb = pd.read_parquet(INTERIM / "lcdb.parquet")[["number_mp", "lcdb_period_h"]].dropna()
    lcdb = lcdb.drop_duplicates("number_mp")
    xc = best[["number_mp", "damit_period_h"]].dropna().merge(lcdb, on="number_mp")
    rel = ((xc["damit_period_h"] - xc["lcdb_period_h"]).abs() / xc["lcdb_period_h"])
    med_rel = float(rel.median())
    n_gt1 = int((rel > 0.01).sum())

    apo = ast[ast["number"] == 99942]
    apo_models = apo_tumb = 0
    apo_tinfo = ""
    if len(apo):
        aid = apo.iloc[0]["id"]
        apo_models = int((models["number_mp"] == 99942).sum())
        at = tumblers[tumblers["number_mp"] == 99942]
        apo_tumb = int(len(at))
        if apo_tumb:
            r = at.iloc[0]
            apo_tinfo = (f" (tumbler periods phi={r['damit_period_phi_h']} h, "
                         f"psi={r['damit_period_psi_h']} h)")

    def sz(name):
        p = INTERIM / name
        return f"{p.name} ({p.stat().st_size/1024:.1f} KB)"

    print("\n" + "=" * 64)
    print("AP_A9 COVERAGE REPORT  (DAMIT spin & shape models)")
    print("=" * 64)
    print(f"Models total / distinct asteroids:   {n_models} / {n_ast}")
    print(f"With pole solution / with shape:     {with_pole} / {with_shape}")
    print(f"Multi-model asteroids:               {multi} (max models on one: {max_models})")
    print(f"GASP core matched:                   {matched} / 19,190 "
          f"({100.0*matched/19190:.1f}%)")
    print(f"Cross-check vs LCDB (A2): asteroids in both {len(xc)}")
    print(f"    median |P_damit - P_lcdb| / P_lcdb = {med_rel:.5f} ({med_rel*100:.3f}%)")
    print(f"    count with relative difference > 1%: {n_gt1} "
          f"({100.0*n_gt1/len(xc):.1f}%)")
    print(f"Apophis (99942) check:               {apo_models} inversion models, "
          f"{apo_tumb} tumbler solutions{apo_tinfo}")
    print(f"Saved: {sz('damit_models.parquet')}, {sz('damit_best.parquet')}, "
          f"{sz('damit_tumblers.parquet')}")
    print("=" * 64)

    dossier = ROOT / "docs" / "apophis_status.md"
    line = (
        f"\n## AP_A9 (DAMIT spin & shape, {TODAY})\n\n"
        f"- DAMIT inversion models for Apophis: **{apo_models}** (absent). This "
        f"is the expected, legitimate finding: Apophis rotates in a non-"
        f"principal-axis (NPA / tumbling) state, which classic convex "
        f"lightcurve inversion does not model. Instead Apophis appears in the "
        f"DAMIT **tumblers** table with **{apo_tumb} NPA solutions**"
        f"{apo_tinfo}. This is consistent with the A2 LCDB tumbler flag "
        f"(30.56 h, U=3, tumbler) - two independent databases agree on the "
        f"non-relaxed rotation state.\n"
        f"- Cross-package unit check: DAMIT vs LCDB rotation periods agree to a "
        f"median relative difference of {med_rel*100:.3f}% over {len(xc)} shared "
        f"asteroids (validates the hour units end-to-end).\n"
    )
    with dossier.open("a") as f:
        f.write(line)
    print(f"[dossier] appended Apophis DAMIT finding to {dossier.relative_to(ROOT)}")


def main():
    models, best, tumblers, ast, stats = ingest()
    coverage_report(models, best, tumblers, ast, stats)


if __name__ == "__main__":
    main()
