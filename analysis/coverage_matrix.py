#!/usr/bin/env python3
"""CM v0 - Coverage matrix + completeness funnel figures.

Pure analysis on the Phase B core catalog (data/final/esox_master_core.parquet,
19,190 objects). The per-object coverage flags were already computed in the
merge; this script only tabulates and visualizes them - it does not recompute
precedence. Two figures back the Phi-Lab application and the website:
  1. coverage_matrix_v0  - property coverage bar chart (the coverage figure)
  2. completeness_funnel_v0 - the AND-cut waterfall that explains the 204
     headline (taxonomy is the binding constraint -> motivates the C4 ML
     classifier).
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "data" / "final"
CORE = FINAL / "esox_master_core.parquet"

# Minimalist-scientific palette (restrained, brand-consistent).
INK = "#2b2b2b"
GREY = "#b9bcc0"
ACCENT = "#c1666b"       # the taxonomy bottleneck
COMPLETE = "#4059ad"     # the complete-profile / final survivor
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.edgecolor": "#888888",
    "svg.hashsalt": "esox-cm-v0",
})

# property -> (flag/value column, is_flag, source package label)
PROPS = [
    ("orbit", "mpc_a", False, "A1 MPC/SBDB"),
    ("h_best", "h_best", False, "A1 SBDB>MPCORB"),
    ("phase_curve", "has_phase_curve", True, "GAPC (HG1G2)"),
    ("gaia_spectrum", "has_gaia_spectrum", True, "GASP (Gaia DR3)"),
    ("diameter_best", "has_diameter_best", True, "A4 AKARI>NEOWISE>IRAS"),
    ("albedo_best", "albedo_best", False, "A4 AKARI>NEOWISE>IRAS"),
    ("family", "has_family", True, "A3 Nesvorny"),
    ("period_best", "has_period_best", True, "A2/A9 LCDB/DAMIT"),
    ("taxon_any", "has_taxon_any", True, "A5b PDS + GASP-Mahlke"),
    ("groundtruth_spectrum", "has_groundtruth_spectrum", True, "A5 SMASS/MITHNEOS/ECAS"),
]
FUNNEL = [("spectrum", "has_gaia_spectrum"),
          ("+ diameter_best", "has_diameter_best"),
          ("+ family", "has_family"),
          ("+ period_best", "has_period_best"),
          ("+ taxon_any", "has_taxon_any")]


def count(core, col, is_flag):
    return int(core[col].astype(bool).sum()) if is_flag else int(core[col].notna().sum())


def build_table(core):
    n = len(core)
    rows = []
    for name, col, is_flag, pkg in PROPS:
        c = count(core, col, is_flag)
        rows.append({"property": name, "n_core": c, "pct_core": round(100 * c / n, 1),
                     "source_package": pkg, "note": ""})
    cp = int(core["complete_profile"].sum())
    rows.append({"property": "complete_physical_profile", "n_core": cp,
                 "pct_core": round(100 * cp / n, 1), "source_package": "derived (Phase B)",
                 "note": "gaia_spectrum AND period_best AND family AND diameter_best "
                         "AND taxon_any"})
    return pd.DataFrame(rows)


def funnel_counts(core):
    mask = pd.Series(True, index=core.index)
    out = []
    for label, col in FUNNEL:
        mask = mask & core[col].astype(bool)
        out.append((label, int(mask.sum())))
    return out


def write_table(df):
    df.to_csv(FINAL / "coverage_matrix_v0.csv", index=False)
    lines = ["# Esox Master Catalog v0 - coverage matrix", "",
             "Property coverage over the 19,190-object GASP spectral core "
             "(retrieval July 2026).", "",
             "| Property | n (core) | % core | Source package | Note |",
             "|---|---:|---:|---|---|"]
    for _, r in df.iterrows():
        lines.append(f"| {r['property']} | {r['n_core']:,} | {r['pct_core']:.1f}% | "
                     f"{r['source_package']} | {r['note']} |")
    (FINAL / "coverage_matrix_v0.md").write_text("\n".join(lines) + "\n")


def fig_coverage(df):
    d = df.sort_values("pct_core", ascending=True).reset_index(drop=True)
    colors = []
    for p in d["property"]:
        colors.append(ACCENT if p == "taxon_any"
                      else COMPLETE if p == "complete_physical_profile" else GREY)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = range(len(d))
    ax.barh(list(y), d["pct_core"], color=colors, height=0.62, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(d["property"])
    for yi, (pct, n) in enumerate(zip(d["pct_core"], d["n_core"])):
        ax.text(pct + 1.2, yi, f"{pct:.1f}%  ({n:,})", va="center", fontsize=8, color=INK)
    ax.set_xlim(0, 100)
    ax.set_xlabel("coverage of the 19,190-object spectral core (%)")
    ax.set_title("Esox Master Catalog v0 - property coverage of the\n"
                 "19,190-object spectral core", fontsize=11, color=INK, loc="left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.xaxis.grid(True, color="#e6e6e6", zorder=0)
    ax.set_axisbelow(True)
    ti = list(d["property"]).index("taxon_any")
    ax.annotate("bottleneck -> ML classifier (in progress)",
                xy=(d["pct_core"][ti] + 2, ti), xytext=(44, 0.75),
                fontsize=8, color=ACCENT, va="center",
                arrowprops=dict(arrowstyle="->", color=ACCENT, lw=0.8))
    fig.text(0.01, -0.02,
             "Sources: A1 MPC/JPL-SBDB; A2/A9 LCDB/DAMIT; A3 Nesvorny families; "
             "A4 AKARI/NEOWISE/IRAS; A5/A5b spectra & taxonomy; GASP (Gaia DR3); "
             "GAPC phase curves. Retrieval July 2026.",
             fontsize=6.2, color="#666666")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FINAL / f"coverage_matrix_v0.{ext}", dpi=150,
                    bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def fig_funnel(fc):
    labels = [l for l, _ in fc]
    vals = [v for _, v in fc]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    y = list(range(len(fc)))[::-1]
    colors = [GREY] * (len(fc) - 1) + [COMPLETE]
    ax.barh(y, vals, color=colors, height=0.6, zorder=3)
    for yi, v, lab in zip(y, vals, labels):
        ax.text(v + max(vals) * 0.012, yi, f"{v:,}", va="center", fontsize=9, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, max(vals) * 1.12)
    ax.set_xlabel("objects surviving each cumulative AND-requirement")
    ax.set_title("Why only 204 complete profiles - the completeness funnel\n"
                 "(19,190-object spectral core)", fontsize=11, color=INK, loc="left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.xaxis.grid(True, color="#e6e6e6", zorder=0)
    ax.set_axisbelow(True)
    drop = vals[-2] - vals[-1]
    ax.annotate(f"taxonomy cut: -{drop:,}\n(the binding constraint)",
                xy=(vals[-1], y[-1]), xytext=(max(vals) * 0.30, y[-1] + 0.35),
                fontsize=8.5, color=ACCENT,
                arrowprops=dict(arrowstyle="->", color=ACCENT, lw=0.9))
    fig.text(0.01, -0.02,
             "Each step adds one AND-requirement to the previous survivors; "
             "order ends on taxonomy, the binding constraint. Retrieval July 2026.",
             fontsize=6.2, color="#666666")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FINAL / f"completeness_funnel_v0.{ext}", dpi=150,
                    bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def apophis_check():
    cols = ["object_key", "has_gaia_spectrum", "has_diameter_best", "has_family",
            "has_period_best", "has_taxon_any", "tumbler_flag"]
    full = pd.read_parquet(FINAL / "esox_master_full.parquet", columns=cols)
    r = full[full["object_key"] == "99942"]
    print("\n=== Apophis (99942) funnel check (full-catalog record) ===")
    if r.empty:
        print("  Apophis not found in full catalog")
        return
    r = r.iloc[0]
    survived = True
    for label, col in FUNNEL:
        ok = bool(r[col])
        survived = survived and ok
        mark = "PASS" if ok else "DROP"
        extra = ""
        if col == "has_period_best" and not ok:
            extra = f" (tumbler_flag={bool(r['tumbler_flag'])} -> no single period, Rule 1.2)"
        print(f"  {label:<18} {mark}{extra}")
        if not ok:
            print(f"  -> Apophis first drops at '{label}'.")
            break
    print("  Note: Apophis is not in the 19,190 spectral core (no Gaia DR3 "
          "spectrum), so it drops at the spectrum step by data; it ALSO fails "
          "period_best by physics (NPA tumbler). Incomplete either way - the "
          "showcase object is incomplete by physics, not only by data gaps.")


def main():
    core = pd.read_parquet(CORE)
    df = build_table(core)
    write_table(df)
    fig_coverage(df)
    fc = funnel_counts(core)
    fig_funnel(fc)

    print(df.to_string(index=False))
    print("\nFunnel survivor counts:")
    for label, v in fc:
        print(f"  {label:<18} {v:,}")
    apophis_check()
    print("\nSaved figures:")
    for f in ("coverage_matrix_v0.png", "coverage_matrix_v0.svg",
              "completeness_funnel_v0.png", "completeness_funnel_v0.svg",
              "coverage_matrix_v0.csv", "coverage_matrix_v0.md"):
        p = FINAL / f
        print(f"  {p}  ({p.stat().st_size/1024:.1f} KB)")
    print("\n=== CM v0 complete ===")
    print("Coverage figure + completeness funnel ready.")
    cp = int(core["complete_profile"].sum())
    print(f"Headline: {cp}/{len(core):,} complete profiles; taxonomy the binding "
          f"constraint at {100*count(core,'has_taxon_any',True)/len(core):.1f}%.")
    print("Ready for C1 (degeneration analysis).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
