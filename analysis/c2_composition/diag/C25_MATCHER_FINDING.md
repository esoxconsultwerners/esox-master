# C2.5 Step 2 — Matcher fix + a deeper finding: Gaia↔RELAB domain divergence

**Bottom line:** The per-group matcher **eliminates the CM density catch-all** (CM 97% → ~0–5% of confident, S-complex→CM = 0) — the primary objective is met. But doing the matcher honestly exposed the *real* root cause, one level below both the passband and the matcher: **the Gaia DR3 asteroid-reflectance core and the RELAB lab-meteorite templates are structurally divergent distributions.** The Gaia core sits ~8× beyond the RELAB self-manifold, so an *honest* confident meteorite-analog set collapses to ~0.1%. The RF's 2,208 confident CMs were manufactured by an unconstrained density catch-all running far outside its training distribution.

**I did NOT overwrite the production `analog_*` columns.** This is a strategic call for Werner (options at the end): the honest result changes what the composition layer can claim.

---

## What was built (Step 2)

`analysis/c2_composition/c25_matcher.py` — a shared-covariance Gaussian classifier (per-group centroid + one pooled Ledoit-Wolf covariance → Mahalanobis score; uniform priors so template count cannot win; temperature-scaled softmax; isotonic confidence calibration; edge bands 374 & 1034 nm dropped). This is a clean, count-unbiased replacement for the RandomForest.

Three variants were tested, each teaching the next (`diag/c25_numbers.json`):

| matcher | self-CV bal.acc | ECE | confident (of core) | CM frac of confident |
|---|---:|---:|---:|---:|
| RF (production C2) | 0.332 | 0.042 | 2,286 (12%) | **97%** |
| per-group QDA (full cov) | 0.268 | 0.321 | 191 | 71% (unstable) |
| shared-cov + temp (relative) | 0.508 | 0.198 | 6,251 (33%) | **1.2%** |
| shared-cov + absolute-fit gate | 0.508 | — | **15 (0.1%)** | **0%** |

The shared-covariance matcher fixes CM (1.2%) and puts S-complex on stony analogs (ordinary_chondrite / lodranite / ureilite / EH), **never CM**. Reddening-augmentation (C1.2 weathering marginalization) was tried and rejected — it stripped the discriminative slope axis and inflated confident to 58% of the core, violating C1's honesty principle.

## The deeper finding — why "confident" is the wrong frame here

A per-group softmax only measures *relative* fit: it will confidently pick the nearest group even if the object matches **no** template. So I added an absolute goodness-of-fit gate (Mahalanobis d² vs the RELAB in-group χ² distribution). Result (`diag/c25_manifold_gap.png`):

- RELAB templates, distance to their **own** group centroid: mean d² = 6.8, q95 = 23.1.
- Gaia core, distance to the **nearest of any** group centroid: **median 188.8** (q25 85, q75 356).
- **Only 5.3% of the Gaia core falls within the RELAB q95 manifold.** 94.7% are, in absolute terms, unlike any meteorite template.

Is it a correctable systematic (one de-bias step)? **No.** Subtracting the global median Gaia−RELAB offset moves the core median d² only 188.8 → 155.5 (within-q95 *drops* to 2.6%). The divergence is high-dimensional and structural — space weathering (band-depth suppression) + phase reddening + the DR3 SSO reflectance reconstruction — not a rigid translation. This is the same class of issue as the passband finding (a real physical/instrumental domain shift), now shown to be **pervasive across all 13 retained bands**, not just the two edges.

## What this means

1. **CM bias: fixed.** Any principled per-group matcher kills it; the RF's 97%-CM was a catch-all artifact on out-of-distribution data. That part of the diagnosis is fully resolved.
2. **Confident Gaia-only meteorite analogs are not honestly supportable at the group level.** With an honest absolute-fit gate, ~0.1% of the core matches. Relaxing to relative confidence (1.2% CM, spread across stony groups) re-commits the RF's original sin — false confidence on OOD spectra — just distributed differently.
3. **The external-NIR path is unaffected.** Apophis (MITHNEOS NIR, not Gaia) still returns ordinary chondrite; ground-truth-spectrum matching lives on firmer ground than the Gaia core.

## Recommendation (Werner decides — three strategic options)

1. **(Recommended) Complex-level composition + honest unresolved core.** Ship the composition layer at the level that survives scrutiny: C4's complex taxonomy (57%, 81% held-out agreement) as the primary product; emit the meteorite-analog *distribution* as a soft indicator but mark the Gaia-only core `analog_status` overwhelmingly **unresolved/low_snr**, with confident group-level analogs reserved for objects with ground-truth NIR spectra (Apophis-style). Truthful, and it doesn't overclaim.
2. **Domain-adaptation research (the real Stage-1 science).** Build a Gaia↔RELAB space-weathering/phase forward-model that maps templates into the Gaia observational frame before matching. Weeks of work, own package — this is what would legitimately unlock per-object Gaia analogs. Park the catalog columns until then.
3. **NIR-only confident analog path.** Restrict confident meteorite analogs to the A5 ground-truth / MITHNEOS subset (where matching is defensible) and label the Gaia core complex-level-only.

Until Werner picks, the production `analog_*` columns are **unchanged** (still the RF C2 output). `c25_matcher.py` and this finding are committed for the record; no `esox_master_core.parquet` bytes were modified.
