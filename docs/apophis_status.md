# Apophis (99942) — Showcase Dossier

Status of (99942) Apophis across the Esox Master Catalog ingestion
packages. The Apophis line is mandatory in every coverage report from
AP_A2 onward; this file collects the findings per package.

| Package | Source | Finding |
|---------|--------|---------|
| A1 | MPCORB + SBDB | Present in orbit backbone (NEO flag, PHA flag set) |
| A2 | LCDB V4.0 | Period **30.56 h**, U = 3, amplitude 0.30–1.14 mag, Notes `T` → **tumbler / NPA rotator** |
| A3 | Nesvorný families V2.0 | Absent, as expected (NEO; families are main-belt) |
| A4 | AKARI AcuA 1.0 | Absent (too small/faint for the survey) |
| A4 | IRAS SIMPS V6.0 | Absent (too small/faint for the survey) |
| A5 | SMASSII (VIS) | Absent |
| A5 | **MITHNEOS (NIR)** | **PRESENT: 3 measured spectra, 931 points, 0.770–2.485 µm** (`99942_*` in `data/interim/groundtruth_spectra.parquet`) — measured NIR ground truth for the showcase; literature class Sq (Binzel et al. 2009), not distributed in the PDS4 bundle (`gt_taxon` null, see provenance) |
| A5 | ECAS | Absent (1985 survey predates discovery) |
| A5b | PDS Asteroid Taxonomy V1.1 | absent (compilation predates wide NEO coverage; Sq label per Binzel et al. 2009 remains literature-only) |
| A6 | NEOCC + Sentry risk lists | Absent from both, as expected — removed 2021 after the radar campaign ruled out impacts for 100+ years (2004: Torino 4, highest ever measured) |
| A6 | JPL NHATS | PRESENT in NHATS: min dv 6.049 km/s (354 d), 272,327 viable trajectories; radar windows 2028-09 (Arecibo-class) / 2029-04 (Goldstone) archived in raw response |
| A7 | RELAB (lab reference set) | **91 LL-chondrite spectra** available as the laboratory counterpart to Apophis' literature Sq classification (Sq ↔ LL connection, Binzel et al. 2009) |

Next expected touchpoints: A6 (NHATS/NEOCC risk data — Apophis should
appear), Phase C1 (RELAB validation), Phase C4 (ML classifier ground
truth: the MITHNEOS NIR spectra are the anchor).

## AP_A8 (Meteoritical Bulletin + orbits, 2026-07-08)

- LL falls fraction: **8.6%** of classified MetBull falls (91/1063) are LL chondrites. This is the base rate behind the Apophis Sq<->LL story: Apophis is an Sq-type, and Sq is the spectral bridge to LL ordinary chondrites - the most probable meteorite analogue. LL is a minority of falls, so a delivered Apophis fragment would be a comparatively rare, diagnostic sample.
- Observed-fall orbits with published references: 60 events (incl. Chelyabinsk LL5, source region traceable per row).
