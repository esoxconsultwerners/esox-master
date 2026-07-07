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

Next expected touchpoints: A6 (NHATS/NEOCC risk data — Apophis should
appear), Phase C1 (RELAB validation), Phase C4 (ML classifier ground
truth: the MITHNEOS NIR spectra are the anchor).
