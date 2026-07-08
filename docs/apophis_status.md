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

- LL falls fraction: **8.0%** of classified MetBull falls (85/1063; pure LL - L/LL intermediates now split into the L-LL soft group) are LL chondrites. This is the base rate behind the Apophis Sq<->LL story: Apophis is an Sq-type, and Sq is the spectral bridge to LL ordinary chondrites - the most probable meteorite analogue. LL is a minority of falls, so a delivered Apophis fragment would be a comparatively rare, diagnostic sample.
- Observed-fall orbits with published references: 60 events (incl. Chelyabinsk LL5, source region traceable per row).

## AP_A9 (DAMIT spin & shape, 2026-07-08)

- DAMIT inversion models for Apophis: **0** (absent). This is the expected, legitimate finding: Apophis rotates in a non-principal-axis (NPA / tumbling) state, which classic convex lightcurve inversion does not model. Instead Apophis appears in the DAMIT **tumblers** table with **2 NPA solutions** (tumbler periods phi=27.38 h, psi=263.0 h). This is consistent with the A2 LCDB tumbler flag (30.56 h, U=3, tumbler) - two independent databases agree on the non-relaxed rotation state.
- Cross-package unit check: DAMIT vs LCDB rotation periods agree to a median relative difference of 0.028% over 5582 shared asteroids (validates the hour units end-to-end).

## Assembled dossier v0 (Phase B merge)

The showcase record every A1-A9 layer has been building toward - (99942) Apophis, assembled from the merged master catalog:

- **Orbit (A1):** a = 0.922 au, e = 0.191, i = 3.34 deg; SBDB class ATE, NEO=True, PHA=True, MOID = 0.000108 au. H = 19.1 [SBDB].
- **Rotation (A2 + A9):** non-principal-axis **tumbler** (tumbler_flag=True) - so period_best is deliberately null (Rule 1.2). LCDB reports ~30.56 h (U=3, note 'T'); DAMIT carries the NPA solution in its tumblers table. Two independent databases agree on the non-relaxed rotation state.
- **Size/albedo (A4/SBDB):** diameter_best = 0.34 km [NEOWISE], albedo_best = 0.35 [NEOWISE].
- **Accessibility & risk (A6):** NHATS min dv = 6.05 km/s (2.72e+05 trajectories); NEOCC/Sentry impact probability = -/- (removed from risk lists after the 2004-2021 observation arc).
- **Taxonomy (literature):** mahlke=-, tholen=-, bus=-, bus_demeo=-, consensus=- - S-complex / Sq, the spectral bridge to LL ordinary chondrites.
- **Spectra (A5 + GASP):** gaia_spectrum=False; MITHNEOS NIR spectra provide the 0.77-2.49 um coverage. The laboratory analogue is the RELAB LL-chondrite set (A7).
- **Complete physical profile:** False (a tumbler has no single period_best by design, so Apophis is intentionally not a 'complete profile' object - the flag behaves correctly).

## Composition analog (external-spectrum path)

Via C1.4: Apophis has no Gaia spectrum, so its composition analog is recovered from MITHNEOS NIR spectra resampled to the 7 Gaia bands >= 770 nm. Top analog **L** (p=0.23); ordinary-chondrite mass 0.54, LL p=0.18. Consistent with the literature Sq -> LL story; the external-spectrum showcase mechanism works. (Lower bound: no visible slope, no 2 um band at Gaia resolution.)


- **taxon_esox = N/A** (C4): Apophis has no Gaia spectrum, so the model-derived complex classifier does not run for it. Its composition comes only via the external-spectrum path (C1.4). This is the correct behaviour, not a gap.

## Composition analog (production)

Via C2 (external-spectrum path, nir-7band): Apophis analog distribution top **ordinary_chondrite** (conf 0.56), ordinary-chondrite mass 0.56. Per C1 honesty, the ordinary-chondrite subgroup is reported **unresolved** (H/L/LL are not separable at Gaia/NIR resolution) - Apophis gets 'ordinary chondrite, subgroup unresolved', not 'LL'.

