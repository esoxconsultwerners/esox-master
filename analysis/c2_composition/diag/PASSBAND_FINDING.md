# C2.5 Step 0/1 — Passband acquisition + STOP finding

**Status: Step 0 done (authoritative passbands acquired & verified). Step 1 STOPPED before touching the resampler/matcher — the passband swap as scoped is (a) not executable with the authoritative file and (b) demonstrably not the cause of the diagnosed slope offset. Awaiting Werner decision.**

This is the brief's own Step-1 gate firing: *"If the offset does NOT shrink, the passbands aren't the whole story — stop and report before touching the matcher."* It does not shrink, and here is why, with evidence.

---

## What was acquired (Step 0 — success)

Official **Gaia (E)DR3 photometric passbands**, Riello+ 2021 (A&A 649, A3; CDS J/A+A/649/A3), *version 2* (revised). Downloaded from ESA Cosmos, checksummed, unzipped to `data/raw/gaia_passbands/`. Provenance: `data/provenance/gaia_passbands.json` (source URL, version, both SHA256s). The EDR3 passbands **are** the DR3 photometric system (carried from EDR3 into DR3) — this is the correct, non-DR2 file.

## Blocker 1 — granularity/category mismatch (fix not executable as written)

The brief says: *"convolve each RELAB spectrum with the actual BP/RP response per GASP band (replacing the 44 nm Gaussian)."* This is not possible with the authoritative file:

- `passband.dat` provides **three broadband integrated curves** — G, BP, RP — sampled 320–1100 nm.
- The 16 GASP "bands" (374…1034 nm, 44 nm spacing) are **not broadband photometry**. They are the sample points of the **Gaia DR3 SSO asteroid *reflectance spectrum*** product (Galluccio+2023, A&A 674, A35): a low-resolution BP/RP *spectrum* internally calibrated and divided by a solar analog, sampled at 16 wavelengths.
- Convolving a lab spectrum with the BP (or RP) broadband curve yields **one** number (the band flux), not 16. Three curves cannot produce a 16-element feature vector.

The correct per-sample instrument model for the 16-band product is the Gaia XP dispersion/sampling kernel (GaiaXPy) *plus* a replica of Galluccio's SSO reflectance pipeline — a separate research task, **not a drop-in replacement** for `gaussian_resample()`. Forcing the broadband curves into a per-band convolution would fabricate an instrument model, which the project rules forbid.

## Blocker 2 — the slope offset is NOT passband-caused (evidence)

The diagnosis's Check-4 offset (RELAB 0.967 vs Gaia 1.582 @1034 nm) compared the *whole RELAB library* against the *whole Gaia core* — different populations. I re-tested with **taxonomy-matched populations** (`c25_matched_offset.png/.csv`):

| matched pair | Gaia @1034 | RELAB analog @1034 | offset |
|---|---:|---:|---:|
| Gaia **S** vs RELAB **ordinary_chondrite** | 1.578 | 0.988 | **+0.590** |
| Gaia **C** vs RELAB **CM** | 1.436 | 1.035 | +0.401 |
| Gaia **V** vs RELAB **HED** | 1.367 | 0.855 | +0.511 |

**The offset persists at full strength within matched taxonomy** (+0.590 for S-vs-its-own-OC analog ≈ the +0.615 global figure). So it is *not* population composition — it is a **systematic Gaia-SSO-vs-RELAB domain shift**, present for every complex.

And it is **concentrated at the two extreme bands**: for S-vs-OC, mean |offset| is **0.026** across mid-bands (462–726 nm) but **0.516** at the edge bands (374, 1034) — a 20× concentration. Those two edge bands are exactly the ones Gaia DR3 SSO reflectance spectra are documented to reconstruct least reliably (BP/RP edge, solar-analog division). A resolution-kernel change (44 nm Gaussian → true kernel) on already-smooth reflectance **cannot** generate a +0.5 edge offset. Therefore swapping the passband would leave this offset essentially untouched — the Step-1 success criterion fails by construction.

## What this means for the CM bias

The domain shift is real (good — the diagnosis was right that Gaia lives off the RELAB manifold), but its cause is **edge-band SSO reconstruction + real surface effects (space weathering / phase reddening)**, not the bandpass shape. The actual cure for false-confident CM is therefore the **matcher change (per-group aggregate likelihood)** plus **handling the two unreliable edge bands** — not a passband swap.

## Recommendation (Werner decides — three options)

1. **(Recommended) Matcher-first, edge-band-aware.** Skip the passband rewrite. Do C2.5 Step 2 (per-group Mahalanobis/GMM likelihood — the diagnosed cure, demotes CM in 45/50) and additionally **down-weight or drop the 374 & 1034 nm edge bands** (or model the systematic lab→Gaia offset as a nuisance, like C1.2 weathering). This attacks both the density catch-all and the real domain shift, using only verified data.
2. **Full instrument forward-model (separate research task).** Reproduce the Gaia DR3 SSO XP→reflectance pipeline (GaiaXPy sampling + solar-analog division) to forward-model RELAB templates into the true 16-band basis. Scientifically ideal, weeks of work, own package — not part of this fix.
3. **Provide a real per-band SSO response** if one exists as an auxiliary DR3 product Werner can upload; then a genuine per-band convolution becomes possible. (The broadband Riello file is not it.)

I did **not** modify `gaussian_resample()`, the matcher, or any `analog_*` column. The acquired passbands + provenance + this finding are committed for the record.
