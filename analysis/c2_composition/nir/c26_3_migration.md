# C2.6 Step 3 - external path migrated to the C2.5 matcher

`match_external_spectrum()` now delegates to `c26_external` (C2.5 per-group shared-covariance matcher + RELAB-manifold gate), for both the full-band and NIR paths. **The RandomForest external path is retired for composition.** taxon_esox / C4 is a separate classifier and is untouched.

## NIR ground-truth under the migrated matcher (268 objects)

- within-manifold: **93.7%** (251/268)
- analog_status: unresolved 195, degenerate 68, **ok (confident) 5**
- the 5 confident are all HED: **2 diogenite, 2 howardite, 1 eucrite** - the distinctive pyroxene 1 um band clears the bar even in the truncated 770-1034 nm range.

**Key nuance: in-manifold != confident.** Unlike the Gaia core (out-of-distribution), NIR ground-truth is domain-valid (93.7% inside) but the NIR-7 space (truncated at 1034 nm, the Gaia grid edge) lacks the resolution to confidently discriminate most groups; only spectrally-distinctive HED clear the 0.55 threshold. Documented follow-up: matching external NIR spectra on the FULL MITHNEOS range (to ~2.5 um, beyond the Gaia 16-band grid) against full-range RELAB templates would raise NIR confidence materially - a separate matcher, not part of this migration.

## Apophis (migrated)

Path nir-6band: manifold d2 **0.54** (q95 13.41) -> **INSIDE** the manifold, leading indicated analog **ordinary_chondrite**, but calibrated confidence 0.18 < 0.55 -> honest status **indicative OC, unresolved** (OC ~0.13, ureilite ~0.13, CV ~0.12 near-tied). This SUPERSEDES the old RF line (conf 0.56), which was density-catch-all overconfidence. No confident analog is forced.
