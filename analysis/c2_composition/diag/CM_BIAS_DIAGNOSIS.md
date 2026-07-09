# CM-bias diagnosis (C2 composition matcher) — no pipeline change

**Dominant cause (one sentence):** The RandomForest matcher behaves as a template-*density* estimator — CM, the largest and spectrally-flattest separable group (278 RELAB spectra), blankets feature space so that low-feature Gaia spectra, pushed off the RELAB manifold by a large RELAB↔Gaia NIR red-slope mismatch, collapse onto CM as a catch-all; the C4 prior only masks this for the S/V/A minority while **unclassified (flat prior) and C-complex (CM is plausible) objects flow straight through**, which is why 2,208 of 2,286 confident analogs are CM.

The premise to test was "are S-complex objects getting confident CM?" — **the data falsifies it.** The prior is working; the leak is elsewhere.

---

## Check 1 — analog_top (confident) × taxon_esox complex

```
analog_top            A    C   D   S   V   X  unclassified
CM                    0  988  34   1   0  62          1123
brachinite            1    0   0   0   0   0             0
eucrite               0    0   0   0   1   0             0
howardite             0    0   0   0   5   0             0
ordinary_chondrite    0    0   0  69   1   0             0
ureilite              0    0   0   1   0   0             0
```

**No smoking gun where the brief expected one.** Only **1** S-complex object gets confident CM; S-complex confident analogs correctly go to `ordinary_chondrite` (69). The confident-CM flood is **unclassified 1,123** (flat prior — no down-weight) + **C-complex 988** (CM is a plausible C analog — no down-weight). The soft C4 prior suppresses S→CM as designed; it has no lever on these two populations.

Upstream evidence — raw flat-prior argmax across all 19,190 objects: **CM 12,590 (65.6%)**, ordinary_chondrite 2,255, EH 2,023, CI 965, EL 815, … The bias lives in the **likelihood itself**, before any prior.

## Check 2 — template count vs wins

`analysis/c2_composition/diag/check2_templates_vs_wins.png` · `check2_template_vs_wins.csv`

Pearson r(n_templates, times_confident) = **0.52** (all), 0.51 (separable). But it is not pure count:

| group | n_templates | times_confident |
|---|---:|---:|
| ordinary_chondrite | 383 | 70 |
| **CM** | **278** | **2,208** |
| howardite | 152 | 5 |
| eucrite | 138 | 1 |
| ureilite | 73 | 1 |

`ordinary_chondrite` has **more** templates than CM yet 70 confident wins vs CM's 2,208. So the driver is not template count alone — it is CM's combination of many templates **and** spectral flatness, which lets a single CM template sit near almost any low-contrast spectrum with high RF vote purity.

## Check 3 — centroid vs nearest-template (mechanism)

The requested "50 S-complex confident-CM" sample does not exist (only 1 such object — the prior removed the rest). Re-run on a deterministic 50-sample of the **actual** confident-CM population:

- CM is the **group-centroid-closest**: **5 / 50**
- CM is only the **nearest single template** (centroid NOT closest): **25 / 50**
- What a nearest-**centroid** matcher would pick instead: **aubrite 26, brachinite 16, CM 5, EH 2, CI 1**

**Verdict: mechanism (b) — density artifact.** CM wins because 278 templates give it more chances to own the nearest neighbour, not because the CM group centroid matches. A per-group aggregate likelihood demotes CM in **45 / 50** of these objects. (Note the alternatives are aubrite/brachinite — also flat groups — meaning these low-contrast objects are genuinely ambiguous and the honest label is *unresolved*, not confidently anything.)

## Check 4 — normalization / feature-scale

550 nm reference normalization is **consistent** (GASP `refl_550 ≡ 1.0`; RELAB normalized at 550). But two systematic issues:

1. **CM is the low-feature attractor.** Median band-contrast (std across 16 bands): confident-CM **0.132** vs rest **0.192**. CM templates are the flattest group (median contrast 0.074). Featureless spectra default to the featureless-and-densest group.
2. **RELAB↔Gaia NIR slope domain shift** — the real root. Per-band mean, RELAB template vs Gaia core:

   | band | RELAB | Gaia | Δ |
   |---:|---:|---:|---:|
   | 374 | 0.669 | 1.097 | **+0.428** |
   | 550 | 1.000 | 1.000 | 0.000 |
   | 858 | 0.957 | 1.145 | +0.188 |
   | 1034 | 0.967 | **1.582** | **+0.615** |

   RELAB templates roll over past ~770 nm; Gaia core spectra keep rising red to 1.58 at 1034 nm. The Gaia spectra sit **outside the RELAB training manifold**, and out-of-distribution points fall into whatever region dominates by density — CM. This is the same 44 nm Gaussian-bandpass approximation flagged as Werner open-task #2 (verify against real Gaia DR3 BP/RP response), now shown to have a downstream consequence.

## Check 5 — prior realized effect on CM

Where the prior applies it works but only masks: S-complex CM candidates have raw p(CM) mean **0.726** → post-prior **0.603** (down-weight ×0.25, renormalized) — still above the 0.55 threshold. Crucially the prior has **no effect** on the two dominant contributors: unclassified (flat prior) and C-complex (CM plausible → weight 1.0). The prior is a band-aid on a wound in the likelihood.

---

## Recommendation (Werner decides)

**Primary fix — switch to a per-group aggregate likelihood.** Replace the density-driven RF/nearest-template behaviour with a per-group centroid / Mahalanobis or per-group GMM likelihood, so groups cannot win by sheer template count. Check 3 shows this alone demotes CM in 45/50 confident cases; expect a large share of today's 2,208 confident-CM to fall to **degenerate/unresolved** — which is the honest C1 outcome, not a regression. Pair with template class-balancing (subsample or weight by 1/n_templates) to remove residual count bias.

**Necessary precondition — reconcile the RELAB↔Gaia NIR slope (Check 4) first.** Until templates and Gaia spectra share a slope convention (the Gauss-passband / BP-RP-response task), the matcher runs out-of-distribution and *any* classifier will find a catch-all. Fix this before re-running C2.

**Rejected candidates.** Hard C4-prior consistency filter — cannot touch the two dominant contributors (unclassified 1,123 has no prior; C-complex 988 finds CM plausible), so it barely moves the number. Tighter confidence threshold — cosmetic: CM raw p ≈ 0.73 sits well above any reasonable cut, and raising it discards the few genuine calls too.

Apophis is unaffected (NIR-7-band path, `ordinary_chondrite` conf 0.556) — the bias is a visible-band, template-density effect on the Gaia core path.
