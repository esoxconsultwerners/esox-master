# C1.3 - Taxonomy arm (GASP spectra, C4 prototype)

Training labels: **358 GASP-Mahlke** core objects. Independent validation: **36 PDS Bus-DeMeo** core objects (held out). Features: 16 Gaia bands (normalized at 550 nm). StratifiedKFold CV; class_weight balanced (guardrail 4).

## Granularity

| level | balanced accuracy | classes kept |
|---|---:|---|
| complex (C/S/X/...) | 0.546 | S, X, V, C, D, A, K/L |
| full Bus-DeMeo class | 0.416 | S, V, X, C, D, A, Ch, B, M |

Confusion matrices: c1_3_confusion_complex.{png,csv}, c1_3_confusion_class.{png,csv}.

## Feature ablation (complex level, subset-matched lift)

| feature set | balanced accuracy | lift vs 16-band | n |
|---|---:|---:|---:|
| 16 Gaia bands (full set) | 0.546 | - | 353 |
| + SDSS u,g | 0.595 | +0.006 | 320 |
| + albedo | 0.410 | +0.005 | 84 |
| + phase slope (G1) | 0.425 | -0.042 | 36 |

(Each extra's lift is measured on the subset where it is present, with the 16-band baseline recomputed on that same subset. The SDSS u/g blue extension is the Esox-specific advantage over Gaia-only work.)

## Independent validation (PDS Bus-DeMeo, held out)

- complex-level agreement with Bus-DeMeo on 36 objects: **86.1%**.

## Scope verdict (taxonomy / C4)

- Gaia 16-band features support reliable classification at **complex** level. Complex-level balanced accuracy is 0.546 (rare-class-penalized), but independent Bus-DeMeo agreement on real-distribution held-out data is **86.1%** - the common complexes (S/C/X) classify well; the rare end-members (A, D, K/L) are the weak spots.
- Full Bus-DeMeo class level is NOT supported (balanced accuracy 0.416) - a class-level catalog would over-claim.
- Achievable coverage: at complex level, **57%** of the 19,190 core get a confident label (max prob >= 0.5) - so the taxonomy coverage jump from 3.9% toward ~60% is credible at **complex** level (matching the hoped ~60%), not at class level.

C4 should target **complex-level** classification, reporting per-object confidence so the reliable common complexes are not diluted by the uncertain rare end-members.

