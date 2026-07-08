# C1.1 - Composition information ceiling

Target: 23 viable RELAB **primary** meteorite groups (1486 spectra, 932 distinct meteorites). GroupKFold by meteorite (guardrail 1); 550 nm normalization (guardrail 3); Gaussian-bandpass resampling to the 16 Gaia bands (guardrail 6, see method note in c1_common.py).

## Overall separability

| metric | noise-free (ceiling) | with Gaia noise |
|---|---:|---:|
| balanced accuracy (LDA) | 0.412 | 0.174 |
| balanced accuracy (RandomForest) | 0.332 | 0.174 |
| mean pairwise ROC-AUC | 0.912 | 0.796 |

## Per-class ROC-AUC (with Gaia noise, one-vs-rest)

| group | AUC (noisy) | pairwise-separable fraction (AUC>=0.75) |
|---|---:|---:|
| diogenite | 0.935 | 95% |
| howardite | 0.909 | 95% |
| lodranite-acapulcoite | 0.899 | 91% |
| CM | 0.895 | 77% |
| eucrite | 0.875 | 86% |
| brachinite | 0.870 | 68% |
| EL | 0.867 | 68% |
| CO | 0.855 | 50% |
| aubrite | 0.832 | 68% |
| EH | 0.792 | 55% |
| ureilite | 0.770 | 64% |
| martian | 0.767 | 77% |
| CR | 0.750 | 55% |
| CK | 0.749 | 36% |
| H | 0.743 | 59% |
| L | 0.733 | 73% |
| CV | 0.725 | 59% |
| LL | 0.724 | 55% |
| lunar-meteorite | 0.720 | 77% |
| R | 0.704 | 41% |
| pallasite | 0.594 | 36% |
| CI | 0.490 | 27% |
| mesosiderite | 0.441 | 5% |

## Pairwise separability (the scientific core)

- mean off-diagonal pairwise AUC: 0.796 (noisy), 0.912 (ceiling).
- groups separable (AUC>=0.75 vs a majority of others): **17/23** (74%).
- least separable groups (collapse at Gaia resolution): mesosiderite, CI, CK, pallasite, R.
- most separable: diogenite, howardite, lodranite-acapulcoite, eucrite, CM.

Full pair matrix in pairwise_auc.csv.

## SDSS u/g blue-extension lift

Adding synthetic SDSS u (354 nm) + g (477 nm) changes the mean pairwise AUC by **+0.029** (1486 spectra cover the blue). This quantifies the Esox-specific advantage over Gaia-only work.

## KILL CRITERION (composition)

Threshold: catalog framing supported if >=30% of primary groups are pairwise-separable (AUC>=0.75 vs most others) WITH Gaia noise. Observed: 74%.

**VERDICT: the data support a CATALOG paper (confident meteorite analogs for a well-defined subset).** A large fraction of primary groups stay distinct at 16-band Gaia resolution with noise; a matching pipeline is justified as a separate later package (not built here).

No matching pipeline is built in C1 regardless of verdict (by design).

