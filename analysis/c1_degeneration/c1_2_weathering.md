# C1.2 - Space weathering as a nuisance parameter

Reddening model: **W(lambda) = exp(-C / lambda)** anchored at 550 nm; strength swept over C = [0.0, 0.1, 0.2, 0.3] (um). Marginalized into the classifier by training across all weathering strengths (weathered copies of a meteorite stay in its GroupKFold fold - no leakage).

Citation: **TODO: verified reference (exponential reddening continuum)**

## Separability under weathering (mean pairwise ROC-AUC, Gaia noise)

| | mean pairwise AUC |
|---|---:|
| weathering-free | 0.855 |
| weathering-marginalized | 0.885 |
| **mean change** | **+0.030** |

Mean separability is essentially unchanged: marginalizing weathering acts as data augmentation and does not collapse the average. The real signal is per-pair: among the 45 pairs that degrade, the mean drop is -0.030 and the worst is CO vs EL at -0.225.

## Most weathering-sensitive group pairs (largest AUC drop)

| pair | AUC before | AUC after | delta |
|---|---:|---:|---:|
| CO vs EL | 0.841 | 0.616 | -0.225 |
| CR vs R | 0.717 | 0.533 | -0.183 |
| CK vs aubrite | 0.986 | 0.857 | -0.129 |
| L vs lodranite-acapulcoite | 0.838 | 0.734 | -0.104 |
| CV vs ureilite | 0.751 | 0.662 | -0.089 |
| howardite vs martian | 0.903 | 0.837 | -0.066 |
| EH vs ureilite | 0.861 | 0.814 | -0.047 |
| R vs martian | 0.878 | 0.837 | -0.042 |

Full table in c1_2_auc_delta.csv.

## Estimated weathering strength (would-be catalog column)

Per-spectrum reddening slope C, fit as ln(refl) vs (1/lambda - 1/0.55). This is a derivable science product (a per-object weathering-strength column):

- median C = 0.109, IQR [0.003, 0.269], range [-0.540, 1.070] over 1486 spectra.

## Summary

Mean pairwise separability is robust to weathering once marginalized (change +0.030 AUC - augmentation compensates). The honest cost is concentrated in specific weathering-ambiguous pairs (worst CO vs EL -0.225); those are the pairs a composition catalog must flag rather than assert.

