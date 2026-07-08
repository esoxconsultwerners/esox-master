# C4.1 - Train & calibrate

Features: 16 Gaia bands only (C1.3 showed SDSS u/g adds only +0.006 at complex level, so the classifier is kept purely spectral). Complex scheme (from C1.3): ['S', 'C', 'X', 'V', 'D', 'A', 'K/L']. Training 311 GASP-Mahlke objects, leakage-free (the 36 PDS Bus-DeMeo objects are excluded - they are 100% inside the Mahlke set, so C1.3's 86% was leakage-inflated).

## Model selection (CV balanced accuracy)

| model | balanced acc |
|---|---:|
| LogReg(balanced)+sigmoid | 0.434 |
| LogReg+sigmoid | 0.438 |
| LogReg+isotonic | 0.494 |
| RF+isotonic | 0.525 |

Chosen: **LogReg+isotonic**. Calibrated (sigmoid) balanced acc 0.494.

## Calibration

Reliability diagram: calibration.png. Expected calibration error **ECE = 0.050** (lower is better; the calibrated confidence is trustworthy where the curve tracks the diagonal).

## Operating point

Confidence threshold **0.544**, chosen to reproduce the C1.3 ~57% confident-coverage point. Achieved coverage 57.0%.
