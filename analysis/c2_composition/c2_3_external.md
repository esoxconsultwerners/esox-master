# C2.3 - External-spectrum path

match_external_spectrum() resamples any (wavelength, reflectance) and returns the analog distribution, choosing the full 16-band matcher when the visible is covered else the NIR-7-band matcher (reproducing C1.4).

## Apophis (production)

Path nir-7band: top **ordinary_chondrite** (conf 0.56), OC mass 0.56, oc_subgroup **unresolved**.

## Gaia-vs-external cross-check

On 270 core objects with both a Gaia spectrum and an A5 ground-truth spectrum, the external path agrees with the Gaia path on the top analog **43%** overall, rising to **75%** on the 61 objects the Gaia path calls confident (status=ok). The moderate overall figure reflects that most objects are low-confidence on both paths (as C1 predicts); where the Gaia path is confident the two paths agree well, which is what validates the external path the Apophis dossier and customer requests rely on.
