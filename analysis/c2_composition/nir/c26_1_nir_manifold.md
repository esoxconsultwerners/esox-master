# C2.6 Step 1 - NIR ground-truth manifold test (MEASURED, not assumed)

Before pruning the NIR path with the manifold gate, we MEASURED where NIR ground-truth spectra fall relative to the RELAB manifold - in the NIR-7 feature space actually used by the external matcher (bands 770/814/858/902/946/990/1034 nm, normalized at 770, 6 log features). The 1034 nm edge-drop from the Gaia core matcher does NOT apply here: for real MITHNEOS NIR spectra 1034 nm is a legitimate measured point, not a Gaia-SSO reconstruction artifact.

- RELAB NIR self-distance: mean d2 3.35, **q95 = 13.41**.
- 268 NIR ground-truth objects: **93.7% INSIDE** the manifold (d2 <= q95), median d2 **1.5**.
- **Apophis: d2 0.54 -> INSIDE** (near the RELAB median).

**Contrast with the Gaia core: ~5% inside (median d2 188.8 vs q95 23.1, visible-band).** The NIR ground-truth path is therefore NOT the Gaia out-of-distribution case: lab-NIR and observed-NIR share wavelength range and absorption bands, and 93.7% of NIR ground-truth spectra are genuinely within the RELAB manifold. Confident analogs on the NIR path are domain-valid; the gate still prunes the ~6% that fall outside. Figure: nir_manifold.png. Matcher (NIR-6band): bal_acc 0.407, ECE 0.216, T 1.5.
