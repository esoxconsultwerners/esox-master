# C4.4 - Validation & honesty

Independent held-out PDS Bus-DeMeo set (36 objects, never trained on). **Complex-level agreement: 80.6%** (leakage-free; C1.3's 86% was inflated by 100% train/val overlap - this is the trustworthy number). Among the 26 confident held-out predictions (conf >= 0.544), purity is **92.3%**. The agreement is carried by the dominant S complex; rare complexes (K/L, D) are the weak spots - see per-complex below and treat them with caution.

## Per-complex agreement (correct/total)

| complex | agree |
|---|---:|
| S | 24/25 |
| X | 0/1 |
| V | 4/4 |
| D | 1/1 |
| K/L | 0/5 |

## Apophis

taxon_esox = N/A for Apophis (99942): it has no Gaia spectrum (not in the GASP core), so no Gaia-feature prediction is produced. Its composition comes via the external-spectrum path (C1.4), not C4. The 44 nm Gaussian bandpass approximation (C1) applies here too.
