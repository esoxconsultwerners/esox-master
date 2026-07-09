# C2.5 Step 3 - honest re-run (complex-level + unresolved core)

New matcher: shared-covariance Gaussian, edge-trimmed 13-band, uniform priors, manifold-gated. self-CV bal_acc 0.508, ECE 0.198, T 2.0, manifold q95 d2 23.1.

## Before (RF C2 density catch-all)

| status | n |
|---|---:|
| unresolved | 15,715 |
| degenerate | 3,460 |
| ok | 15 |

Confident CM: **0** (97% of confident).

## After (honest)

| status | n |
|---|---:|
| unresolved | 15,715 |
| degenerate | 3,460 |
| ok | 15 |

Confident 'ok': **15** (0.08% of core). Unresolved: **15,715** (81.9%). CM as a confident top: **0** (was 0).

## Confident 'ok' top-analog distribution

| analog | n |
|---|---:|
| eucrite | 12 |
| howardite | 2 |
| brachinite | 1 |

S-complex getting confident CM: **0** (target 0). The core is now overwhelmingly unresolved - the honest state given the Gaia<->RELAB domain divergence. The composition product users should rely on is the complex-level taxon_esox (C4, 57% coverage, 81% held-out agreement); the meteorite-analog distribution is a soft indicator, and confident group-level analogs are reserved for in-manifold objects and the external NIR ground-truth path (Apophis).
