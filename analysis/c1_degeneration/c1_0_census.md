# C1.0 - Census & band grid

## GASP / Gaia DR3 16-band grid

Centers (nm): [374, 418, 462, 506, 550, 594, 638, 682, 726, 770, 814, 858, 902, 946, 990, 1034]
Spacing: 44 nm uniform. Normalization reference: **550 nm** (refl_550 == 1.0 for all GASP spectra; RELAB is normalized identically in every stage).

## RELAB primary-group viability (>=8 spectra AND >=3 meteorites)

**23 primary groups survive** the viability filter.

| group | n_spectra | n_meteorites |
|---|---:|---:|
| CM | 285 | 126 |
| L | 185 | 149 |
| howardite | 152 | 101 |
| eucrite | 138 | 83 |
| H | 109 | 79 |
| LL | 89 | 74 |
| CV | 73 | 42 |
| ureilite | 73 | 46 |
| diogenite | 67 | 42 |
| lunar-meteorite | 63 | 18 |
| lodranite-acapulcoite | 41 | 26 |
| EL | 37 | 24 |
| EH | 29 | 18 |
| aubrite | 29 | 21 |
| pallasite | 22 | 4 |
| martian | 20 | 15 |
| CI | 16 | 6 |
| brachinite | 15 | 15 |
| CK | 13 | 9 |
| R | 13 | 10 |
| CR | 11 | 8 |
| CO | 11 | 11 |
| mesosiderite | 11 | 5 |

Primary groups excluded (below threshold): CB(5), CH(6), angrite(2), winonaite(4)

soft groups (L-LL, H-L) and coarse groups (C-ungrouped, iron, OC-ung, E) are reported separately in later stages and never counted as clean separability wins (guardrail 5).

## Taxonomy arm training set (GASP core, 19,190 objects)

GASP-Mahlke labels in core: **358** objects.

By complex:

| complex | n |
|---|---:|
| S | 129 |
| X | 61 |
| V | 61 |
| C | 53 |
| D | 25 |
| A | 13 |
| K/L | 11 |
| other | 5 |

By class (top):

| class | n |
|---|---:|
| S | 126 |
| V | 61 |
| X | 38 |
| C | 36 |
| D | 25 |
| A | 13 |
| Ch | 9 |
| B | 8 |
| M | 8 |
| Xk | 6 |
| K | 6 |
| Z | 5 |
| P | 5 |
| L | 5 |
| Q | 3 |
| Mk | 2 |
| Xe | 1 |
| Pk | 1 |

Independent validation labels in core (PDS, A5b):

- Bus-DeMeo: **42**  |  Bus: 273  |  Tholen: 145  |  S3OS2: 196

(The ~371-object PDS Bus-DeMeo set is mostly outside the GASP spectral core; only 42 overlap and serve as the held-out validation set.)

## Decision gate C1.0

- primary RELAB groups viable: 23 (need >= 6) -> PASS
- labeled GASP-core objects (Mahlke): 358 (need >= ~150) -> PASS

**GATE PASSED - proceed to C1.1.**

