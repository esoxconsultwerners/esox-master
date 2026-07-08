# Esox Master Catalog — Precedence Rules for Conflicting Values

## Purpose

When multiple A-series sources provide a value for the same physical property
of the same object, these rules define **which source populates the derived
`*_best` column**. Original per-source values are **never overwritten**: every
`*_best` column has a paired `*_source` column that names the winning source
for that object, and all contributing source values are retained in their own
prefixed columns (`akari_*`, `sbdb_*`, `iras_*`, `lcdb_*`, `damit_*`, …).

This file is the **authoritative input to `pipeline/build_master.py`**: the
merge reads these rules as ground truth and executes them mechanically. It is
also written to be **directly citable as the methods basis** in the catalog
data paper — each rule states its decision, its rationale, and its verified
citation where an external one applies.

---

## Rule 1.1 — Diameter and albedo (NEOWISE/SBDB, AKARI, IRAS)

**Decision.** Derive `diameter_best` and `albedo_best` by a fixed source
precedence **AKARI > WISE/NEOWISE > IRAS**. Do **not** average across sources.
Retain all source values in their `akari_*`, `sbdb_*` (NEOWISE), and `iras_*`
columns; add `diameter_source` and `albedo_source` naming the winning survey
per object.

**Rationale.** This precedence, and the decision not to average, follow the
published merged-catalog convention of Usui et al. (2014, PASJ 66, 56), who
compared the three infrared surveys over 1993 commonly detected asteroids,
found agreement within ±10 % in diameter and ±22 % in albedo (1σ), and
prioritized AKARI over WISE over IRAS because AKARI carries less uncertainty
for the largest asteroids (no detector saturation on bright targets).

Known limitation: for small / dark objects NEOWISE has the larger, more
sensitive sample; `diameter_source` is therefore recorded **per object** so
that any later stage can override the default on a per-object basis. Albedo
accuracy for all surveys is limited by the ~0.3 mag dispersion in MPC *H*
values (documented as the dominant albedo error term in the literature), so
`albedo_best` inherits that floor regardless of source.

**Citation.** Usui, F., Hasegawa, S., Ishiguro, M., Müller, T. G., Ootsubo, T.
(2014). A comparative study of infrared asteroid surveys: IRAS, AKARI, and
WISE. *PASJ* 66(3), 56. DOI [10.1093/pasj/psu037](https://doi.org/10.1093/pasj/psu037).

---

## Rule 1.2 — Rotation period (DAMIT vs LCDB)

**Decision.** `period_best` = DAMIT inversion period where a DAMIT model
exists; otherwise the LCDB period with quality **U ≥ 2**. On the ~14 % of
shared objects where the two disagree (period aliases: ½×, 2×, pole-mirror
fitted-period differences), **DAMIT wins**. `period_source` records `DAMIT` or
`LCDB`. The LCDB value is always retained in its `lcdb_*` column.
Non-principal-axis rotators (tumblers) are flagged and **excluded from
single-period best selection** — for those, both the LCDB tumbler entry and the
DAMIT tumbler solutions are kept **without** a `period_best` (Apophis is the
reference case).

**Rationale.** DAMIT is treated in the literature as the most accurate list of
ground-based period and pole determinations (e.g. the Gaia DR3 Solar System
survey uses it as the reference standard); the LCDB itself deliberately does
not ingest DAMIT results, so the two are genuinely independent sources rather
than a circular reference. The A2×A9 cross-check in this project found a median
relative period difference of **0.028 % over 5,582 shared objects**, confirming
end-to-end unit and parsing consistency; the residual tail is astrophysical
aliasing, not error.

**Citation.** Ďurech, J., Sidorin, V., Kaasalainen, M. (2010). DAMIT: a
database of asteroid models. *A&A* 513, A46.
DOI [10.1051/0004-6361/200912693](https://doi.org/10.1051/0004-6361/200912693).

---

## Rule 1.3 — Taxonomy (Mahlke vs PDS multi-system)

**Decision.** Do **not** fuse. Keep `taxon_mahlke` (from GASP) and the PDS
per-system columns (`taxon_tholen`, `taxon_bus`, `taxon_bus_demeo`,
`taxon_s3os2` from A5b) side by side. Derive `taxon_literature_consensus`
**only** where two or more independent systems agree at the simple
complex / letter level; leave it null otherwise. The future ML classifier
(Phase C4) will add an independent `taxon_esox` column; three independent views
are retained by design rather than collapsed.

**Rationale.** Taxonomic systems are complementary frameworks, not
interchangeable measurements; the A5b cross-check found **83.3 % agreement**
between Mahlke and Bus-DeMeo, with disagreements concentrated at
adjacent-complex boundaries (X/P, X/C, Ch/B) — i.e. inter-system drift rather
than error. Forcing a single class would destroy information the boundaries
themselves carry.

*No external citation required (project-internal evidence).*

---

## Rule 1.4 — Absolute magnitude H

**Decision.** `h_best` = SBDB (JPL-maintained consolidated value); MPCORB *H*
as fallback. `h_source` records the origin. Both retained.

**Rationale.** JPL SBDB curates the consolidated *H*; trivial precedence.

---

## Rule 1.5 — Canonical object key

**Decision.** `object_key` = `str(number_mp)` for numbered objects, else the
principal designation normalized to MPC spacing. This is the single join key
across all layers. Establishes as a project-wide standard the normalization
already implemented in A6 (compact `"2023VD3"` and spaced `"2023 VD3"`
normalize identically).

**Rationale.** Consistent joins across numbered and unnumbered objects;
prerequisite for later LSST / DR4 matching.

---

## Rule 1.6 — Provenance granularity

**Decision.** Every derived `*_best` field carries a paired `*_source` column;
no `*_best` value may be non-null while its `*_source` is null (enforced as a
validation gate in `build_master.py`). This file itself is the global
precedence record and is versioned in git.

---

## Realized outcomes

*Filled by `pipeline/build_master.py` after the merge (counts over the full catalog).*

| Best field | Source | Objects (n) |
|---|---|---|
| `diameter_best` | AKARI | 5,086 |
| `diameter_best` | NEOWISE | 130,521 |
| `diameter_best` | IRAS | 0 |
| `albedo_best` | AKARI | 5,086 |
| `albedo_best` | NEOWISE | 129,600 |
| `albedo_best` | IRAS | 0 |
| `period_best` | DAMIT | 10,741 |
| `period_best` | LCDB | 22,687 |
| `period_best` | (tumbler, no best) | 148 |
| `h_best` | SBDB | 895,910 |
| `h_best` | MPCORB | 656,758 |
| `taxon_literature_consensus` | >=2 systems agree | 671 |
