# Esox Master Catalog

Successor infrastructure to [GASP](https://doi.org/10.3847/2515-5172/ae5e45).
Builds a single, versioned, citable asteroid master catalog by ingesting
independent public sources into one spine (the MPC orbit catalog) and
enriching it package by package.

Maintained by Esox Consult GmbH, Vienna — https://esoxspace.com

## Ingestion contract

Every ingestion package ships as `pipeline/ingest_<source>.py` and must:

1. **Download** its source and cache the raw file(s) under `data/raw/`
   (dated filenames, SHA256 recorded). Static, versioned, citable
   downloads are preferred over APIs.
2. **Normalize keys** — the join key is `number_mp` (IAU minor planet
   number, nullable Int64).
3. **Prefix columns** per source (e.g. `mpc_`, `sbdb_`, `lcdb_`) so
   provenance of every value is visible in the column name.
4. **Write provenance** JSON to `data/provenance/<source>.json`
   (url, retrieval date, sha256 where applicable, row counts, citation,
   license note).
5. **Have a contract test** under `tests/test_<source>.py` (pytest).
6. **Print a coverage report** against the GASP core
   (19,190 objects, key column `number_mp`,
   reference list: `data/interim/gasp_core_keys.parquet`).

Raw source files are cached locally and never redistributed; the project
publishes only derived tables.

## Packages

| # | Package | Source(s) | Status |
|---|---------|-----------|--------|
| A1 | Orbit backbone | MPC MPCORB extended + JPL SBDB | done |
| A2 | Rotation periods | LCDB (PDS4 bundle V4.0, doi:10.26033/j3xc-3359) | done |
| A3 | Asteroid families | Nesvorný HCM PDS4 V2.0 (doi:10.26033/5hyq-6k90) | done |
| A4 | Diameters & albedos | AKARI / IRAS (SIMPS) | planned |
| A5 | Reflectance spectra | SMASS / MITHNEOS + ECAS | planned |
| A6 | Phase curves | GAPC (HG parameters) | planned |
| A7 | Taxonomy | Bus–DeMeo / SDSS / Gaia | planned |
| A8 | Binaries, pairs & satellites | Johnston / literature | planned |
| A9 | Accessibility & close approaches | NHATS / CNEOS | planned |

(Planned names are the current roadmap and may be refined per package.)

## Layout

```
pipeline/          ingestion scripts (one per source)
data/raw/          cached raw downloads (not in git)
data/interim/      normalized per-source tables (not in git)
data/final/        merged master catalog releases (not in git)
data/provenance/   provenance JSON per source (in git)
tests/             contract tests (pytest)
```

## Usage

```
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python pipeline/ingest_mpcorb.py
.venv/bin/pytest tests/ -v
```
