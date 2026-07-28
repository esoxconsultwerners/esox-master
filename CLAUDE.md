# CLAUDE.md — Esox Master Catalog

## Projektzweck

Nachfolge-Infrastruktur zu [GASP](https://doi.org/10.3847/2515-5172/ae5e45).
Baut **einen** versionierten, zitierbaren Asteroiden-Master-Katalog: unabhängige
öffentliche Quellen werden in ein Rückgrat (MPC-Bahnkatalog) eingespeist und
Paket für Paket angereichert. Maintainer: Esox Consult GmbH, Wien — https://esoxspace.com

**Dieses Repo ist öffentlich (MIT).** Es werden nur *abgeleitete* Tabellen
publiziert; rohe Quelldateien bleiben lokal (`data/raw/`) und werden nie
weiterverteilt. **Keine internen Markt-/Preis-/Strategie-Dokumente in dieses
Repo committen** — auch nicht gitignored (ein versehentliches `git add -f`
leakt sie). Private Strategie gehört in ein privates Repo.

---

## Single Source of Truth

- **Projektstand:** `docs/PROJECT_STATUS.md` — der ehrliche Gesamtstand
  (Composition-Realität nach C2.5/C2.6). Bei Statusfragen **immer zuerst hier
  lesen**, nicht aus dem Gedächtnis antworten. Keine konkurrierende Zweit-SSOT
  anlegen; dieses Dokument fortschreiben.
- **Apophis-Dossier:** `docs/apophis_status.md` — Showcase-Status, ehrlich.
- **Merge-Regeln:** `data/provenance/precedence_rules.md` — autoritativer Input
  für `pipeline/build_master.py`, zugleich zitierbare Methods-Basis fürs Paper.

---

## Ingestion-Contract (aus README, verbindlich)

Jedes Paket ist `pipeline/ingest_<source>.py` und muss:

1. **Download** der Quelle, Rohdaten unter `data/raw/` cachen (datierte
   Dateinamen, SHA256 protokolliert). Statische, versionierte, zitierbare
   Downloads vor APIs bevorzugen.
2. **Keys normalisieren** — Join-Key ist `number_mp` (IAU minor planet number,
   nullable Int64).
3. **Spalten prefixen** pro Quelle (`mpc_`, `sbdb_`, `lcdb_`, `akari_`, …), damit
   die Provenance jedes Werts im Spaltennamen sichtbar ist.
4. **Provenance schreiben** nach `data/provenance/<source>.json` (url, Abrufdatum,
   sha256, row counts, citation, license note).
5. **Contract-Test** unter `tests/test_<source>.py` (pytest).
6. **Coverage-Report** gegen den GASP-Kern (19.190 Objekte, key `number_mp`,
   Referenz `data/interim/gasp_core_keys.parquet`).

**`*_best`-Prinzip:** Original-Quellwerte werden nie überschrieben. Jede
`*_best`-Spalte hat eine gepaarte `*_source`-Spalte, die die Gewinnerquelle pro
Objekt nennt; alle Beiträger bleiben in ihren prefixed Spalten erhalten.

---

## Wissenschaftliches Leitprinzip

**Ehrlichkeit vor großen Zahlen.** Jede konfidente Aussage muss
mannigfaltigkeitsgedeckt sein. Kein Matcher fabriziert Analoge außerhalb seiner
Gültigkeit (C2.5-Befund: Gaia-DR3- und RELAB-Verteilungen divergieren
strukturell — 95 % des Gaia-Kerns liegen außerhalb der RELAB-Mannigfaltigkeit).
Das Composition-Produkt ist die **C4-Komplex-Taxonomie** plus konfidente
Analoge nur dort, wo domänen-gültig (NIR-Pfad). Details: `docs/PROJECT_STATUS.md`.

---

## Struktur

| Pfad | Inhalt |
|---|---|
| `pipeline/` | `ingest_*.py` pro Quelle + `build_master.py` (liest precedence_rules) |
| `analysis/` | `c1_degeneration/`, `c2_composition/` (`diag/`, `nir/`), `c4_classifier/` |
| `data/raw/` | roh, lokal, nie in git (`.gitignore`) |
| `data/interim/` | gitignored außer den freigegebenen CSVs (`!`-Regeln in `.gitignore`) |
| `data/final/*.parquet` | zu groß für git; öffentlich über Zenodo-DOIs |
| `data/provenance/` | pro Quelle JSON + `precedence_rules.md` |
| `tests/` | pytest Contract-Tests pro Quelle |
| `docs/` | Status-Dokumente (SSOT) |

---

## Stack & Konventionen

- **Python 3.13**, `.venv`. Deps: `requirements.txt` (pandas ≥2, numpy, pyarrow,
  scikit-learn, scipy, matplotlib, pytest).
- **Tests grün halten** — der Stand ist „159 grüne Tests, alles auf main". Vor
  einem Commit `pytest` laufen lassen; jede neue Quelle bringt ihren Test mit.
- **Doku-Sprache:** Status-Docs auf Deutsch, Code/Pfade/README auf Englisch.
- **Versionskontrolle:** GitHub-Remote `esoxconsultwerners/esox-master`.

---

## Roadmap (Kurzform)

- **Jetzt–November:** Katalog nutzen (Phi-Lab-Antrag mit den ehrlichen Zahlen,
  erste Dossiers, Website-Abbildungen). Optional Full-Range-NIR-Matcher, wenn ein
  Showcase-Dossier ansteht.
- **DR4 (2. Dezember):** Density Catalog (Stage 2), Perioden-/Familien-Refresh,
  Domänen-Adaptions-Forwardmodell (der C2.5-Befund ist seine Rechtfertigung).
- **Danach:** Source Attribution (Stage 3, braucht vollständige MetBull),
  Thermal Inertia (Stage 4), LSST-Enrichment-Feed.
