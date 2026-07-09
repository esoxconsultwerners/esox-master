# Esox Master Catalog — Projektstand nach Composition-Arbeit
**Konsolidierung | Stand: 8. Juli 2026**
**Repo:** github.com/esoxconsultwerners/esox-master · 159 grüne Tests · alles auf main

Dieses Dokument ersetzt den Stand vom Vormittag. Es hält die *ehrliche* Composition-Realität fest, nachdem die C2.5/C2.6-Arbeit einen fundamentalen Befund zutage gefördert hat.

---

## Der Katalog heute

Zwei Tiers: **1.552.834 Objekte** (Bahn-Backbone) und ein dichter **19.190-Objekt-Spektralkern**. Der Kern trägt:

| Ebene | Deckung (Kern) | Status |
|---|---|---|
| Bahn + H | 100 % | fest |
| Gaia-Spektrum (16 Band) | 100 % | fest |
| Phasenkurve | 95,4 % | fest |
| Durchmesser + Albedo | 69,2 % | fest (Usui-Präzedenz) |
| Familie | 39,9 % | fest |
| Rotationsperiode | 33,0 % | fest |
| Literatur-Taxonomie | 3,9 % | fest |
| **Modell-Taxonomie (Komplex)** | **~57 % konfident** | **das Composition-Produkt** |
| Meteoriten-Analog (Gaia-Kern) | **0,08 % konfident, 81,9 % ehrlich unresolved** | Indikator, nicht Produkt |
| Meteoriten-Analog (NIR-Pfad) | 5 konfident (alle HED), 93,7 % in-manifold | verteidigbar, resolution-limited |

**Vollständige physikalische Profile: 1.258** (6,6 % des Kerns; von 204 durch C4).

---

## Der zentrale wissenschaftliche Befund (C2.5/C2.6)

Die ursprüngliche Annahme — jeder Asteroid bekommt einen Meteoriten-Analog aus seinem Gaia-Spektrum — ist **empirisch widerlegt**, und das ist ein wertvolles Resultat, kein Rückschlag:

- **Gaia-DR3-Asteroidenspektren und RELAB-Laborspektren sind strukturell divergente Verteilungen.** 95 % des Gaia-Kerns liegen außerhalb der RELAB-Mannigfaltigkeit (Median-Mahalanobis d² 188,8 vs. RELAB-Selbst-q95 23,1). Der ursprüngliche RandomForest-Matcher hat 2.208 konfidente „CM"-Analoge *fabriziert*, indem er weit außerhalb seines Trainingsbereichs extrapolierte.
- **Der NIR-Pfad ist dagegen domänen-gültig:** 93,7 % der MITHNEOS-Ground-Truth-Spektren liegen *innerhalb* der Mannigfaltigkeit (Median d² 1,5 vs. q95 13,4). Labor-NIR und beobachtetes NIR teilen dieselben Absorptionsbanden.
- **Konsequenz:** Konfidente Gaia-only-Meteoriten-Analoge sind ohne Domänen-Adaptionsmodell nicht ehrlich stützbar. Das Composition-Produkt ist die **C4-Komplex-Taxonomie** (auf echten Gaia-Spektren trainiert und mit 81 % gegen unabhängiges PDS validiert), nicht der Analog-Layer.

**Warum das ein Gewinn ist:** Ein Katalog, der nur die Dinge konfident nennt, die stimmen, ist einem vertrauenswürdig, der große, aber fabrizierte Zahlen liefert. Die konfidenten Analoge, die überleben, sind astrophysikalisch makellos (Gaia-Kern: die 15 sind fast alle V→HED/Vesta; NIR: 5 HED). Wenn diese Pipeline konfident ist, hat sie recht.

---

## Apophis — Showcase-Status, ehrlich

- In-manifold (d² 0,54 — besteht ein Gate, das 95 % des Gaia-Kerns ablehnt), aber Konfidenz 0,18 → **„indicative ordinary chondrite, unresolved"**. Die Methode darf angewandt werden; die Auflösung im auf 1034 nm beschnittenen NIR-Bereich reicht nicht, OC von Ureilit/CV sicher zu trennen.
- Vollständiges Dossier (docs/apophis_status.md): Aten-PHA-Bahn, NPA-Tumbler aus zwei unabhängigen Quellen, von den Risk-Lists nach dem 2004→2021-Bogen, erreichbar bei Δv 6,05 km/s, S-Komplex-Literatur, gemessene MITHNEOS-NIR-Spektren, HED-nahe LL-Laborbrücke.

---

## Zwei getrennte Follow-up-Projekte (nicht verwechseln)

Der Composition-Befund erzeugt zwei *verschiedene* künftige Arbeiten, die verschiedene Probleme lösen:

1. **Full-Range-NIR-Matcher** (kleiner, für das Showcase-Geschäft wertvoll): Externe NIR-Spektren auf dem *vollen* MITHNEOS-Bereich (bis 2,5 µm) gegen volle RELAB-Templates matchen, statt am Gaia-1034-nm-Gitter abgeschnitten. Hebt die NIR-Pfad-Konfidenz — nützt Apophis und jede Kundenanfrage mit gemessenem Spektrum. Eigenes Paket, wenig Aufwand.
2. **Domänen-Adaptions-Forwardmodell** (groß, Stage-1-Wissenschaft nach DR4): RELAB-Templates über ein Space-Weathering-/Phasen-/Instrument-Modell in den Gaia-Beobachtungsrahmen abbilden, *bevor* gematcht wird. Würde konfidente Gaia-only-Analoge *legitim* freischalten. Wochen-Arbeit, richtiger Moment ist die DR4-SSO-Pipeline. Der Mannigfaltigkeits-Befund ist seine wissenschaftliche Rechtfertigung.

---

## Offene kleine Werner-Aufgaben

1. **meteorite_group_properties.csv** — Zitationen noch alle TODO. Jetzt niedrigere Priorität (fast nichts bekommt konfidenten Analog), erst vor einem konkreten Kundendossier nötig.
2. **Gauß-Passband-To-do ist obsolet** — C2.5 hat gezeigt, dass der Versatz nicht passband-verursacht ist (edge-band-konzentriert). Die zwei Schrott-Randbänder (374, 1034 nm) werden im Matcher gedroppt. Kein Passband-Tausch nötig. (Der authoritative DR3-Passband-Satz liegt trotzdem checksummed im Repo, falls je gebraucht.)

---

## Roadmap

- **Jetzt–November:** Katalog nutzen. Phi-Lab-Antrag mit den ehrlichen Zahlen (1.258 Profile, 57 % Komplex-Taxonomie als Produkt, konfidente Analoge nur mannigfaltigkeitsgedeckt). Erste Dossiers (Apophis, ehrlich). Website (Abbildungen liegen: coverage_matrix, funnel v0/v1, calibration, nir_manifold). Optional: Full-Range-NIR-Matcher, wenn ein Showcase-Dossier ansteht.
- **2. Dezember (DR4):** Density Catalog (Stage 2), plus Perioden-/Familien-Refresh (die neuen Engpässe aus dem C4-Funnel). Natürlicher Moment für das Domänen-Adaptions-Forwardmodell.
- **Danach:** Source Attribution (Stage 3, braucht vollständige MetBull — die A8-Census ist auf ~2013 eingefroren, Vorbedingung), Thermal Inertia (Stage 4), LSST-Enrichment-Feed.

---

## Der ehrliche Satz zum Gesamtstand

Der Composition-Layer ist end-to-end ehrlich und konsistent: Gaia-Kern und NIR-Pfad folgen demselben Manifold-Gate, kein Matcher fabriziert Analoge außerhalb seiner Gültigkeit, jede konfidente Aussage ist mannigfaltigkeitsgedeckt. Das Composition-Produkt ist die C4-Komplex-Taxonomie plus konfidente Analoge, wo verteidigbar. Das ist ein Katalog, den du einem Phi-Lab-Gutachter oder einem Planetary-Defense-Kunden zeigen kannst, ohne eine einzige Aussage zurücknehmen zu müssen — und das ist mehr wert als jede große, aber fabrizierte Zahl.
