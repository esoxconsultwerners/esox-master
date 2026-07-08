# C1 - Degeneration analysis: findings & decision memo

*Esox Master Catalog, Phase C / step 1. How much compositional and taxonomic
information is reconstructible from 16 Gaia bands (374-1034 nm) + SDSS u/g?
This memo synthesizes stages C1.0-C1.4 and is written to be liftable into two
paper introductions. All separability is measured with GroupKFold by meteorite
(no sample leakage), realistic per-band Gaia noise, identical 550 nm
normalization, viability filtering (>=8 spectra & >=3 meteorites), group_kind
separation, and a documented Gaussian-bandpass resampling to the Gaia grid.*

---

## Top line

- **Composition:** a **catalog paper is supported for a well-defined separable
  subset** (74% of 23 viable primary meteorite groups stay pairwise-separable
  at AUC >= 0.75 under Gaia noise), **on the explicit condition that the
  degenerate groups are reported as degenerate**. It is a catalog with an
  honest confidence/degeneracy map - which is simultaneously the methods
  contribution. It is *not* "every asteroid gets a confident, unique analog".
- **Taxonomy / C4:** reliable to **taxonomic-complex level** (C/S/X + common
  end-members), **not** to full Bus-DeMeo class. Independent validation against
  PDS Bus-DeMeo agrees at **86%**. Achievable confident coverage of the 19,190
  core is **~57%** at complex level - the 3.9% -> ~60% jump is credible at
  complex level, and only there.

---

## 1. Composition (RELAB meteorite groups at Gaia resolution)

**Verdict: catalog-viable for the separable subset; degeneracy map is the
methods core.**

- 23 viable primary groups; **mean pairwise ROC-AUC 0.912 noise-free
  (information ceiling) -> 0.796 with realistic Gaia noise**.
- **17 / 23 groups (74%)** remain separable (AUC >= 0.75 vs a majority of other
  groups) with noise - well above the 30% kill threshold.
- **Separate strongly** (the confident-analog set): the HED achondrites
  (diogenite 0.94, howardite 0.91, eucrite 0.88), primitive achondrites
  (lodranite-acapulcoite 0.90), and carbonaceous CM 0.90 / CO 0.86.
- **Collapse at 16 bands** (must be flagged, never asserted): the stony-irons
  (mesosiderite 0.44, pallasite 0.59) and CI 0.49; and critically the
  **ordinary-chondrite subtypes H/L/LL blend into each other** (each ~0.72-0.74
  one-vs-rest) - they are separable *from other classes* but not *from one
  another*.
- **SDSS u/g blue extension** lifts mean pairwise AUC by **+0.029** - a real,
  quantified advantage over Gaia-only work (largest help for the carbonaceous
  types where the blue drop-off discriminates).

*Single most important limitation:* the ordinary chondrites - by far the most
common meteorite falls (~80% of falls, per A8) - are the least internally
separable at Gaia resolution. A composition catalog can confidently say
"ordinary chondrite" but usually not "LL vs L vs H". Any framing must own this.

## 2. Space weathering (nuisance parameter)

- Reddening model W(lambda) = exp(-C / lambda), anchored at 550 nm, strength C
  swept and **marginalized** into the classifier.
- **Mean separability is robust**: marginalizing weathering changes mean
  pairwise AUC by **+0.030** (augmentation compensates). The cost is
  concentrated in specific weathering-ambiguous pairs - worst **CO vs EL
  -0.225**, CR vs R -0.183, CK vs aubrite -0.129.
- A **per-object weathering-strength column** is derivable (reddening slope C,
  median 0.109) - itself a catalog science product.

*Single most important limitation:* the reddening model is a placeholder form;
its literature reference is deliberately left as a TODO for verification before
publication (no fabricated citation). The named weathering-ambiguous pairs are
the ones a catalog must flag rather than assert.

## 3. Taxonomy / C4 (real GASP spectra)

**Verdict: complex-level reliable; class-level over-claims.**

- Trained on **358 GASP-Mahlke** core labels; validated on **36 held-out PDS
  Bus-DeMeo** core objects.
- **Complex level** (S/C/X/V/D/A/K-L): balanced accuracy 0.546
  (rare-class-penalized), but **86% agreement with independent Bus-DeMeo** on
  real-distribution data - the common complexes (S/C/X) are reliable; the rare
  end-members (A, D, K/L) are the weak spots.
- **Full Bus-DeMeo class level**: balanced accuracy **0.416** - not supported.
- Feature ablation (subset-matched): SDSS u/g +0.006, albedo +0.005, phase
  slope -0.042 - the 16 bands carry almost all the taxonomic signal at complex
  level.
- **Achievable coverage: ~57%** of the 19,190 core get a confident complex
  label (max prob >= 0.5). This is the honest basis for the 3.9% -> ~60%
  coverage jump - credible **at complex level only**.

*Single most important limitation:* only 358 labeled training objects, heavily
S-dominated; the rare complexes are label-starved, so C4 must report per-object
confidence and not present a uniform class-level catalog.

## 4. Apophis external-spectrum demo (showcase mechanism)

- Apophis has no Gaia spectrum; served from MITHNEOS NIR spectra resampled to
  the 7 Gaia bands >= 770 nm (normalized at 770 nm).
- The matcher assigns **54% probability mass to ordinary chondrites** (top
  analog L 0.23, LL 0.18, H 0.13) - **correctly recovering an ordinary-chondrite
  analog**, consistent with the literature Sq classification and validating the
  external-spectrum mechanism.
- But it **cannot resolve LL specifically** vs L/H - the exact H/L/LL
  degeneracy measured in section 1, now demonstrated on the showcase object.

*Single most important limitation:* NIR-only, 7 bands, no visible slope and no
2 um pyroxene band at Gaia resolution - a lower bound on what a full external
spectrum resolves. The mechanism works; the Gaia-grid resolution is the limit.

---

## Decision

Both papers are supported and honest:

- **Composition** -> a catalog paper that assigns confident meteorite analogs
  to the separable subset (achondrites, CM/CO carbonaceous) **and** publishes
  the pairwise degeneracy map as its methods spine, explicitly flagging the
  ordinary-chondrite and stony-iron degeneracies. If the emphasis is the
  degeneracy map itself, the identical evidence is a methods/limits paper. The
  data support either framing; they do **not** support "every asteroid gets a
  unique confident analog".
- **Taxonomy / C4** -> proceed, targeting **complex-level** classification with
  per-object confidence, ~57% credible coverage. Independent of the composition
  verdict.

**Next steps by verdict:** the composition matching pipeline is justified *only
for the separable subset* and is a separate later package (not built in C1).
C4 (complex-level taxonomy classifier) is greenlit.
