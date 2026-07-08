# C2.1 - Composition matcher (separable subset)

Analog target classes (21): ['CI', 'CK', 'CM', 'CO', 'CR', 'CV', 'EH', 'EL', 'R', 'aubrite', 'brachinite', 'diogenite', 'eucrite', 'howardite', 'lodranite-acapulcoite', 'lunar-meteorite', 'martian', 'mesosiderite', 'ordinary_chondrite', 'pallasite', 'ureilite']. Separable subset from C1 pairwise AUC (>=0.75 vs majority); H/L/LL merged into **ordinary_chondrite** (C1: mutually ~0.72, never emit a confident OC subgroup). Degenerate sink classes (reported, never confidently asserted): ['CI', 'CK', 'CO', 'R', 'mesosiderite', 'pallasite'].

RandomForest on RELAB templates (16 Gaia bands, 550 nm norm). Self-classification GroupKFold-by-meteorite: balanced acc 0.332, reliability ECE 0.042 (c2_reliability.png).

## C4 complex prior

Complex -> plausible analog groups (implausible down-weighted x0.25, not zeroed): {'S': {'ureilite', 'lodranite-acapulcoite', 'brachinite', 'ordinary_chondrite'}, 'C': {'CK', 'CO', 'CR', 'CV', 'CM', 'CI'}, 'X': {'pallasite', 'aubrite', 'EH', 'mesosiderite', 'EL'}, 'V': {'diogenite', 'eucrite', 'howardite'}, 'A': {'pallasite', 'brachinite', 'ureilite'}, 'D': {'CR', 'CM', 'CI'}, 'K/L': {'CO', 'CK', 'CV'}}. Where taxon_esox is 'unclassified', a flat prior is used. **The prior reduces mean analog entropy by 8.3%** - the quantitative reason C4 came before C2.
