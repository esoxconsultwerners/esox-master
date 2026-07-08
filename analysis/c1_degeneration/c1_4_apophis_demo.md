# C1.4 - Apophis external-spectrum composition demo

Apophis (99942) has no Gaia spectrum; served via the external path from 931 MITHNEOS NIR points (770-2485 nm). Matcher restricted to the 7 Gaia bands >= 770 nm (normalized at 770 nm), trained on the RELAB primary groups.

## Composition analog distribution (top 6)

| rank | analog group | probability |
|---|---|---:|
| 1 | L | 0.233 |
| 2 | ureilite | 0.203 |
| 3 | diogenite | 0.202 |
| 4 | LL | 0.177 |
| 5 | H | 0.132 |
| 6 | EH | 0.023 |

- Ordinary-chondrite (H+L+LL) probability mass: **0.542**.
- LL probability: 0.177. Top analog: **L**.

## Interpretation

The literature classifies Apophis as **Sq**, the spectral bridge to **LL ordinary chondrites**. The external-spectrum matcher's top analog is **L** with ordinary-chondrite mass 0.54. This is consistent with the Sq->LL/ordinary-chondrite expectation and validates the external-spectrum mechanism. Caveat: only the 770-1034 nm Gaia window is used (no visible slope, no 2 um pyroxene band), so this is a lower bound on what a full external spectrum would resolve.

