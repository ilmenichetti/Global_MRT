# Storyline outline — zonality of soil-carbon turnover (working doc)

*Operational bones for the rewrite. Numbers are the production (500-tree, MAIN
6-layer biology) values unless noted. Companion: `results_reconciliation_2026-06-26.md`
(number map + figure dispositions); every figure below is regenerable from the
listed script.*

---

## Thesis (one sentence)

Beyond climate, soil-carbon turnover carries a **coherent, attributable, mappable
spatial organisation** — modest in size but real — that Earth System Models
**over-represent and disagree on**; we offer it as an interpretive framework and an
observational benchmark.

## OCAR spine

- **O (opening):** Turnover is treated as climate-driven, the rest as noise. But on a
  *spatially-honest* (block-CV) basis, climate explains only **~⅓** of the predictable
  variance.
- **C (challenge):** Is the other **~⅔** just noise — or is it *organised*,
  *attributable to processes*, and *mappable*? And do ESMs get that organisation right?
- **A (action):** Common-row Shapley attribution + a spatial-organisation test
  (Moran's I) + the zonality map + mechanistic ALE curves, all wrapped in an honesty
  apparatus (spatial block-CV, noise floor, provenance control, biological parsimony,
  block-size invariance).
- **R (resolution):** Non-climate context imposes a coherent, ecologically-interpretable
  zonality — **edaphic-led, with biology a distinct second axis** — modest in magnitude
  but spatially organised and ecological (not a sampling artifact). ESMs amplify it
  **~3.6× too steeply** and disagree among themselves by an order of magnitude.

## The four pillars (what carries the paper)

| Pillar | Question | Headline result | Script |
|---|---|---|---|
| Shapley attribution | *how much?* | Climate 35 / Edaphic 33 / Biology 22 / LandUse 10 % | 13c |
| Zonality map | *where?* | non-climate bends τ coherently (Amazon/Congo longer; drylands shorter) | 15b |
| ALE response curves | *how?* | τ ↑ w/ seasonality, snow, elevation; ↓ w/ sand, bulk density | 18 |
| ESM contrast | *so what?* | obs polar/tropical 2.5× vs ESM ~9× (3.6× over-zonalization) | 20 / 20b |

## Key numbers (for the text)

- Spatially-honest full-model **R² = 0.38**; climate alone **0.31**; **beyond-climate ≈ 7 pp**
  (the old "~10 pp" was a footprint artifact — climate scored on a larger, less
  predictable sample).
- **Singles:** Climate 0.307 ≈ Edaphic 0.307 > Biology 0.227 ≫ LandUse 0.110;
  Σ singles 0.95 vs full 0.38 → **overlap 0.57** ("correlated, not redundant").
- **Shapley** (sums exactly to R²): Climate 35.1, Edaphic 32.8, Biology 22.3, LandUse 9.7 %.
  Robust to biological-set choice and to CV block size (shares move ≤ 2.4 pp over 1–5°).
- **Tier-1:** beyond-climate residual Moran's I = **0.25** (p ≪ 0.001) → spatially organised.
- **Noise floor:** variogram nugget/sill ≈ 0.56, range ≈ 16 km → ~⅔ unstructured noise;
  map only at meso (2°) scale.
- **Provenance:** only **2.4 %** of the residual is between source databases → ecological.
- **Biology:** a distinct axis; MAIN set = fungal + 5 SPUN (96 % coverage); the Barceló
  root layers dropped (inert: +0.001 unique R²).
- **Maps:** abiotic modulation near-global over soil land (~90 % between 60°S–50°N);
  real gap only in boreal high latitudes (SoilGrids-sparse — and where ESMs over-zonalize).
- **Abiotic vs biology, spatially:** correlate when each is taken over climate
  (r ≈ 0.4–0.5) but are ~orthogonal once biology is conditional on abiotic (r ≈ 0) —
  the overlap result, expressed spatially.
- **ESM over-zonalization:** polar/tropical τ — RF 2.5×, ESM median 9.2× (range 0.6–14×),
  ensemble/RF ≈ 3.1–3.6×; MPI even inverts the gradient.

## Figure plan

**Main (Results):**
1. **Zonality map (symmetric)** — `15b` `zonality_map_dual_symmetric.png`: (a) abiotic
   log(M5/M1), (b) biology log(M4/M1), (c) abiotic-vs-biology scatter. Near-global
   abiotic panel doubles as the headline. → *where*.
2. **Variance decomposition** — `13c` `13c_decomposition.png`: singles → build-up →
   Shapley. → *how much* + the overlap.
3. **ALE response curves** — `18` `18D_ALE_top10_MRT.png` (+ `18C_ALE_by_latband`). → *how*.
4. **ESM structural comparison** — `20` `fig_biome_relative_tau.png`. → *so what*.
   (+ one-line over-zonalization stat from `20b`.)

**Tier-1 (Results or Appendix):** `13e` block-mean residual map + Moran scatter.

**Appendix:** scale/noise-floor `13f`; provenance table `13g`; block-size `13h`;
biological parsimony `13d`; nested map + abiotic/biology correlation `15b`
`zonality_map_dual.png`; relative latitude heatmap `20`
`fig_latitude_heatmap_relative.png`; **absolute latitudinal profile `20`
`fig_latitudinal_profile.png`** (used by the Discussion side-arc below).

**Methods/Data:** `17` `latitudinal_sampling_coverage.png` (+ optionally the SOC
density heatmap).

**Dropped/obsolete:** step-17 M1–M7 prediction maps, all-models comparison,
histogram-by-latitude; step-20 absolute heatmap.

## Two deliberate narrative devices

- **Own the incomparability (ESM section):** state up front that absolute τ is not
  comparable across products (depth / pool / reference-period), so we compare *normalized
  spatial structure*. This pre-empts the "apples to pears" critique.
- **Resolved-tension side-arc (Discussion):** raise the absolute-magnitude caveat
  (Appendix abs. profile) → resolve it via the structural comparison. Small, for
  robustness; not central.

## Status / what's left

- **Settled:** Shapley-primary lead; tool-first framing; symmetric map main / nested
  appendix; correlation → brief Discussion; ESM structural reframe; step dispositions.
- **Open (small):** map min-fill coverage threshold; optional per-pixel Shapley map;
  final 500-tree refresh of the 300-tree block-size sweep.
- **Next:** the prose pass — write `manuscript.tex` section by section per the
  reconciliation memo (then a code-polish pass after the story is frozen).
