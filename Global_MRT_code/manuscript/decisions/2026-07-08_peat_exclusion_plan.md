# Peatland (+ open water): full exclusion over the pipeline (execution plan)

**Date:** 2026-07-08 (water exclusion added 2026-07-09) **Status:** Planned
(supersedes the flag-not-mask decision of 2026-07-07 once executed)

> ⚠ **TEXT UPDATE PENDING (Lorenzo, 2026-07-09).** Once the exclusion re-run lands,
> the MANUSCRIPT TEXT must be swept to match: every quoted number (R^2, beyond-climate
> gap, Shapley shares, Moran's I, n=, ESM ratios), the peat paragraph (flag-not-mask ->
> excluded), a new sentence on open-water masking, figure captions, and the data/methods
> citation of GPM 2.0. This is Phase 5 below — do NOT consider the work done until the
> text is reconciled. Flagged here so it is not forgotten between sessions.

**Phase 0 status (2026-07-09): DONE.** GPM 2.0 downloaded + inspected, licence checked,
D6 (water) added. Ready for D1-D6 lock then Phase 1.

**Phase 1 status (2026-07-09): DONE.** New `08b_extract_peat.R` written + run on branch
`peat-water-exclusion`. Per-obs flag (6,922 obs = 3.49% flagged peat; GPM catches ~2.7x
more than modal-Histosol) -> `outputs/08b_peat_flag.rds`. 0.1deg peat/water fraction +
combined exclusion mask -> `outputs/MRT_predictions/{peat_frac,water_frac,exclusion_mask}_0p1deg.tif`.
SANITY GATE PASSED: all 5 benchmark peatlands render correctly (Congo/W.Siberia/HudsonBay/
Fennoscandia/SE-Asia), overlay at `plots/step_08b_peat/08b_peat_sanity_overlay.png`. NOTE:
total flagged peat area 9.04 M km^2 (tier1 6.92 + tier2 2.12); tier2 matches the quoted
2.08, tier1 above the web-quoted 3.72 (that ref was a web summary, not the readme) -- GMC
uses an inclusive definition; over-inclusive is the safe direction for exclusion. Map mask
= 40,551 land cells.

**Phase 2 status (2026-07-09): DONE.** Peat exclusion is now the MAIN DEFAULT (Option A,
Lorenzo's call): `MRT_13C_PEAT=drop` + `MRT_PEAT_SOURCE=gpm` are the defaults, writing the
UNSUFFIXED headline files so maps + numbers follow with no repointing. Shared filter in new
`Global_MRT_code/peat_filter.R` (`apply_peat_filter()`), sourced by the 4 fitting scripts
that read 12b directly: 13, 13c, 13e, 13g (13h chains off 13c; 14/18 off the 13 model).
Suffix guards: histosol cross-check -> `_nopeat`, peat-kept -> `_peatkept` (neither can
clobber the headline). Smoke-tested: gpm drops 6922, histosol drops 2573 (=old count), keep
returns all. Reversibility = git (feature branch); no `_peatkept` re-run needed unless asked.
NEXT = Phase 3 heavy refit (13, 13c 500-tree Shapley overnight, 13e, 13g, 13h) then NUMBERS GATE.

## Why we are switching

The flag-not-mask decision was correct on the *numbers* (dropping Histosols barely
moves anything), but Lorenzo is making the **framing call**: scope the paper
explicitly to mineral soils rather than carry peat as a flagged caveat. This is the
"Revisit if" clause of `2026-07-07_peat_flag_not_mask.md`. It is a narrative /
elegance choice, not a numbers-driven one. It also upgrades the filter quality
(a real peat map catches thin/degraded peat the 0.1-degree modal Histosol class misses),
so it is a strictly better answer to Aleksi, not only cosmetic.

Expected effect (known from the sensitivity): beyond-climate gap 7.3 -> ~6.7 pp,
climate-only R^2 up, Moran's I ~flat, Shapley shares stable to <1 pp. Full exclusion
moves the flagship number slightly in the CONSERVATIVE direction.

## Source (decided)

**Global Peatland Map 2.0 (GPM 2.0)**, Greifswald Mire Centre / Global Peatlands
Initiative, released 2024-12. Two tiers: "peat-dominated" and "peat in soil mosaic".
Native 0.5 km, freely downloadable at 1 km. Single authoritative product; no second
map needed (its two tiers give the robustness handle internally).

## Decisions to lock BEFORE any code (5 min each)

- D1 Tier(s) to exclude for the fit filter: default = BOTH tiers (peat-dominated +
  mosaic). Report peat-dominated-only as the robustness variant if a reviewer pushes.
- D2 0.1-degree map-mask threshold: mask a cell if peat fraction >= X. Default X = 0.5
  (majority-peat cell). Tie to D1.
- D3 Observation filter rule: drop an observation if GPM at its coordinate is peat
  (tier per D1). Point sampling of the GPM raster at obs lon/lat.
- D4 Map visual for masked pixels: default = grey ("not interpreted"), caption note.
  Alternative = hatch. Cosmetic, decide once.
- D5 Keep the old modal-Histosol filter as a secondary cross-check? Cheap (already
  built). Default = yes, mention agreement in one sentence.
- D6 (added 2026-07-09) Exclude OPEN WATER from the map product too, for completeness.
  Source = REUSE `landcover/landcover_fractions_0p1deg.tif` **band 8 "water"** (already
  at 0.1 deg, no download). D6b threshold: mask cell if water fraction >= 0.5 (same cutoff
  as D2). D6c: merged into the same grey "not interpreted" mask as peat, caption WORDED so
  the two are not conflated (peat = excluded-scope; water = no soil exists). MAP-ONLY: no
  soil obs sit on open water, so this does NOT touch the fit/attribution. Note: large water
  is likely already blank via SoilGrids NA -> this makes it explicit + cleans edge cells.

## Phases

### Phase 0 - acquire + inspect GPM 2.0  (Lorenzo + me; ~30 min)
- Lorenzo: download GPM 2.0, check licence for our use, drop in `Datasets/peat/`.
- Me: inspect CRS, resolution, tier encoding, extent; confirm two-tier coding.

### Phase 1 - new peat-layer script `08b_extract_peat.R`  (me; ~half day)
- Reproject GPM to WGS84 if needed.
- Aggregate to the 0.1-degree map grid -> peat fraction per cell (+ tier); write
  `Datasets/peat/peat_frac_0p1deg.tif` and thresholded `peat_mask_0p1deg.tif` (D2).
- Point-extract GPM at each observation -> `peat_flag` per site.
- SANITY GATE: overlay against known peat regions (Congo Cuvette Centrale, West
  Siberian Lowlands, Hudson Bay Lowlands, Fennoscandia, SE Asian coastal peat).
  Do not proceed until the layer looks right.

### Phase 2 - wire into data assembly  (me; ~1 hr)
- Join `peat_flag` into the model-ready table (step 12 / 12b assembly). Keep both the
  GPM flag and the existing `soilclass_wrb_name` Histosol name (D5).
- Extend the filter switch: `MRT_13C_PEAT=drop` repointed to the GPM flag as MAIN
  (add `MRT_PEAT_SOURCE=gpm|histosol` so both are reproducible).

### Phase 3 - re-run production fit + attribution, peat dropped as MAIN  (me; background compute, ~overnight)
- Re-run: 13 (RF fit), 13c (500-tree Shapley; long pole, hours), 13e (Moran),
  13g (provenance), 13h (block-size). Capture new headline numbers.
- NUMBERS GATE: confirm direction (gap down ~0.6 pp, clim R^2 up, Moran ~flat). If
  anything moves materially more than the sensitivity predicted, STOP and investigate
  (the GPM filter is broader than modal-Histosol, so a somewhat larger move is
  possible and fine; a sign flip or big jump is not).

### Phase 4 - maps re-run + mask  (me; ~half day + compute)
- 14 extrapolation (or reuse prediction) then apply the 0.1-degree peat mask -> peat
  pixels set NA (D4).
- Re-render every map with the mask: 15b zonality (symmetric / nested / dominant),
  17 global tau map, 20 latitudinal profile. 18 ALE re-run on the refit model (heavy,
  cached model reload).
- 20 ESM comparison: apply the SAME mask to both RF and ESM fields for a fair
  comparison; capture the shifted ratios.

### Phase 5 - manuscript sweep  (me + Lorenzo; ~half day)
- I do a full grep-and-update of every quoted number: R^2, beyond-climate gap,
  Shapley shares, Moran's I, n= (154,675 -> new), ESM ratios in Discussion, and all
  figure captions citing numbers + the tables.
- Rewrite the peat paragraph flag-not-mask -> excluded (Methods sec:robust); convert
  `tab:peat_sensitivity` to an "effect of exclusion" framing or retire it.
- Resolve the inline OPEN QUESTION comment block (peat) -> decided/executed.
- Retire `2026-07-07_peat_flag_not_mask.md` (mark Superseded) and finalise THIS file
  as the adopted record.
- COMPILE GATE: pdflatex clean, no undefined refs, table/figure counts balance, no
  orphaned "flag" / "Histosol-dominated areas" language left behind.

### Phase 6 - commit  (Lorenzo approves)
- Commit locally on a branch (see below). No push, no AI co-author trailer.
- Lorenzo syncs Overleaf + pushes; co-author re-review round (all numbers moved).

## Safety / rollback
- Do the whole thing on a **feature branch** (we are on `main`, which currently holds
  the committed flag-not-mask state). Keeps main intact until we are happy; trivial to
  abandon if the exclusion proves worse than expected.
- Existing `_nopeat` sensitivity outputs are preserved and act as a cross-check that
  the new GPM-based main run lands near the old modal-Histosol sensitivity.

## Two-day shape
- Day 1: Phases 0-3 (acquire, build layer, wire in, launch fit re-runs; Shapley runs
  overnight unattended).
- Day 2: Phases 4-6 (maps + mask, manuscript sweep, commit). Compute overlaps.
- Lorenzo hands-on time: the D1-D5 decisions, the two sanity gates, and reading the
  reworked text. A few hours total, spread out.

## Related
- Supersedes [[2026-07-07_peat_flag_not_mask]].
- Input-side circularity reasoning (why land-cover-structured filters are safe):
  see the input-term decision record / memory.
