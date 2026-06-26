# Results reconciliation — 13c/13d/13e vs current manuscript

**Status:** analysis done; **manuscript prose NOT yet edited** (awaiting narrative
discussion). This memo maps every number that should change, old → new, with line
references, and separates *mechanical corrections* from *framing decisions*.

Production sources:
- `outputs/13c_decomposition.rds` — MAIN (fungal+SPUN, n=154,675), 500 trees
- `outputs/13c_decomposition_allbio.rds` — sensitivity (all-9, n=101,400)
- `outputs/13_model_comparison.csv` — the *current* table source (nested M1–M7)
- `outputs/13d_*` — biological parsimony (appendix; run pending/just-finished)
- `outputs/13e_*` — Tier-1 Moran's I (new content; run pending/just-finished)

---

## A. Mechanical correction — the "10 percentage points" is a footprint artifact

The current text (L826–827, comment L998–999) says climate explains 28 % and a
"further **10 percentage points**" is non-climatic. That increment compares
**different samples**: the manuscript even states this (L645–649):

| Model | n_obs | R²_CV |
|---|---|---|
| M1 Climate | 178,332 | 0.280 |
| M7 Full | 101,400 | 0.379 |

0.379 − 0.280 = 0.099 ≈ 10 pp, but climate is scored on n≈178k (a larger, *less
predictable* sample) and the full model on n≈101k. On a **common row set** the
climate baseline rises and the increment shrinks:

| Footprint | Climate R² | Full R² | Beyond-climate |
|---|---|---|---|
| MAIN common (n=154,675) | 0.307 | 0.380 | **+0.073 ≈ 7 pp** |
| all-9 common (n=101,400) | 0.307 | 0.378 | +0.071 ≈ 7 pp |

**Recommended edit (safe regardless of framing):** change "10 percentage points"
→ "**~7 percentage points**" (L827, L998–999), and add one clause that the
increment is measured on a *common footprint* so it is not inflated by the
M1-on-a-larger-sample artifact. Full R² is essentially unchanged (0.379 → 0.380),
so the headline performance number stands.

---

## B. Framing decision (FOR DISCUSSION) — nested M1–M7 → common-row Shapley

Current Results (L626–643) report **order-dependent nested increments** on
**footprint-varying** samples. 13c replaces this with one common-row set and a
symmetric Shapley split. Qualitative story is **unchanged**; numbers and rigour improve.

**13c MAIN decomposition (n=154,675, full R²=0.380):**

| Domain | Single R² | Shapley share | Shapley R² | Unique |
|---|---|---|---|---|
| Climate | 0.307 | 35.1 % | 0.134 | +0.036 |
| Edaphic | 0.307 | 32.8 % | 0.125 | +0.022 |
| Biological | 0.227 | 22.3 % | 0.085 | +0.010 |
| LandUse | 0.110 | 9.7 % | 0.037 | −0.001 |

Σ singles = 0.950 vs full = 0.380 → **overlap 0.57** ("correlated, not redundant").

Decision needed: **(i)** lead Results with Shapley (singles-first opening), using
the nested ladder only to *explain* order-dependence; **(ii)** whether to keep an
M1–M7 table, replace it with a coalition/Shapley table, or show both.

**Note the framing shift for biology:** nested-last biology is only +0.010, but
its Shapley share is 22.3 %. We must lead with Shapley or biology looks spuriously
negligible. (Land-use is the opposite: genuinely redundant — unique ≈ 0.)

---

## C. Number map (old → new), by line

| Line(s) | Current | New (13c MAIN) | Type |
|---|---|---|---|
| 626–628 | climate 28.0 %, R²_CV 0.280 | climate single R² 0.307 (common rows) | framing/B |
| 631 | edaphic ΔR² +0.076 (M2−M1) | edaphic single 0.307; Shapley 32.8 % | framing/B |
| 632 | biological ΔR² +0.057 (M4) | biological single 0.227; Shapley 22.3 % | framing/B |
| 633 | land-use ΔR² +0.030 (M3) | land-use single 0.110; Shapley 9.7 % | framing/B |
| 634–637 | M6/M7 R²_CV 0.379 | full R² 0.380 (n=154,675) | minor |
| 640–643 | M7 OOB 0.596 vs CV 0.379 | unchanged (13c has no OOB; keep from step 13) | keep |
| 711–712 | LandUse unique ≈ 0 (M6 vs M7) | confirmed: LandUse unique = −0.001, Shapley 9.7 % | confirm |
| 826–827 | climate 28 %, +10 pp | climate ~⅓ of attributable; +~7 pp beyond | **A (fix)** |
| 833–835 | "Edaphic largest non-climate" | HOLDS (Shapley 32.8 % > bio 22.3 % > LU 9.7 %) | confirm |
| 866 | "Biological second largest ΔR²" | HOLDS among non-climate (2nd: 22.3 %) | confirm |
| 878–880 | M6 vs M7 ΔR² < 0.001 | HOLDS (LandUse unique ≈ 0) | confirm |
| 998–999 | comment "~10 pp" | "~7 pp" | A (fix) |
| Table L1177–1184 | M1–M7 nested table | regenerate or replace w/ coalition+Shapley | framing/B |

---

## D. New appendix content (additive — not reconciliations)

- **Biological parsimony (13d):** per-layer coverage (fungal 99.8 %, SPUN 96.2 %,
  Barceló 62 %), collinearity, single-layer skill-vs-coverage, subset comparison.
  Justifies MAIN = fungal+SPUN (drop Barceló-3, unique ΔR² ≈ +0.001). Figs:
  `plots/appendix/biological_{collinearity,single_layer_skill}.png`. Slots into the
  existing appendix TODO block (~L1365). *(production numbers from 13d_run.log)*
- **Tier-1 Moran's I (13e):** beyond-climate residual Moran's I ≈ 0.25 at 2° and 5°
  (p ≪ 0.001) → spatially organised, not noise. Complements the weak/scale-dependent
  ICC. Figs: `plots/step_13c_commonality/13e_{block_residual_map,moran_scatter}.png`.
  New short Results/appendix paragraph for Tier-1 rigour. *(numbers from 13e_run.log)*
- **Provenance table (L1349–1353):** still the 300-tree numbers; separate
  `provenance_control.R` (not 13c). Regenerate only if we want 500-tree consistency.

---

## E. Open block-size sensitivity (proposed, not yet run)

Reviewer-proofing: re-run the 13c Shapley split at 1°/2°/5° CV blocks to show the
shares are stable (kills "why 2°?"). Cheap via a `BLOCK_SIZE` sweep in 13c.
