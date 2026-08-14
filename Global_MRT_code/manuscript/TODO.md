# Manuscript TODO

Last updated: 2026-07-09. Covers everything outstanding **except** Bucket 2
(reviewer figure/caption tweaks), which is being handled this session. Update as
items close.

---

## 0. Current revision (peat/water exclusion + Lorenzo comments, 2026-07-09)

DONE and pushed to `main` (peat-excluded results + comment round + Fig 3). Below closed.

- [x] **Full number sweep** to the peat-excluded headline (gap 7.31->7.12 pp;
  R^2 0.38->0.39, clim 0.31->0.32; Shapley 35.4/32.7/22.2/9.6; Moran 0.253/0.261;
  n 154,675->150,863; tau mean/median 23.5/23.0, range 3.9-72.8; ESM 2.45x vs 9.19x).
- [x] **Rewrite the peat text** flag-not-mask -> excluded (scoped to mineral soils);
  open-water masking note added; GPM 2.0 cited (Greifswald Mire Centre 2022,
  CC BY-NC-SA). NB: the "signal survives peat exclusion" robustness point stays OUT
  of the MS (coauthor reply only, Lorenzo's call).
- [x] **ESM comparison domain asymmetry** caveat added to the ESM Discussion (our tau
  field is coverage-limited vs the spatially complete ESMs, so part of the poleward-
  steepness difference may be a domain mismatch, not model behaviour alone).

---

## 1. Reviewer comments needing data rework (Bucket 3)

**RESOLVED 2026-07-07** without a pipeline rebuild. The input-term critique was
answered with a carbon-input sensitivity sweep (three conventions + a
literature-anchored harvest h-curve) showing the relative/zonality results are
invariant; peat handled flag-not-mask with a sensitivity table. See
`13k_input_sensitivity_appendix.R`, Fig.~\ref{fig:input_sensitivity},
Tables~\ref{tab:input_sensitivity}/\ref{tab:peat_sensitivity}, Methods
Sec.~\ref{sec:robust}, and `decisions/2026-07-07_peat_flag_not_mask.md`.

These came up independently from more than one co-author, so they are not optional.

**Framing (from 2026-07-03 discussion).** Steady-state tau can be written from the
input side (tau = Cs/input, what we do, input ~ BNPP) or the output side
(tau = Cs/Rh, respiration, what Varney et al. 2020 do); same quantity, different
observable, different error structure. Useful context for the limitation, but it
positions the choice rather than resolving the critique. Key point: the paper's
claim is about *relative/spatial* structure (zonality), so what matters is whether
an input error is spatially smooth (mostly rescales magnitude, fairly harmless) or
spatially structured (can masquerade as zonality). On that test the three comments
are NOT equal weight:
  - Aboveground litter + litter quality -> largely magnitude / steady-state
    caveats; probably a careful paragraph, not a re-run.
  - **Harvest removals -> the one with teeth**: concentrated in managed
    boreal/temperate forests, so it is spatially organised and overlaps the
    beyond-climate signal. This likely warrants an actual managed-forest
    sensitivity test, not just prose. (Open question to settle: do the
    literature BNPP fractions themselves carry enough biome-structured error to
    matter even for the "smooth" cases?)

- [x] **Carbon input = belowground NPP only.** Answered via the carbon-input
  sensitivity sweep (belowground / total NPP / total$-$harvest, same NPP field):
  relative results invariant, only LandUse share responds. Belowground defended as
  a conservative boundary (root-C preferential stabilisation) at the BNPP
  definition; justification + sensitivity in Methods Sec.~\ref{sec:robust}.
- [x] **Harvest removals.** Addressed with a crude land-cover harvest export on
  tree-dominant pixels, swept across the literature range of $h$ (0.2--0.6): the
  beyond-climate gap does not fall (Fig.~\ref{fig:input_sensitivity}D). A dedicated
  managed-forest layer (Hansen/Lesiv) was judged unnecessary; the effect lands on
  the LandUse axis and the crude sweep already brackets it.
- [x] **Peatland exclusion.** Flag-not-mask: dropping Histosol-classified obs
  leaves the gap and Moran's $I$ essentially unchanged (Table~\ref{tab:peat_sensitivity});
  records retained in the fit, flagged on the maps. Decision record filed.

## 2. Manuscript placeholders to fill

The live target is now `main_Nature.tex` (+ `main_Nature_SI.tex`);
`manuscript.tex` is the superseded long-format draft.

- [ ] **Front/back matter**, the only `\PH` left in `main_Nature.tex`:
  additional authors + affiliation 4, `[DOI]` (data), `[DOI/repository]` (code),
  further acknowledgements, CRediT.
- [ ] **Author list.** Candidates noted on the draft: Todd-Brown, Lim, Abramoff,
  Schwarz, Mäkipää. Bruni's affiliation still to confirm with her.
- [ ] **AI-use disclosure** required by Springer Nature (see target-journal note).
- [ ] **Prose/narrative polish** of the restored arc, plus Lorenzo's note asking
  the synthesis to state more openly that the framework avoids treating emergent
  properties reductionistically.
- [x] **Conclusions.** No longer a stub: the Nature format ends on the synthesis
  and modelling-consequence beats.
- [x] **Zonality map.** Finalised as the combined three-panel Fig. 3; the
  "draft figure" marker is gone.
- [x] **Figure in-panel labels.** MRT -> tau sweep done and re-rendered
  (17/18/13e/20, 15b C/D); nothing says "MRT" in a label now.

## 3. Missing citations (\PH{[cite: ...]}) --- DONE 2026-08-14

All 14 citation placeholders in `main_Nature.tex` are wired; `references.bib`
went 27 -> 55 entries, every new entry generated from Crossref rather than typed.
Three pre-existing entries had wrong metadata and were corrected (Schmidt2011
author list, Yu2022 title/authors, Varney2020 authors).

- [x] Intro, classical models: Jenkinson1990; Parton1987; Coleman1996 (RothC).
- [x] Intro, ESM inheritance: ToddBrown2013; Wieder2013.
- [x] Intro, persistence-as-ecosystem-property: Kleber2007; LehmannKleber2015.
- [x] Intro, pedological zonality: Dokuchaev1883; Jenny1941.
- [x] Methods, root-derived carbon: Rasse2005; Jackson2017.
- [x] Results, edaphic mechanisms: Hassink1997; Kleber2007; Rasmussen2018.
- [x] Results, pH / microbial: Rousk2010; Malik2018.
- [x] Results, mycorrhizal / high-latitude: Averill2014; Clemmensen2013;
  Frey2019; Tedersoo2014; Steidinger2019.
- [x] Discussion, forest-management intensity: Suvanto2025.

Open for Lorenzo's approval (my picks, not from the placeholders):
Davidson2006 + Carvalhais2014 (first-order climate control); Levin1998
(emergence resists reconstruction); StricklandRousk2010 + Frey2019 (higher
fungal proportion shortens tau --- the one that cuts against the conventional
"fungi = slower turnover" reading).

In the bib but uncited, deliberately: Batjes2024 (WoSIS --- must stay uncited,
it is not our data source), Hijmans2024 (terra; could go in Methods software),
Six2002, Rowley2018 (Ca-OM; no carbonate paragraph in the Nature format),
Fierer2007.

## 4. Decisions / notes

- [x] **Title.** Changed to "transit times" for consistency with the body.

---

### Addressed this session (for reference)
- Abstract reframed to lead with the zonality rediscovery; edaphic/biological
  defined inline; litter-quality nod; observation-based-benchmark clause.
- Terminology standardised on transit time tau; the 10 stray "MRT" replaced;
  intro sentence added justifying transit time (answers Shoji).
- Zonality-map caption "how to read it"; horizon-standardisation clarification;
  Limitations forward-look on unpredictable variance.
- Bucket 2 (complete): Moran per-Koppen fits + caption; ALE shared y-axis +
  caption; 16 km driver sentence; temperature-tau appendix scatter (new Step 22
  script, Fig. A5, cf. Varney et al. 2020 Fig 3a; apparent Q10 ~ 1.35).
- Varney et al. 2020 (spatial emergent constraint, tau = Cs/Rh, effective
  Q10 ~ 2.5) integrated: bib entry added; cited in intro (sets up the
  "is the spatial tau-T gradient really thermal?" question), in the ESM
  discussion (our R2=0.06 tau-T fit + Moran's I=0.25 qualify their space-for-time
  premise; note their respiration-based tau vs our input-based tau), and the
  Fig. A5 caption.
