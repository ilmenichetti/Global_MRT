# Manuscript TODO

Last updated: 2026-07-03. Covers everything outstanding **except** Bucket 2
(reviewer figure/caption tweaks), which is being handled this session. Update as
items close.

---

## 1. Reviewer comments needing data rework (Bucket 3)

**To tackle next week (week of 2026-07-06).**

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

- [ ] **Carbon input = belowground NPP only.** Aleksi (twice) and Shoji both
  question using the belowground fraction of NPP as the soil input term. Open
  points: aboveground litter is ignored; litter *quality* / chemistry is not
  represented; the implicit assumption that all belowground NPP becomes litter.
  Action: write a defensible justification paragraph in Methods, and/or revise
  the input term. (manuscript.tex, inline comment near BNPP definition.)
- [ ] **Harvest removals.** In managed forests (e.g. Finland, Sweden) a large
  fraction of NPP is exported to mills and never enters the soil. Aleksi suggests
  a managed-forest layer (Hansen Global Forest Change,
  https://earthenginepartners.appspot.com/science-2013-global-forest) and a
  fraction-of-NPP-removed term, possibly biome-specific. Requires new input layer
  and a pipeline re-run.
- [ ] **Peatland exclusion.** Peat observations may distort interpretation.
  Filter peat sites from the training data (global peat maps exist); the maps
  themselves can stay as is. Requires re-running the pipeline on filtered obs.

## 2. Manuscript placeholders to fill

- [ ] **Conclusions section** is a stub (four bullet drafts exist) and needs to
  be written into prose.
- [ ] **Zonality map (Fig.)** flagged "draft figure --- to refine together";
  finalise design.
- [ ] **Authors / affiliations:** author list slot 3, affiliations 2 and 3.
- [ ] **Data / code availability DOIs** and repository link.
- [ ] **Acknowledgements** and **Author contributions (CRediT)**.
- [ ] **Figure in-panel labels still say "MRT".** After standardising the text on
  transit time (tau), figure titles/axes baked into the PNGs are now inconsistent
  (e.g. step_18 ALE panels: "Global Mean MRT", "ALE effect on MRT (years)"; likely
  the map figures too). Sweep the plotting scripts' labels to transit time / tau
  and regenerate. (Deferred to avoid repeated heavy re-runs; batch with the next
  pipeline pass.)

## 3. Missing citations (\PH{[cite: ...]})

Grouped by location:

- [ ] Intro, classical models: Jenkinson; Parton/CENTURY; RothC lineage.
- [ ] Intro, ESM inheritance: Todd-Brown 2013; Wieder 2013.
- [ ] Intro, persistence-as-ecosystem-property: Kleber 2007/2021; Lehmann &
  Kleber 2015; Kallenbach 2016; Cotrufo (MEMS); Averill 2014.
- [ ] Discussion, edaphic mechanisms: Six 2002; Hassink 1997; Rasmussen 2018;
  Kleber 2007.
- [ ] Discussion, pH / microbial: Rousk 2010; Malik 2018; Fierer 2006.
- [ ] Discussion, geochemical (Ca-OM): Rowley 2018; Rasmussen 2018; Rothe 2019.
- [ ] Discussion, mycorrhizal / high-latitude: Averill 2014; Clemmensen 2013;
  Tedersoo 2014; Steidinger 2019.

## 4. Decisions / notes

- [ ] **Title.** Body now standardised on mean transit time (tau); title still
  reads "mean residence time." Lorenzo is fine with this for now; revisit only if
  a reviewer flags the mismatch.

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
