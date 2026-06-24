# Response to Reviewers — Round #1

**Manuscript:** *palimpsestR: An R Package for the Identification and Probabilistic Decomposition of Archaeological Palimpsests*
**PCI Archaeo #1019** — Preprint DOI: 10.5281/zenodo.19969277
**Recommender:** Alessio Palmisano · **Reviewers:** Giacomo Bilotti, Lise Bellanger

> **STATUS: DRAFT.** Each item below pairs the reviewer's comment with our planned/implemented response.
> Tags: `[DONE]` change made · `[TODO]` change pending · `[CODE]` package change · `[PAPER]` manuscript change · `[ALREADY]` already present in the current package (v0.23.0) and only needs to be surfaced in the paper.
> Section/line references will be finalised against the revised manuscript before submission.

---

## Reply to the Recommender (A. Palmisano)

We thank the recommender and both reviewers for a constructive and encouraging assessment. We are glad that both reviewers found the package usable and the topic relevant, and that Reviewer 1 reproduced the analyses from both the paper and the vignette without issue. We have addressed every point below. The two central themes — (i) a fuller, better-contextualised methodological exposition and (ii) a clearer manuscript structure — have guided a substantial revision of the methods and introduction sections. We note that, since the submitted preprint (package v0.17.1), the statistical core has continued to develop (current v0.23.0) in directions that directly answer several of the reviewers' concerns (mixed-type categorical modelling of cultural class, per-find chronological-uncertainty propagation, a recording-resolution advisor); these were partly under-described in the submitted text and are now made explicit.

### Mandatory PCI sections
- **Funding** — `[TODO][PAPER]` A *Funding* section has been added. *(state actual funding or "no specific funding")*
- **Conflict of interest** — `[TODO][PAPER]` Corrected to the plural form ("The authors declare…") to match the three-author byline.
- **Data, script and code availability** — `[DONE]` Code and the `villa_romana` dataset are archived on Zenodo (concept DOI 10.5281/zenodo.19881542) and cited in the reference list (Cocca 2026). `[TODO]` We confirm the dataset + reproducibility script are deposited with a DOI and metadata (not GitHub only), per PCI's perenniality requirement.
- **Supplementary information** — `[TODO][PAPER]` The methodological vignette (`introduction.Rmd`) is deposited as Supplementary Information with its own DOI and referenced in a dedicated section and in the reference list.

---

## Reviewer 1 — Giacomo Bilotti

**R1.1 — Detailed explanation & literature contextualisation of SEF / SEI / Diagonal Gaussian mixture, with examples.**
`[TODO][PAPER]` We expanded the *Methodological framework* with (a) a conceptual introduction to the SEF, (b) a formal definition and interpretation of the SEI, and (c) a pedagogical presentation of the Gaussian-mixture/EM model, each illustrated with a small example and positioned against related approaches (model-based clustering, spatial point-pattern analysis, Harris-matrix recording). This also clarifies *why* palimpsestr is useful in practice.

**R1.2 — BIC, ICL, PDI are not explained in the text or in the function help pages; "typically chosen at the BIC minimum" needs a reference; emphasise combining statistics with expert archaeological judgement.**
`[DONE][PAPER][CODE]` Added the BIC definition (\(-2\log\hat L + d\log n\); Schwarz 1978), the ICL definition (\(\mathrm{ICL}=\mathrm{BIC}+2\sum_i H_i\); Biernacki et al. 2000), and a PDI definition at the model-selection passage, which now stresses weighing statistical criteria against archaeological plausibility (the case study retains K=4 over the BIC-favoured K=6). Help pages of `compare_k()` and `pdi()` document all three metrics and cite Schwarz (1978) and Biernacki et al. (2000). **Statistical-correctness fix surfaced while documenting ICL:** the implementation computed `icl = bic - 2*sum(entropy)`, which (with lower-is-better BIC) *rewarded* fuzzy partitions — the opposite of the criterion's purpose. Corrected to `bic + 2*sum(entropy)` (regression test added); this changes `compare_k()` ICL values only and leaves BIC, PDI, and the case study unaffected (ICL appears in no manuscript figure). We disclose this correction as part of the revision.

**R1.3 — Diagnostic procedures explained more clearly in `introduction.Rmd` could be expanded and incorporated into the paper.**
`[TODO][PAPER]` Key diagnostic explanations from the vignette have been folded into the methods/visualisation sections; the vignette is retained and cited as Supplementary Information (see above).

**R1.4 — `archaeo_sim()`: use `floor()` so `tmin`/`tmax` are integer years.**
`[TODO][CODE]` `archaeo_sim()` now floors simulated date bounds to integer years, matching archaeological convention. Test + NEWS updated.

**R1.5 — `harris_from_contexts()`: the default verticality rule conflates genuinely problematic cases with normal features (ubiquitous violations); same concern for multi-period fills (pits/ditches). Address methodologically and discuss.**
`[DONE][CODE][PAPER]` `harris_from_contexts()` gains an `exclude_contexts` argument (default `NULL`, behaviour unchanged): the listed contexts are exempted from the verticality (depth-rank) penalty — for the fills of cuts and multi-period contexts where depth does not track chronology. The help page now carries a *Verticality assumption* section documenting the rule and its failure modes, and the manuscript's *Limitations* note records that the rule is frequently violated and that an explicitly recorded Harris matrix (`read_harris()`) should be preferred when relations are known. (A flag was added rather than changing the default, per author preference, to avoid disrupting existing users.)

**R1.6 — Beyond phase overlap, discuss three further sources of relative-chronological uncertainty (within-phase, phase-assignment, phase-boundary). Cite Crema & Kobayashi 2020; Crema 2024.**
`[DONE][PAPER]` Added a subsection *Sources of relative-chronological uncertainty* that names phase overlap as the principal challenge already modelled, then treats the three further types **honestly**: **phase-assignment** uncertainty is represented *directly* by the soft posteriors p_ik / entropy / noise probability; **within-phase** uncertainty is *partly* handled by `chrono_uncertainty` (interval-width propagation, variance (t_max−t_min)²/12) but the model does not yet date events within a phase; **phase-boundary** uncertainty is *not yet* modelled (phases are latent clusters, not OxCal-style Bayesian-dated intervals) and is flagged as future work, with the `chronology_from_rcarbon()`/`chronology_from_oxcal()` adapters as the current bridge. Crema & Kobayashi (2020) and Crema (2024) are cited in text and added to the reference list (real titles fetched via Crossref).

**R1.7 — Case study: discuss the detected intrusions archaeologically (field/material-culture interpretation; are they genuine and plausible?).**
`[DONE][PAPER]` Added an archaeological reading of the noise-component outliers to the case study (Section "Intrusion detection"). The ten flagged finds cluster in five stratigraphic units (US 21, 173, 204, 276, 302) and are predominantly non-ceramic — building material, an amphora fragment, and an isolated metal object — standing out within the ceramic-dominated assemblages of their macro-events; here it is material-class atypicality and spatial isolation, not chronological mismatch, that drives the flag (consistent with redeposited construction debris and residual items in accumulation deposits). Examined against the excavators' field knowledge (R.M., excavation director; M.C.), these flags were judged archaeologically coherent. We also make explicit that, because the chronology is unit-tied, the directional residual (older-than-context) vs latent-feature (younger-than-context) classification from `detect_intrusions()` reports all dated finds as in-context here, and becomes informative only with per-find typological dating.

---

## Reviewer 2 — Lise Bellanger

### Major

**R2.1 — Manuscript structure.** Numbered sections + cross-references; reorganise to separate *archaeological context / theoretical framework / software implementation / case study / discussion+conclusion*; restructure the list-like *Limitations and future developments*; explicitly discuss using a GMM when variables are not all continuous; in *Perspectives*, compare with alternative approaches; add references; **number the equations**; use appendices to lighten the text (e.g. move some case-study figures).
`[TODO][PAPER]` Adopted numbered sectioning with the five-part structure requested; equations are numbered; *Limitations* rewritten as structured prose; a *Perspectives* paragraph compares the SEF with alternatives (mixture-of-categorical/latent-class models, spatial clustering, Bayesian chronological modelling); selected case-study figures moved to an Appendix; additional references added.

**R2.2 — Conceptual framework.** Define early, with concrete examples: *find*, *cultural class*, *unit purity*, *typology*, *Excavation Stratigraphic Energy*. Specific lines: palimpsest stratigraphic vs interpretative (55–59); *find* = artefact/ecofact/aggregated unit (74–79); *decomposition* (79); *cultural class* — typological/functional/material (80–82); *low purity* (88); *typology* (91–92); SEI matrix as regulariser + ESE definition (127–130).
`[DONE][PAPER]` Added a new subsection *Key concepts and definitions* at the head of the Methodological framework, defining *find* (artefact/ecofact, with the unit-resolution fallback), *stratigraphic unit*, *cultural class*, *typology*, *decomposition*, and *unit purity*, each with an archaeological example, and making the strictly-stratigraphic vs interpretative/landscape senses of "palimpsest" explicit. ESE is now formally defined (Eq. ESE) and tied to the SEI affinities and to the fact that it is computed independently of the fitted phases (see R2.3).

**R2.3 — Statistical modelling.** State and justify model assumptions; give SEI a formal definition (range, interpretation, behaviour) and clarify whether it is novel or from the literature; rename "Diagonal Gaussian Mixture EM" → "Stratigraphic Entanglement Field (SEF)"; state modelling objectives (clustering / unsupervised / phase-number estimation); introduce the GMM pedagogically with its probabilistic formulation; **explicitly explain handling of categorical variables**; justify the diagonal covariance (problematic for $p \gg 20$); justify k-means (line 136) and its role; add a schematic of the modified EM; describe synthetic-data generation in detail (lines 217–218, assumptions + procedure).
`[DONE][PAPER]` The methods section is renamed *The Stratigraphic Entanglement Field (SEF) model* and reorganised into numbered subsubsections: *Modelling objective* (states the unsupervised-clustering / soft-assignment / K-estimation goals), *Probabilistic formulation* (the mixture density, Eq. mixture, with both structural assumptions stated), *Treatment of the categorical class variable*, *Estimation by Expectation-Maximisation* (with **Algorithm 1**, the requested schematic of the modified EM), *Optional refinements*, and *Why a diagonal covariance?*. SEI is now explicitly stated to be **newly proposed** (not from the literature), with its [0,1]-scaled components and interpretation. **Categorical handling** is foregrounded as the per-phase multinomial (Gaussian × multinomial); we add the point that, because class is multinomial rather than one-hot Gaussian, the numeric block stays low-dimensional ($p \sim 5$–8), which directly answers the $p \gg 20$ / diagonal-covariance concern — now justified for the operative regime and its limit acknowledged. The role of **k-means** (initialisation only, n_init restarts) is clarified. The **synthetic-data generation** (`archaeo_sim()`) is documented in full (box, Gaussians, Beta taph. score, mixing perturbation) with a pointer to a Supplementary sensitivity analysis (ARI vs mixing).

**R2.4 — "Diagnostic statistics" presents four indices, not three; make notation consistent with the SEF framework.**
`[TODO][PAPER]` Confirmed and corrected: the section now consistently presents the indices it lists (SEI, ESE, PDI, Shannon entropy, and the intrusion/noise posterior), with unified notation. *(In the submitted text the prose said "three" while four bullets followed.)*

### Minor

**R2.5 — Line 59–63: long sentence; split for readability.** `[TODO][PAPER]` Done.
**R2.6 — Line 81: GMM introduced without detail; refer the reader to Methods.** `[TODO][PAPER]` Added a forward reference.
**R2.7 — Line 126: clarify "normalised" (type; global vs per variable/domain).** `[TODO][PAPER]` Clarified: each of the four SEI components is min–max scaled to [0,1] within the dataset (per domain) before weighting.

---

## Summary of changes by type
- **Code (package):** R1.4 `archaeo_sim()` floor; R1.5 `harris_from_contexts()` default + docs; R1.2 help pages for BIC/ICL/PDI.
- **Paper (exposition/structure):** R1.1, R1.2, R1.3, R1.6, R1.7, R2.1, R2.2, R2.3, R2.4, R2.5–R2.7; mandatory PCI sections.
- **Co-author input needed:** R1.7 (archaeological interpretation of intrusions — R.M./M.C.); Funding statement.
