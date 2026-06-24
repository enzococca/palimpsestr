# Revision Roadmap — PCI Archaeo #1019, Round #1

Execution plan mapped to the reviewer points in `response-to-reviewers.md`.
Target deliverables for resubmission: (1) revised preprint on Zenodo; (2) point-by-point reply (PDF); (3) tracked-changes manuscript.

Legend: **P** = paper (`docs/palimpsestr-paper.Rmd`) · **C** = code/package · **D** = docs/help · **X** = needs co-author input.

---

## Phase 0 — Mandatory PCI sections (fast, do first)
| # | Item | Type | Where |
|---|------|------|-------|
| 0.1 | Add **Funding** section | P, X | new section before *Conflict of interest* |
| 0.2 | "The author declares" → **"The authors declare"** | P | *Conflict of interest* |
| 0.3 | Add **Supplementary information** section + deposit `introduction.Rmd` with DOI | P | new section + Zenodo |
| 0.4 | Verify dataset + repro script archived with DOI/metadata (not GitHub-only) | P | *Data and code availability* |
| 0.5 | Update package version string in paper (0.21.0 → release used for resubmission) | P | several |

## Phase 1 — Quick, certain fixes
| # | Item | Type | Where |
|---|------|------|-------|
| 1.1 | **Number all sections** (1, 1.1, …) | P | YAML `number_sections: true` + restructure headers |
| 1.2 | **Number equations** (SEI, PDI, entropy, EM updates) | P | `\begin{equation}` + `\label` |
| 1.3 | **Fix "Three statistics" → consistent count** (R2.4) | P | *Diagnostic statistics* |
| 1.4 | `archaeo_sim()` floor dates (R1.4) | C | `R/simulate.R` + test + NEWS |
| 1.5 | Split long sentence L59–63 (R2.5) | P | Introduction |
| 1.6 | Forward-ref GMM at first mention (R2.6) | P | Introduction |
| 1.7 | Clarify "normalised" = per-domain min–max to [0,1] (R2.7) | P | SEI section |

## Phase 2 — Conceptual definitions & framework (R2.2)
| # | Item | Type |
|---|------|------|
| 2.1 | Early *Definitions* paragraph: find, cultural class, unit purity, typology, ESE — each w/ example | P |
| 2.2 | Make explicit: strictly-stratigraphic (taphonomic) vs interpretative/landscape "palimpsest" | P |
| 2.3 | Formally define **ESE** and its link to the SEI matrix and the model | P |

## Phase 3 — Methods rewrite (R1.1, R2.1 structure, R2.3) — the core
| # | Item | Type |
|---|------|------|
| 3.1 | Reorganise into 5 numbered parts: archaeological context / theoretical framework (SEF) / software / case study / discussion+conclusion | P |
| 3.2 | Rename "Diagonal Gaussian Mixture EM" → **"Stratigraphic Entanglement Field (SEF)"**; state objectives (clustering / unsupervised / K estimation) | P |
| 3.3 | Pedagogical GMM intro: probabilistic formulation + assumptions, stated explicitly | P |
| 3.4 | **SEI**: formal definition, [0,1] range, interpretation, behaviour; state it is newly proposed | P |
| 3.5 | **Categorical-variable handling**: foreground the per-phase multinomial (Gaussian×multinomial) as the answer to heterogeneous variables — already in code | P |
| 3.6 | Justify **diagonal covariance** (p∼5–8 regime); acknowledge p≫20 limitation | P |
| 3.7 | Justify **k-means** role (init only, restarts) — L136 | P |
| 3.8 | Add **schematic diagram of the modified EM** algorithm | P (new figure) |
| 3.9 | Document **synthetic-data generation** in detail (assumptions + procedure) — L217–218 | P |
| 3.10 | Define/reference **BIC, ICL, PDI** in text (R1.2); reference "BIC minimum"; stress stats+expert judgement | P |
| 3.11 | Fold key diagnostic explanations from `introduction.Rmd` into the paper (R1.3) | P |

## Phase 4 — Chronological uncertainty (R1.6)
| # | Item | Type |
|---|------|------|
| 4.1 | Paragraph on 3 uncertainty types: within-phase / phase-assignment / phase-boundary | P |
| 4.2 | Map each to package features (`chrono_uncertainty`, soft posteriors/entropy, `type_longevity`/`chronology_from_*`) | P |
| 4.3 | Cite **Crema & Kobayashi 2020** (10.1016/j.jas.2020.105136) + **Crema 2024** (10.1111/arcm.12984) | P |

## Phase 5 — harris_from_contexts redesign (R1.5)
| # | Item | Type |
|---|------|------|
| 5.1 | Decide design: default that does not flag normal verticality configs / new option | C |
| 5.2 | Document assumption + failure modes (verticality violations, multi-period fills) | D |
| 5.3 | Expand *Limitations* to treat these as expected, not pathological | P |

## Phase 6 — Discussion, perspectives, appendices (R1.7, R2.1)
| # | Item | Type |
|---|------|------|
| 6.1 | **Archaeological interpretation of detected intrusions** in case study | P, **X** (R.M./M.C.) |
| 6.2 | Rewrite *Limitations & future developments* as structured prose (not a list) | P |
| 6.3 | *Perspectives*: compare SEF with alternatives (latent-class/mixture-of-categorical, spatial clustering, Bayesian chrono modelling) | P |
| 6.4 | Move selected case-study figures to an **Appendix** | P |
| 6.5 | Add supporting references throughout | P |

## Phase 7 — Help pages (R1.2)
| # | Item | Type |
|---|------|------|
| 7.1 | Document BIC/ICL/PDI in `compare_k()`, `fit_sef()`, `pdi()` roxygen | D |
| 7.2 | Document `harris_from_contexts()` assumption (from 5.2) | D |
| 7.3 | `roxygen2::roxygenise()` + regenerate manual | C |

## Phase 8 — Finalise
| # | Item |
|---|------|
| 8.1 | Knit revised paper; proofread; check equation/section cross-refs |
| 8.2 | Produce **tracked-changes** version (export to .docx or latexdiff) |
| 8.3 | Finalise `response-to-reviewers.md` ("will add" → "have added, §X, lines Y") → export PDF |
| 8.4 | Bump version, NEWS, tag/release; deposit new preprint on Zenodo → get new DOI |
| 8.5 | Update PCI: edit article data + upload reply + tracked manuscript + resubmit |

---

## Open decisions for you (Enzo)
1. **Funding** — what do we state? (grant/none)
2. **harris_from_contexts** (5.1) — change the default, or add a flag and keep the default? Affects existing users.
3. **Supplementary** — deposit `introduction.Rmd` as SI with DOI? (recommended)
4. **R1.7 intrusions** — get R.M./M.C. note on the 10 flagged Poggio Gramignano finds.
5. **Tracked changes** — latexdiff on the .Rmd→.tex, or maintain a parallel .docx? (latexdiff cleaner for this LaTeX paper)
