# Portal metadata — paste-ready

> Use the abstract of the **recommended version** (below) so the portal matches the reviewed
> manuscript. Keywords, authors and statements follow.

## Title
palimpsestR: An R Package for the Identification and Probabilistic Decomposition of Archaeological Palimpsests

## Authors (order + ORCID + affiliation)
1. **Enzo Cocca** (corresponding) — Independent Researcher, Italy — ORCID 0000-0002-8096-2810 — enzo.ccc@gmail.com
2. **Roberto Montagnetti** — Università degli Studi dell'Aquila, Dipartimento di Scienze Umane, L'Aquila, Italy — ORCID 0009-0004-2089-3122
3. **Maurizio Cattani** — Università di Bologna, Dipartimento di Storia Culture Civiltà, Bologna, Italy — ORCID 0000-0003-2607-3657

## Keywords
archaeological methods; palimpsest analysis; probabilistic clustering; mixed-type mixture model; Gaussian mixture model; Harris Matrix; R package; Shiny application; post-excavation stratigraphy

## Abstract (recommended version — plain text)
This paper presents palimpsestr, an R package and Shiny application dedicated to the identification and probabilistic decomposition of archaeological palimpsests — deposits in which material from multiple occupation phases is superimposed and partially intermixed. Such deposits represent one of the most persistent analytical challenges in field archaeology, and existing tools (Harris Matrix recording systems, k-means spatial clustering, GIS overlays) typically address only one evidence domain at a time. palimpsestr implements the Stratigraphic Entanglement Field (SEF) framework, which integrates four evidence domains — horizontal coordinates, vertical elevation, chronological range, and cultural class — into a single diagonal-covariance Gaussian mixture model fitted by Expectation-Maximisation. The model is augmented with optional taphonomic weighting and stratigraphic entanglement penalties derived from the Harris Matrix. Three interpretable diagnostics — the Stratigraphic Entanglement Index (SEI), Excavation Stratigraphic Energy (ESE), and Palimpsest Dissolution Index (PDI) — allow practitioners to assess deposit coherence, detect potentially redeposited finds, and evaluate the reliability of chronological attribution at the find, unit, and phase levels. The framework supports archaeologists in three principal interpretive tasks: assessing the reliability and coherence of stratigraphic units (distinguishing genuine palimpsests from recording errors); identifying intrusive finds with respect to both their typology and their stratigraphic position; and providing the basis for future statistical estimation of type longevity (the duration of use of pottery and material categories). The package supports CSV/TSV/Excel/SQLite/PostgreSQL data import, GIS export to GeoPackage via the sf package, publication-quality plots via ggplot2 and interactive plots via plotly, and includes a built-in Shiny dashboard for non-programmatic use. Since the original release the statistical core has been substantially refined: the cultural class is now modelled as a per-phase categorical distribution (a Gaussian-times-multinomial mixed-type mixture) rather than as one-hot Gaussian columns, so that stratigraphic units are no longer split across phases by typology; the spatial and vertical similarity kernels are bounded and scale-invariant; the domain weights enter the likelihood and can be cross-validated; a uniform noise component yields a genuine posterior probability of intrusion and shields phase estimates from outliers; and optional treatments propagate per-find dating uncertainty and apply the stratigraphic constraint as a dynamically updated Neighborhood-EM field. The intrusion diagnostic also distinguishes the direction of the chronological mismatch (residual vs. latent-feature), and a companion helper (recommend_setup()) inspects a dataset and reports when its recording resolution limits the achievable inference. We present an overview of the application's functions and demonstrate its use through a case study at the multi-period Roman villa of Poggio Gramignano (Lugnano in Teverina, Italy). We also discuss the methodological assumptions of the SEF framework — in particular the implicit assumption of horizontal stratigraphy and the dependence of phase resolution on the spatial and chronological recording resolution of the data — and outline planned developments.

## Statements
- **Funding:** This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.
- **Conflict of interest:** The authors declare no conflict of interest relating to the content of this article.
- **Data and code availability:** palimpsestr (v0.24.0) is available at https://github.com/enzococca/palimpsestr (MIT licence) and archived on Zenodo, concept DOI 10.5281/zenodo.19881542. The villa_romana case-study dataset is bundled with the package.

## Author contributions
E.C. designed and implemented the palimpsestr package, the Stratigraphic Entanglement Field framework, and the Shiny dashboard, and wrote the manuscript. R.M. directed the archaeological excavations at Poggio Gramignano (Lugnano in Teverina), provided the dataset and the taphonomic scores per stratigraphic unit, validated the model output against on-site stratigraphic interpretation, and contributed to the discussion of methodological limitations and feature-space improvements. M.C. provided methodological supervision throughout the development of the framework, contributed to the archaeological interpretation of the diagnostics, highlighted the implicit assumption of horizontal stratigraphy, introduced the residual vs. latent-feature intrusion distinction, and identified the future dendrochronological validation context.

## Recommendation / linkage (for the portal fields)
- Recommended by: **PCI Archaeo**
- Recommendation DOI: **10.24072/pci.archaeo.101019**
- Recommender: **Alessio Palmisano**
- Recommended preprint DOI (ver.2, the version recommended by PCI): **10.5281/zenodo.20848849**
  (latest Zenodo version, v3, with PCI badge: 10.5281/zenodo.21534452)
