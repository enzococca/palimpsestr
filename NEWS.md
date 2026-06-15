# palimpsestr 0.22.1

## Per-US taphonomic score from the chronology table

- `read_pyarchinit()`: when `taf` is `NULL` (the default) and the
  `chronology_table` (`palimpsest_chronology`) carries a per-US `taf` column,
  those values now populate `taf_score` (units without a value default to 0.5).
  This lets the pyArchInit QGIS algorithms (fit / intrusions / report) honour a
  per-US taphonomic score stored alongside the absolute chronology, without any
  change to the calling `.rsx`.

# palimpsestr 0.22.0

## Per-US absolute chronology (OxCal)

- `read_pyarchinit()` gains `chronology_table` (default `"palimpsest_chronology"`):
  an optional per-US table of absolute dating with columns `sito`, `area`, `us`,
  `start`, `end` (calendar years, BCE negative) — e.g. OxCal-calibrated ranges
  from `chronology_from_oxcal()`. When present it overrides the free-text
  `datazione` for the matching units (envelope of multiple dates per US), joined
  on `(sito, area, us)` with a `(sito, us)` fallback.

## PostgreSQL / PostGIS in the QGIS algorithms

- The three Processing R scripts (`r:palimpsestrfit`, `r:palimpsestrintrusions`,
  `r:palimpsestrreport`) accept an optional `PG_connection` parameter (a libpq
  DSN) and read from PostgreSQL/PostGIS when it is set, falling back to the
  SQLite/Spatialite `Database_file` otherwise. `read_pyarchinit()` itself was
  already database-agnostic (plain DBI). The pyArchInit tab passes its active
  connection string, so PostgreSQL pyArchInit projects are supported.

# palimpsestr 0.21.0

## Narrated report export

- New `export_sef_report()`: assembles a complete, narrated report from a
  `sef_fit` — the `report_sef()` interpretive text, all applicable `gg_*`
  diagnostic plots, and diagnostic tables (US summary, phase-transition matrix,
  per-US phase assignments, model statistics) — rendered to PDF and/or DOCX from
  a shipped RMarkdown template (`inst/rmarkdown/sef_report.Rmd`), in Italian or
  English. A markdown narrative sidecar (`<file>.md`) is always written; when
  pandoc/LaTeX are unavailable the function degrades gracefully to that
  markdown plus PNG figures. `tinytex` added to Suggests.

## Corrected pyArchInit data sourcing

- `read_pyarchinit()` now sources elevation per archaeological practice: US
  elevation from `pyarchinit_quote` (mean `quota_q` per unit, with a
  `(sito, us)` join fallback), material elevation from `quota_usm`, and pottery
  inheriting its US elevation. The deprecated `us_table.quota_abs` form field is
  no longer used.
- New `source = "both"|"materials"|"pottery"` argument: reads
  `inventario_materiali_table` and/or `pottery_table` (pottery `class` taken
  from the first non-empty of `ware`/`material`/`form`); a `find_source` column
  is added when `source = "both"`.

## QGIS Processing

- New `r:palimpsestrreport` algorithm (`palimpsestr_report_db.rsx`): database →
  fit → narrated PDF/DOCX report, with Source / Language / Format parameters.
- `r:palimpsestrfit` and `r:palimpsestrintrusions` gain a `Source` parameter.
- The `villa_romana` example database builder now populates `quota_usm`,
  `pottery_table`, and `pyarchinit_quote` to exercise the corrected sourcing.

# palimpsestr 0.20.0

## pyArchInit database bridge

- New `read_pyarchinit()`: reads the `inventario_materiali_table` (finds) and
  `us_table` (stratigraphic units) of a pyArchInit database (SQLite/Spatialite
  or PostgreSQL/PostGIS) and assembles a `fit_sef()`-ready data.frame. Each
  find inherits its stratigraphic unit's coordinates (centroid of the US
  polygon supplied via `us_geometry`, e.g. the `pyunitastratigrafiche` layer),
  elevation, chronology, and context. Free-text period strings (e.g. "II sec.
  a.C.", "I/II sec. d.C.", "età romana") are resolved to numeric `date_min`/
  `date_max` by a built-in archaeological-date parser (extendable via
  `date_labels`); units without a digitised polygon fall back to synthetic
  per-unit coordinates so the analysis can still run at unit resolution.
  Validated against real pyArchInit Spatialite databases.

# palimpsestr 0.19.0

## New features

- `chronology_from_oxcal()`: OxCal counterpart of `chronology_from_rcarbon()`.
  Accepts either an `oxcAAR::oxcalCalibrate()` result
  (`oxcAARCalibratedDatesList`, reduced by HPD / median+IQR / weighted mean) or
  a generic data.frame of `start`/`end` calibrated ranges, and returns the
  `date_min`/`date_max`/`date_mid` columns expected by `fit_sef()`. `oxcAAR`
  is in Suggests.
- `detect_intrusions()` gains an `intrusion_threshold` argument and an
  `intrusion_type` factor column that couples the intrusion magnitude (the
  noise-component posterior) with the chronological direction: `not_flagged`,
  `residual` (flagged and older-than-context), `latent_feature` (flagged and
  younger-than-context), or `outlier_in_context`. Existing columns are
  unchanged.
- `gg_outliers()` now colours each find by `intrusion_type`.

# palimpsestr 0.18.0

## New diagnostic plots

- `gg_phase_composition()`: per-phase class profile from the multinomial model
  (`cat_prob`), as a heatmap (default) or stacked bars.
- `gg_direction()`: diverging lollipop of the directional intrusion diagnostic
  (residual / older-than-context vs latent-feature / younger-than-context),
  with a graceful placeholder when every find is in context.
- `gg_outliers()`: non-spatial ranking of the intrusion probability, using the
  noise-component posterior when available, otherwise the heuristic composite.
- `gg_unit_coherence()`: within-stratigraphic-unit phase coherence, with an
  optional comparison fit (e.g. the legacy Gaussian-class model).

All four are exposed in a new "Model profile" tab of the Shiny app and are
compatible with `as_plotly()`.

# palimpsestr 0.17.1

Bug fixes for noise-component model selection (addresses code review on the
0.17.0 pull request).

- **Held-out scoring now includes the noise component.** `cv_sef()` and
  `optimize_weights()` scored only the Gaussian phases of a `noise = TRUE`
  fit, so the uniform component was omitted and outlier held-out points were
  mis-scored (their mass was missing and the Gaussian weights no longer summed
  to one). The fit now stores the noise mixing weight (`noise_weight`), its
  log-density (`noise_logdens`), and the global class frequencies
  (`cat_global`), and both evaluators add the noise component to the held-out
  likelihood.
- **BIC/ICL now count the noise mixing weight** as one extra parameter, so
  `compare_k(..., noise = TRUE)` no longer gives noise models a free parameter
  relative to non-noise fits.

# palimpsestr 0.17.0

## Noise / outlier mixture component

- **New `fit_sef(noise = FALSE, noise_prior = 0.05)`.** When `TRUE`, a uniform
  background component is added to the mixture (Fraley & Raftery, 1998). Finds
  that fit no Gaussian phase accumulate posterior probability on it, so the fit
  gains a `noise_prob` vector that is a **genuine posterior probability of
  being an outlier/intrusion** — unlike the previous composite score, it is not
  forced onto `[0, 1]` by rescaling, so a clean assemblage need not contain a
  probability-1 intrusion. As a side effect the noise component absorbs extreme
  finds, keeping them out of the phase centroid and variance estimates
  (robust phase estimation). The reported `phase_prob` is then the distribution
  over phases conditional on a find not being noise.

- **`detect_intrusions()` now returns the model-based posterior** in its
  `intrusion_prob` column whenever the fit was made with `noise = TRUE`,
  falling back to the heuristic composite (rescaled entropy + energy +
  inverse local SEI) otherwise. The directional columns are unchanged.

- The option is propagated through `bootstrap_sef()` and exposed in the Shiny
  app.

# palimpsestr 0.16.0

This release adds two model-based treatments of evidence that the EM engine
previously simplified away: per-find dating uncertainty, and the
stratigraphic constraint. Both are opt-in and off by default, so existing
fits are unchanged.

## Per-find chronological uncertainty

- **New `fit_sef(chrono_uncertainty = FALSE)`.** When `TRUE`, a find's dating
  interval width is propagated into the likelihood instead of collapsing the
  interval to its mid-point. The mid-date is modelled as observed with a
  measurement variance `(date_max - date_min)^2 / 12` (a uniform prior over
  the interval), added to the temporal dimension's component variance in the
  E-step and removed from the component variance estimate by method-of-moments
  deconvolution in the M-step. Finds with wide dating ranges then rely less on
  chronology and carry their dating uncertainty into the phase probabilities.
  Unlike `chrono_precision` (which adds `1/tspan` as a feature), this is a
  model-based treatment of the uncertainty. The internal `diag_log_density()`
  gains an `extra_var` argument.

## Dynamic (Neighborhood-EM) stratigraphic constraint

- **New `fit_sef(strat_dynamic = FALSE, strat_beta = 1)`.** The stratigraphic
  context constraint was previously a static penalty frozen from the k-means
  initialisation: if the EM moved assignments, the penalty still reflected the
  starting clusters. With `strat_dynamic = TRUE` it becomes a Neighborhood-EM /
  hidden-Markov-random-field term recomputed from the current phase posteriors
  at every E-step (Ambroise & Govaert, 1997): each find is rewarded for
  sharing the phase of its stratigraphic-unit neighbours, so units stay
  coherent as the fit evolves. `strat_beta` sets the field strength; any
  `harris` penalty continues to be applied statically.

- Both options are propagated through `bootstrap_sef()` and exposed in the
  Shiny app.

# palimpsestr 0.15.0

## Mixed-type likelihood for the class label

- **The class label is now modelled as a per-phase categorical distribution**
  (a Gaussian-times-multinomial mixed-type mixture), controlled by the new
  `class_model` argument to `fit_sef()`. `"multinomial"` is the default;
  `"gaussian"` reproduces the previous one-hot behaviour and is kept for
  backward compatibility.

  Up to v0.14.0 the class label was one-hot encoded into the Gaussian feature
  block. Because zero/one dummies make finds of different classes almost
  perfectly separable, this had two undesirable effects: per-component
  variances on the dummy columns collapsed to the floor (over-confident
  assignments), and — more importantly — the class label could **split a
  single stratigraphic unit across several phases**, even though a unit is one
  depositional event. Under the multinomial model the numeric evidence
  (space, depth, chronology) drives the phase assignment and the class label
  only tilts it softly, so units with finds of mixed classes stay coherent.

- **New `class_smoothing` argument** (default 1): a Dirichlet pseudo-count,
  shrunk towards the global class frequencies, that keeps every per-phase
  class probability strictly positive. `0` removes the shrinkage.

- **New `phase_composition()` accessor**: returns the estimated per-phase
  class profile (row `j` = probability that a find in phase `j` belongs to
  each class). The model-based, directly interpretable counterpart of a
  phase-by-class cross-tabulation, stored on the fit as `$cat_prob`.

- `cv_sef()`, `optimize_weights()`, and `bootstrap_sef()` propagate the
  class model and include the categorical term in the held-out likelihood;
  `reorder_phases()` reorders the per-phase class profiles; BIC/ICL use the
  correct mixed-type parameter count (`k(L-1)` categorical parameters instead
  of treating the dummies as Gaussian).

Note: on datasets recorded at stratigraphic-unit resolution (finds inheriting
US-centroid coordinates and US-tied chronology, e.g. the bundled
`villa_romana`), phase assignments remain near-certain under either model
because there are only as many distinct numeric points as there are units;
that certainty reflects the data resolution, not the model. The multinomial
model's benefit there is keeping units coherent and yielding interpretable
per-phase class profiles.

# palimpsestr 0.14.0

This release is a statistical-correctness pass on the core engine. Several
parameters that were silently inert are now active, and three diagnostics
that depended on the measurement scale or were distorted by single outliers
have been made scale-invariant and robust. Fits produced with non-default
weights, and the output of `optimize_weights()`, will differ from 0.13.0;
default-weight fits are affected only by the more careful BIC/ICL accounting
and the new multi-start default.

## Statistical corrections

- **Domain weights now enter the model.** `weights = c(ws, wz, wt, wc)`
  previously affected only the SEI matrix, leaving `fit_sef()` — and hence
  the cross-validated objective of `optimize_weights()` — completely
  invariant to them. The weights now scale the corresponding standardised
  feature dimensions (`ws`: planar coordinates, `wz`: depth, `wt`:
  chronology-derived features, `wc`: the one-hot class block), so they
  genuinely influence the mixture likelihood. Weights must be strictly
  positive; they are stored on the fit and propagated by `bootstrap_sef()`,
  `cv_sef()`, and `optimize_weights()`.

- **New `var_structure` argument to `fit_sef()`** (`"diagonal"`, default, or
  `"spherical"`). Free diagonal variances absorb any per-dimension rescaling,
  so domain weights are not identifiable under a diagonal fit; the
  `"spherical"` structure (one shared variance per component) makes them
  identifiable. `optimize_weights()` now selects weights under
  `var_structure = "spherical"` and corrects each held-out log-likelihood by
  the Jacobian of the weighting (`n_test * sum(log(w_d))`) so configurations
  with different weights are directly comparable.

- **Taphonomic down-weighting is no longer applied twice.** `fit_sef()`
  passed `weights_obs = 1 - 0.5 * taf` *and* the EM multiplied by
  `(1 - 0.5 * taf)` again, so a find with `taf = 1` contributed at weight
  0.25 instead of the documented 0.5. It is now applied exactly once.

- **Removed the inert `taf_weight_e` term.** It subtracted a per-row constant
  from every mixture component, which cancels exactly in the softmax and so
  had no effect on the posteriors (only on the reported log-likelihood). The
  dead `taf`, `taf_weight_m`, and `taf_weight_e` arguments were removed from
  the internal EM.

- **BIC and ICL are now built from the true unpenalized mixture
  likelihood.** With a stratigraphic penalty the EM objective is a penalized
  criterion, not a log-likelihood; using it biased comparisons across `K` and
  across context settings. `model_stats$loglik` is the unpenalized
  likelihood; the penalized EM objective is available as
  `model_stats$loglik_penalized`. The `"spherical"` structure also gets its
  own (smaller) parameter count in the BIC.

- **SEI and ESE spatial/vertical components are now scale-invariant and
  robust to duplicate coordinates.** The former `1/d` kernel normalised by
  its maximum was dominated by the single closest pair, so one near-duplicate
  pair (e.g. two finds sharing a grid square) collapsed every other spatial
  affinity towards zero. Both components now use a bounded exponential kernel
  `exp(-d / h)` with a data-driven bandwidth `h` (median pairwise
  separation). `ese()` consequently no longer depends on the measurement unit
  of the site grid. The `eps` and `z_floor` arguments of `sei_matrix()` are
  deprecated and ignored.

- **`detect_intrusions()` directional envelope is now robust.** The
  leave-one-out unit envelope used a strict min/max over the other finds, so
  a single chronological outlier in a unit could mask genuinely
  residual/intrusive neighbours. A new `envelope` argument (default
  `c(0.05, 0.95)`) uses robust quantiles; pass `c(0, 1)` for the previous
  min/max behaviour.

- **`n_init` now defaults to 5** (was 1) and k-means initialisation uses
  `nstart = 5`. The EM objective is multimodal, so a single start risked a
  poor local optimum.

# palimpsestr 0.13.0

## New features

- `detect_intrusions()` now returns two new columns:
  - `direction`: factor classifying each find as `older_than_context`,
    `in_context`, or `younger_than_context` based on interval overlap
    between the find's chronology and the leave-one-out chronological
    envelope of its stratigraphic unit.
  - `chrono_gap`: signed offset in years (negative = residual,
    positive = latent intrusion, 0 = in-context).
  Backward compatible: the existing `intrusion_prob` column is unchanged.

- New `type_longevity()` function: per-class temporal envelope of the
  phases with mean posterior weight above a configurable threshold
  (default 0.1). Returns `longevity_min`, `longevity_max`,
  `longevity_span`, `dominant_phase`, `n_finds`, and a list-column
  `weight_matrix` with the per-phase posterior weights.

- New `gg_longevity()` Gantt-style plot for the `type_longevity()` output.

- New `chronology_from_rcarbon()` adapter: converts an
  `rcarbon::CalDates` object into the `date_min`/`date_max`/`date_mid`
  columns expected by `fit_sef()`. Three reduction methods: HPD
  (default), `median_iqr`, `weighted_mean`. BCE/CE sign convention
  by default. `rcarbon` is in Suggests.

- New `recommend_setup()` helper: inspects a dataset and returns
  recommended `fit_sef()` flags (`class_scale`, `taf_as_feature`,
  `chrono_precision`, `residuality`) based on diagnostic heuristics
  (number of classes, singleton classes, US-centroid coordinates,
  US-tied chronology, taf_score informativeness). Reports caveats
  (e.g. directional intrusions will all be `in_context` when the
  chronology is US-tied) and prints a ready-to-paste `fit_sef()`
  call. Closes #2.

## Documentation

- New vignette subsections demonstrating the three features above.

## Shiny app

- "Intrusions" tab now shows `direction` and `chrono_gap`.
- New "Type longevity" tab with plot, table, and xlsx download.
- "Data import" tab includes an optional rcarbon calibration panel.

## Compatibility

No breaking changes. All v0.12.0 user code continues to work as before.

# palimpsestr 0.10.0

## New features

- New plots: `gg_cv()` (cross-validation diagnostics), `gg_bootstrap()` (confidence interval forest plot), `gg_weights()` (weight sensitivity heatmap).
- `villa_romana` dataset: realistic 300-find Roman villa with 4 phases, bioturbation, construction cuts, and residual pottery.

## Bug fixes and improvements

- **Cross-validation fix**: test data now standardised using training set center/scale, eliminating bias in held-out log-likelihood (`cv_sef()`, `optimize_weights()`).
- **Bootstrap fix**: `reorder_phases()` applied to each replicate before ARI computation, solving label switching in confidence intervals.
- **ESE neighbourhood fix**: when `neighbourhood` is set, energy is divided by actual neighbour count instead of n-1, removing peripheral bias.
- **Numerical stability**: EM now emits a warning when degenerate posteriors occur (instead of silently replacing with uniform).
- **Feature matrix**: `feature_matrix()` now accepts external `center`/`scale` parameters and preserves scaling attributes for proper cross-validation.
- Taphonomic weighting parameters (M-step: 0.5, E-step: 0.15) now configurable in EM engine.
- Context penalty weight (default 0.25) now configurable and documented.
- SEI documentation: explicit warning that absolute values are not comparable across datasets.

# palimpsestr 0.9.0

## New features

- Model validation: `adjusted_rand_index()` and `confusion_matrix()` for comparing estimated vs true phases.
- Structured exports: `us_summary_table()` aggregates diagnostics per stratigraphic unit; `phase_transition_matrix()` reveals vertical phase ordering; `export_results()` writes all results to CSV.
- Harris Matrix tooling: `harris_from_contexts()` auto-generates penalties from depth ordering; `read_harris()` imports CSV edge lists; `validate_phases_harris()` checks phase-stratigraphy consistency.
- New plots: `gg_convergence()` (EM trace), `gg_phase_profile()` (depth vs phase), `gg_confusion()` (confusion matrix heatmap with ARI).
- `bootstrap_sef()` for bootstrap confidence intervals on PDI, entropy, energy, loglik, and ARI.
- `cv_sef()` for k-fold cross-validation of phase count selection.
- `optimize_weights()` for data-driven SEI weight estimation via grid search cross-validation.
- `reorder_phases()` relabels phases by mean depth (phase 1 = deepest/oldest) to solve label switching.
- `sei_sparse()` wrapper for sparse SEI computation on large datasets (n > 1000).

## Improvements

- `fit_sef()` gains `n_init` parameter for multiple random initialisations (default: 1). Best run by log-likelihood is retained.
- `fit_sef()` default `em_iter` increased from 25 to 100 for better convergence.
- EM convergence tracking: `$converged` flag in `sef_fit` objects; warning issued on non-convergence.
- ICL (Integrated Complete-data Likelihood) added to `fit_sef()` and `compare_k()` model selection metrics.
- SEI matrix components normalised to [0, 1] before weighting for consistent dimensional contribution.
- `sei_matrix()` gains `max_dist` parameter for sparse computation (zeroes beyond threshold).
- `sei_matrix()` and `ese()` fully vectorized (50-100x faster on large datasets).

# palimpsestr 0.8.0

## New features

- Interactive plotly support via `as_plotly()` with enriched archaeological tooltips (ID, context, phase, dating, class, entropy, energy, intrusion probability).
- All `gg_*` functions now embed tooltip data for seamless plotly conversion.

## Documentation

- Complete roxygen documentation with `@examples`, `@family`, and `@seealso` for all exported functions.
- Vignette rewritten as paper-ready methodological article with formal SEI/ESE/PDI definitions, worked example, and comparison with traditional approaches.
- Added `_pkgdown.yml` for documentation website.
- Added `cran-comments.md` for CRAN submission.

# palimpsestr 0.6.0

## New features

- Added three built-in demo datasets: `demo_easy`, `demo_moderate`, `demo_compressed`.
- Complete vignette with worked examples and GIS export.
- Improved DESCRIPTION with full package documentation.

## Bug fixes

- Fixed S3 method dispatch for `summary.sef_fit()` and `print.sef_fit()` (NAMESPACE now uses `S3method()` instead of `export()`).
- Fixed `as_sf_phase()` and `as_sf_links()` CRS default (`NA` → `NA_integer_`) for compatibility with sf >= 1.0.
- Fixed startup message to reflect actual package version.
- Corrected `Authors@R` in DESCRIPTION.

## Internal

- Removed unused `LazyData` field from DESCRIPTION.
- Cleaned up Rd documentation.

# palimpsestr 0.5.0

- Added diagonal Gaussian-mixture EM refinement in `fit_sef()`.
- Added likelihood and BIC diagnostics.
- Added stronger Harris validation and retained penalty handling.
- Added `as_sf_links()` and enhanced `as_sf_phase()` with XY/XYZ support.
