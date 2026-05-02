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
