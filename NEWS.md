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
