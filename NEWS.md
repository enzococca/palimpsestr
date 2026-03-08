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
