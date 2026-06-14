# Design: OxCal adapter + intrusion typing (v0.19.0)

## Goal

Implement the two "ready now" roadmap items: an OxCal chronology adapter
(`chronology_from_oxcal()`) symmetric to the existing rcarbon adapter, and a
tighter coupling between the noise-component intrusion magnitude and the
chronological direction, surfaced as an `intrusion_type` column in
`detect_intrusions()` and used to colour `gg_outliers()`.

## Scope

- New exported `chronology_from_oxcal()` in `R/chronology.R`.
- `detect_intrusions()` gains an `intrusion_threshold` argument and an
  `intrusion_type` factor column (existing columns unchanged).
- `gg_outliers()` colours points by `intrusion_type`.
- `oxcAAR` added to Suggests.
- Manual + paper updated; version bump to 0.19.0.

Out of scope: running OxCal itself (the user calibrates with oxcAAR or exports
ranges), new diagnostics, changes to the fitting engine.

## Feature 1 — `chronology_from_oxcal()`

Signature mirrors `chronology_from_rcarbon()`:

```r
chronology_from_oxcal(x,
                      method = c("hpd", "median_iqr", "weighted_mean"),
                      hpd = 0.95, ids = NULL, bce_negative = TRUE)
```

Dispatch on the type of `x`:

- **`oxcAARCalibratedDatesList`** (from `oxcAAR::oxcalCalibrate()`): each
  element carries `$raw_probabilities` (a data.frame with `dates` in
  calendar years — BC negative — and `probabilities`). Reduce each date to an
  interval by `method`:
  - `"hpd"`: the highest-posterior-density region at level `hpd`
    (oldest/youngest covered calendar year);
  - `"median_iqr"`: weighted median plus 25/75 percentiles;
  - `"weighted_mean"`: density-weighted mean ± 2 SD.
  Requires the `oxcAAR` package (Suggests); error with an install hint when
  absent.
- **`data.frame`**: must contain numeric `start` and `end` columns (calendar
  years, BC negative — the typical OxCal range export) and an optional `id`
  column. The interval is taken directly; `method` is not applicable and is
  ignored with a message. `date_mid` is the midpoint.
- **anything else**: informative error.

Output (identical to `chronology_from_rcarbon()`): a data.frame with `id`,
`date_min`, `date_max`, `date_mid`. Bounds are normalised ascending with
`pmin`/`pmax`. `bce_negative = TRUE` (default) returns calendar years with BC
negative; `FALSE` returns calBP (`1950 - year`), so the two adapters share one
output convention. `ids` overrides the identifiers; default `cal_1 … cal_n`.

## Feature 2 — `intrusion_type` in `detect_intrusions()`

New argument `intrusion_threshold = 0.5`. New factor column `intrusion_type`
combining the intrusion magnitude with the chronological direction:

| Condition | `intrusion_type` |
|-----------|------------------|
| `intrusion_prob < threshold` | `not_flagged` |
| flagged & `direction == "older_than_context"` | `residual` |
| flagged & `direction == "younger_than_context"` | `latent_feature` |
| flagged & (`in_context` or `direction` is `NA`) | `outlier_in_context` |

Factor levels (fixed order):
`c("not_flagged", "residual", "latent_feature", "outlier_in_context")`.
The existing `id`, `intrusion_prob`, `direction`, `chrono_gap` columns are
unchanged, so the addition is backward compatible.

## Feature 3 — `gg_outliers()` colours by `intrusion_type`

`gg_outliers()` already computes `detect_intrusions(object)`; it now colours
each lollipop by `intrusion_type` instead of the binary above/below flag,
with a fixed colour-blind-safe mapping
(`not_flagged` grey, `residual` blue, `latent_feature` vermillion,
`outlier_in_context` amber) and `drop = FALSE` so the legend is stable. The
`threshold` reference line and the rest of the layout are unchanged. The
`intrusion_threshold` is passed through so the colouring matches the table.

## Edge cases

- OxCal: a date whose HPD region is disjoint — take the outer covered bounds
  (oldest start, youngest end), mirroring the rcarbon adapter.
- `data.frame` input with `start > end` per row — normalised with pmin/pmax.
- `detect_intrusions()` on a fit with no direction (no context / unit-tied
  dating): flagged finds become `outlier_in_context`; unflagged `not_flagged`.
- `gg_outliers()` when only one `intrusion_type` is present: legend still lists
  the fixed levels (drop = FALSE).

## Testing

- `tests/testthat/test-chronology-oxcal.R`: data.frame path (start/end →
  date_min/date_max, ascending, midpoint); `ids` honoured; `bce_negative`
  conversion; error on unsupported input; oxcAAR-object path skipped when the
  package is unavailable (`skip_if_not_installed("oxcAAR")`).
- `tests/testthat/test-intrusion-type.R`: each `intrusion_type` value arises
  from the right magnitude/direction combination; column is a factor with the
  fixed levels; backward-compatible columns preserved; threshold honoured.
- `tests/testthat/test-gg-model-plots.R` (extend): `gg_outliers()` still
  returns a ggplot and carries an `intrusion_type` mapping.

## Documentation

- Manual: add `chronology_from_oxcal()` to the chronological-import material and
  `intrusion_type` to the intrusion-detection section.
- Paper: move "OxCal integration" and "tighter coupling between the noise
  component and the directional intrusion diagnostics" from Future Development
  into the list of completed methodological developments.
- NEWS entry under `# palimpsestr 0.19.0`; regenerate `man/*.Rd` + NAMESPACE.

## Deliverable / version

`v0.19.0`. Acceptance: full suite green, `R CMD check` 0/0/0 under R 4.6, both
functions exported and documented, manual + paper updated.
