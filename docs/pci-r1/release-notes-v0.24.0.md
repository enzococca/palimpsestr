# palimpsestr 0.24.0

**PCI Archaeo revision (round 1).** This release accompanies the revised manuscript submitted to PCI Archaeology (#1019). It corrects the sign of the ICL model-selection criterion, adds two reviewer-requested refinements to the data-simulation and Harris-matrix helpers, and documents the model-selection metrics in both the help pages and the manuscript. The SEF model itself is unchanged since 0.17.1, so all case-study results reproduce identically.

## ⚠️ Behaviour change

- **ICL sign correction.** The Integrated Completed Likelihood is now computed as `bic + 2 * sum(entropy)` (previously `bic - 2 * sum(entropy)`). With the lower-is-better BIC convention the entropy term must be **added** so that overlapping/fuzzy partitions are *penalised* (Biernacki, Celeux & Govaert 2000); the previous sign rewarded them. **`compare_k()` ICL values change accordingly.** BIC, PDI, and every case-study figure are unaffected (ICL is not used in any manuscript figure), so existing K-selection conclusions stand.

## Bug fixes

- Fixed the ICL sign (see above), with a regression test asserting `icl == bic + 2 * sum(entropy)` and `icl >= bic`.

## Other changes

- **`archaeo_sim()`** now returns the simulated date bounds `date_min` and `date_max` as integer-valued years (via `floor()`), consistent with archaeological dating conventions. *(reviewer suggestion, PCI Archaeo #1019)*
- **`harris_from_contexts()`** gains an `exclude_contexts` argument: a character vector of context labels exempt from the verticality (depth-rank) penalty — for the fills of cuts and multi-period contexts where depth does not track chronology. The default (`NULL`) leaves existing behaviour unchanged. The help page now documents the verticality assumption and its failure modes. *(reviewer suggestion, PCI Archaeo #1019)*
- **Documentation.** The `compare_k()` and `pdi()` help pages now explain the BIC, ICL, and PDI model-selection metrics and reference Schwarz (1978) and Biernacki et al. (2000); the manuscript adds the same definitions and stresses combining statistical criteria with archaeological judgement.

## Installation

```r
# from GitHub
remotes::install_github("enzococca/palimpsestr@v0.24.0")
```

Or install from the source tarball attached to this release. The package is also archived on Zenodo (concept DOI [10.5281/zenodo.19881542](https://doi.org/10.5281/zenodo.19881542)).

## Compatibility

No breaking API changes. The only output difference is the value (not the meaning) of `ICL` returned by `compare_k()`; all other functions, arguments, and return structures are unchanged.

**Full changelog:** see [`NEWS.md`](https://github.com/enzococca/palimpsestr/blob/master/NEWS.md).
