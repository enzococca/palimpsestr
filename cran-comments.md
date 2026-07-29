## Submission

This is a fix for the check ERROR on the r-devel Debian flavours of 0.24.0.

`export_sef_report()` rendered its bundled RMarkdown template in place, and
`rmarkdown::render()` writes its intermediate files (`*.knit.md`, figure
directories) next to the input document. The installed package directory was
therefore written to during rendering, which violates the CRAN Policy rule that
packages must not write anywhere outside the R session's temporary directory,
and failed outright now that the check library is mounted read-only.

The template is now copied into a per-call temporary directory and rendered
there; nothing is written outside `tempdir()`. This was verified by installing
the package into a library made read-only (`chmod -R a-w`) and running the full
test suite against it: 0 failures. Two regression tests were added to
`tests/testthat/test-export-report.R`, asserting that the staged template lies
under `tempdir()` and outside `system.file(package = "palimpsestr")`, and that
`rmarkdown::render()` is called on the staged copy rather than on the installed
file.

No other code in the package writes outside `tempdir()`.

**On the version jump (0.10.0 -> 0.24.x).** Since the 0.10.0 CRAN release, the
package has undergone substantial methodological development in its public
repository, and the manuscript describing it has been peer-reviewed and
recommended by Peer Community in Archaeology
(https://doi.org/10.24072/pci.archaeo.101019). The intermediate versions
(0.11 through 0.23) were distributed only via GitHub and archived on Zenodo;
they were not submitted to CRAN. All changes are documented in NEWS.md.

## Test environments

* local: macOS (aarch64-apple-darwin), R 4.6.0
* local: macOS with the package library mounted read-only, which reproduces the
  Debian check condition
* TeX Live 2025 (vignette and PDF manual build)

## R CMD check results

0 errors | 0 warnings | 2 notes

* NOTE (CRAN incoming feasibility): "Days since last update: 4". This
  submission is solely a fix for the check ERROR reported on the r-devel
  Debian flavours of 0.24.0; no other change is included.
* NOTE (HTML version of manual): "Skipping checking math rendering: package
  'V8' unavailable" — local only, 'V8' is not installed on this machine.

## Downstream dependencies

There are no reverse dependencies.
