# TDD for export_sef_report(): a narrated PDF/DOCX report assembled from an
# RMarkdown template, with a pandoc-free fallback bundle (markdown + PNG figures).

.small_fit <- function() {
  suppressMessages(fit_sef(archaeo_sim(n = 45, k = 3, seed = 1), k = 3,
                           context = "context"))
}

test_that("the no-pandoc fallback bundle writes a markdown narrative + PNG figures", {
  skip_if_not_installed("ggplot2")
  fit <- .small_fit()
  out <- tempfile()
  files <- suppressWarnings(suppressMessages(
    palimpsestr:::.report_fallback_bundle(fit, out, lang = "it",
                                          title = "Test", site = NULL,
                                          intrusion_threshold = 0.5)))
  expect_true(any(grepl("\\.md$", files)))
  expect_true(any(grepl("\\.png$", files)))
  expect_true(all(file.exists(files)))
  md <- files[grepl("\\.md$", files)][1]
  expect_gt(file.info(md)$size, 0)
})

test_that("export_sef_report always writes a markdown narrative sidecar", {
  skip_if_not_installed("rmarkdown")
  fit <- .small_fit()
  out <- tempfile()
  files <- suppressWarnings(suppressMessages(
    export_sef_report(fit, out, format = "docx", lang = "it")))
  md <- files[grepl("\\.md$", files)]
  expect_length(md, 1)
  expect_true(file.exists(md))
  expect_gt(file.info(md)$size, 0)
})

test_that("export_sef_report produces a non-empty DOCX when pandoc is available", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  fit <- .small_fit()
  out <- tempfile()
  files <- suppressWarnings(suppressMessages(
    export_sef_report(fit, out, format = "docx", lang = "it")))
  docx <- files[grepl("\\.docx$", files)]
  expect_length(docx, 1)
  expect_gt(file.info(docx)$size, 0)
})

test_that("export_sef_report narrative differs between Italian and English", {
  skip_if_not_installed("rmarkdown")
  fit <- .small_fit()
  it <- suppressWarnings(suppressMessages(
    export_sef_report(fit, tempfile(), format = "docx", lang = "it")))
  en <- suppressWarnings(suppressMessages(
    export_sef_report(fit, tempfile(), format = "docx", lang = "en")))
  it_md <- readLines(it[grepl("\\.md$", it)][1], warn = FALSE)
  en_md <- readLines(en[grepl("\\.md$", en)][1], warn = FALSE)
  expect_false(identical(it_md, en_md))
})

test_that("the report template is staged in a temporary directory, never rendered in place", {
  skip_if_not_installed("rmarkdown")
  staged <- palimpsestr:::.stage_report_template()
  on.exit(unlink(dirname(staged), recursive = TRUE), add = TRUE)
  expect_true(file.exists(staged))
  # rmarkdown::render() writes its intermediates next to the input file, so the
  # staged copy must live under tempdir() and outside the installed package
  # (the library is read-only on the CRAN check machines).
  expect_true(startsWith(normalizePath(staged), normalizePath(tempdir())))
  expect_false(startsWith(normalizePath(staged),
                          normalizePath(system.file(package = "palimpsestr"))))
})

test_that("export_sef_report renders the staged copy, not the installed template", {
  skip_if_not_installed("rmarkdown")
  pkg <- normalizePath(system.file(package = "palimpsestr"))
  seen <- NULL
  testthat::local_mocked_bindings(
    render = function(input, ..., output_file, output_dir) {
      seen <<- normalizePath(input)
      target <- file.path(output_dir, output_file)
      writeLines("stub", target)
      target
    },
    .package = "rmarkdown"
  )
  fit <- .small_fit()
  suppressWarnings(suppressMessages(
    export_sef_report(fit, tempfile(), format = "docx", lang = "it")))
  expect_true(startsWith(seen, normalizePath(tempdir())))
  expect_false(startsWith(seen, pkg))
})

test_that("export_sef_report errors clearly without rmarkdown is not triggered here", {
  # guard: the function must validate its input type
  expect_error(export_sef_report(list(), tempfile()), "sef_fit")
})
