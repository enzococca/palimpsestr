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

test_that("export_sef_report errors clearly without rmarkdown is not triggered here", {
  # guard: the function must validate its input type
  expect_error(export_sef_report(list(), tempfile()), "sef_fit")
})
