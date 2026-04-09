test_that("chrono_precision adds a column", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  feat_base <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class")
  feat_cp <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                             add_chrono_precision = TRUE)
  expect_equal(ncol(feat_cp), ncol(feat_base) + 1)
})

test_that("residuality adds a column when context present", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  feat_base <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class")
  feat_res <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                              context_col = "context")
  expect_equal(ncol(feat_res), ncol(feat_base) + 1)
})

test_that("residuality skipped when context_col is NULL", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  feat_base <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class")
  feat <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                          context_col = NULL)
  expect_equal(ncol(feat), ncol(feat_base))
})

test_that("class_scale reduces one-hot magnitude", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  feat_raw <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                              class_scale = FALSE)
  feat_sc <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                             class_scale = TRUE)
  n_num <- 5
  class_raw <- feat_raw[, (n_num + 1):ncol(feat_raw)]
  class_sc <- feat_sc[, (n_num + 1):ncol(feat_sc)]
  expect_true(max(abs(class_sc)) < max(abs(class_raw)))
})

test_that("taf_as_feature adds a column", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  feat_base <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class")
  feat_taf <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                              add_taf = TRUE, taf_col = "taf_score")
  expect_equal(ncol(feat_taf), ncol(feat_base) + 1)
})

test_that("subclass overrides class for encoding", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  x$subclass <- paste0(x$class, "_sub", sample(1:2, 40, replace = TRUE))
  feat_class <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class")
  feat_sub <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                              subclass_col = "subclass")
  expect_true(ncol(feat_sub) > ncol(feat_class))
})

test_that("all 5 improvements together in feature_matrix", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  x$subclass <- paste0(x$class, "_", sample(letters[1:2], 40, replace = TRUE))
  feat <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class",
                          add_chrono_precision = TRUE, add_taf = TRUE,
                          taf_col = "taf_score", context_col = "context",
                          class_scale = TRUE, subclass_col = "subclass")
  feat_base <- palimpsestr:::feature_matrix(x, c("x","y","z"), c("date_min","date_max"), "class")
  # Base: 5 numeric + n_class one-hot. New: 5 + 3 (prec,resid,taf) + n_subclass
  expect_true(ncol(feat) > ncol(feat_base))
})

test_that("fit_sef backward compatible with no new params", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2)
  expect_s3_class(fit, "sef_fit")
  expect_false(fit$chrono_precision)
  expect_false(fit$taf_as_feature)
  expect_false(fit$residuality)
  expect_false(fit$class_scale)
  expect_null(fit$subclass)
})

test_that("fit_sef works with all 5 improvements", {
  x <- archaeo_sim(n = 60, k = 2, seed = 1)
  x$subclass <- paste0(x$class, "_", sample(letters[1:3], 60, replace = TRUE))
  fit <- fit_sef(x, k = 2, tafonomy = "taf_score", context = "context",
                 chrono_precision = TRUE, taf_as_feature = TRUE,
                 residuality = TRUE, class_scale = TRUE,
                 subclass = "subclass")
  expect_s3_class(fit, "sef_fit")
  expect_true(fit$chrono_precision)
  expect_true(fit$taf_as_feature)
  expect_true(fit$residuality)
  expect_true(fit$class_scale)
  expect_equal(fit$subclass, "subclass")
})

test_that("bootstrap_sef preserves new params", {
  x <- archaeo_sim(n = 40, k = 2, seed = 1)
  fit <- fit_sef(x, k = 2, tafonomy = "taf_score",
                 chrono_precision = TRUE, taf_as_feature = TRUE)
  bs <- bootstrap_sef(fit, n_boot = 3, verbose = FALSE)
  expect_true(is.data.frame(bs))
})
