test_that("roc_survmat's trapezoidal AUC closely tracks auc_survmat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  y <- survival::Surv(time = veteran$time, event = veteran$status)
  sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1, dimnames = list(NULL, "t=100"))

  r <- roc_survmat(y, predicted = sp, t_star = 100)
  a <- unname(auc_survmat(y, predicted = sp, t_star = 100))

  expect_s3_class(r, "roc_survmat")
  expect_true(all(c("threshold", "sensitivity", "specificity") %in% colnames(r$curve)))
  expect_equal(r$auc, a, tolerance = 0.02)
})

test_that("plot.roc_survmat / plot_roc return ggplots", {
  testthat::skip_on_cran()

  y <- survival::Surv(time = veteran$time, event = veteran$status)
  sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1, dimnames = list(NULL, "t=100"))
  r <- roc_survmat(y, predicted = sp, t_star = 100)

  expect_s3_class(plot(r), "ggplot")
  expect_s3_class(plot_roc(r), "ggplot")
})

test_that("dca_survmat: treat-none is 0, and model/treat-all are finite over thresholds", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  y <- survival::Surv(time = veteran$time, event = veteran$status)
  sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1, dimnames = list(NULL, "t=100"))

  d <- dca_survmat(y, predicted = sp, t_star = 100)
  expect_s3_class(d, "dca_survmat")
  expect_true(all(d$curve$treat_none == 0))
  expect_true(all(is.finite(d$curve$model)))
  expect_true(all(is.finite(d$curve$treat_all)))

  expect_error(dca_survmat(y, predicted = sp, t_star = 100, thresholds = c(0, 0.5)), "strictly between")
})

test_that("plot.dca_survmat / plot_dca return ggplots", {
  testthat::skip_on_cran()

  y <- survival::Surv(time = veteran$time, event = veteran$status)
  sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1, dimnames = list(NULL, "t=100"))
  d <- dca_survmat(y, predicted = sp, t_star = 100)

  expect_s3_class(plot(d), "ggplot")
  expect_s3_class(plot_dca(d), "ggplot")
})
