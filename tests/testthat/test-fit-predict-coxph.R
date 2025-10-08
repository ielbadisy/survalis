test_that("fit_coxph() fits a Cox model and returns expected structure", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")

  df <- survival::veteran

  mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = df)

  expect_true(is.list(mod))
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "coxph")
  expect_true(!is.null(mod$formula))
  expect_true(!is.null(mod$data))
  expect_true(identical(mod$time, "time"))
  expect_true(identical(mod$status, "status"))

  expect_true(!is.null(mod$model))
  expect_s3_class(mod$model, "coxph")
})

test_that("predict_coxph() returns survival probabilities with correct shape and names", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")

  set.seed(123)
  df <- survival::veteran
  n <- nrow(df)
  idx <- sample(seq_len(n), size = floor(0.7 * n))
  train <- df[idx, ]
  test  <- df[-idx, ]

  mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = train)

  times <- c(100, 200, 300)
  pred <- predict_coxph(mod, newdata = test, times = times)

  expect_true(is.data.frame(pred))
  expect_equal(nrow(pred), nrow(test))
  expect_equal(ncol(pred), length(times))
  expect_identical(colnames(pred), paste0("t=", times))
  expect_true(all(is.finite(as.matrix(pred))))
  expect_true(all(as.matrix(pred) >= 0 & as.matrix(pred) <= 1))
  dec_ok <- apply(as.matrix(pred), 1, function(r) all(diff(r) <= 1e-8))
  expect_true(all(dec_ok))
})

test_that("predict_coxph() errors on missing predictors in newdata", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")

  df <- survival::veteran
  mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = df)

  bad_new <- subset(df, select = -age)
  expect_error(
    predict_coxph(mod, newdata = bad_new, times = c(100, 200)),
    regexp = "age|column|predict"
  )
})

test_that("predict_coxph() warns if learner field is not 'coxph'", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")

  df <- survival::veteran
  mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = df)

  mod_bad <- mod
  mod_bad$learner <- "not_coxph"

  expect_warning(
    predict_coxph(mod_bad, newdata = df[1:5, ], times = c(100, 200)),
    regexp = "may not come from fit_coxph"
  )
})
