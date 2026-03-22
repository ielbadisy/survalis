test_that("predict_rpart() returns bounded survival with standard column names", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("rpart")
  testthat::skip_if_not_installed("partykit")

  df <- survival::veteran
  mod <- fit_rpart(survival::Surv(time, status) ~ age + karno + celltype, data = df)

  times <- c(200, 50, 300)
  pred <- predict_rpart(mod, newdata = df[1:6, ], times = times)

  expect_s3_class(pred, "data.frame")
  expect_identical(colnames(pred), paste0("t=", times))
  expect_equal(nrow(pred), 6)
  expect_equal(ncol(pred), length(times))

  M <- as.matrix(pred)
  expect_true(all(is.finite(M)))
  expect_true(all(M >= 0 & M <= 1))

  ord <- order(times)
  dec_ok <- apply(M[, ord, drop = FALSE], 1, function(r) all(diff(r) <= 1e-8))
  expect_true(all(dec_ok))
})

test_that("tune_rpart() returns a fitted model when refit_best is TRUE", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("rpart")
  testthat::skip_if_not_installed("partykit")
  testthat::skip_if_not_installed("rsample")

  df <- survival::veteran

  mod <- tune_rpart(
    formula = survival::Surv(time, status) ~ age + karno + celltype,
    data = df,
    times = c(50, 100),
    param_grid = expand.grid(
      minsplit = c(10, 20),
      cp = c(0.001),
      maxdepth = c(5)
    ),
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 123,
    refit_best = TRUE
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "rpart")
  expect_true(!is.null(attr(mod, "tuning_results")))
})
