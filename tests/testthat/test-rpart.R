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

test_that("fit_rpart(prune =) prunes at the cross-validated cp", {
  skip_if_not_installed("rpart")
  d <- survival::veteran; f <- Surv(time, status) ~ age + karno + diagtime + celltype
  set.seed(1); a <- fit_rpart(f, d, cp = 0.001)
  set.seed(1); b <- fit_rpart(f, d, cp = 0.001, prune = "min")
  set.seed(1); c1 <- fit_rpart(f, d, cp = 0.001, prune = "1se")
  ct <- a$model$cptable; i <- which.min(ct[, "xerror"])
  expect_equal(nrow(b$model$frame), nrow(rpart::prune(a$model, cp = ct[i, "CP"])$frame))
  expect_lte(nrow(c1$model$frame), nrow(b$model$frame))
  expect_lt(nrow(b$model$frame), nrow(a$model$frame))
  set.seed(1); expect_equal(fit_rpart(f, d)$model$frame, fit(f, d, "rpart", seed = 1)$model$frame)
})
