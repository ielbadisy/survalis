
test_that("fit_survsvm() returns a proper mlsurv_model wrapper", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_survsvm(
    survival::Surv(time, status) ~ age + celltype + karno,
    data = df,
    type = "regression",
    gamma.mu = 0.1,
    kernel = "lin_kernel"
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "survsvm")
  expect_identical(attr(mod, "engine"), "survivalsvm")
  expect_true(!is.null(mod$model))
})

test_that("predict_survsvm() returns well-formed survival probabilities (exp + weibull)", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("survival")

  set.seed(123)
  df <- survival::veteran
  mod <- fit_survsvm(
    survival::Surv(time, status) ~ age + celltype + karno,
    data = df,
    type = "regression",
    gamma.mu = 0.1,
    kernel = "lin_kernel"
  )

  newx  <- df[1:5, c("age", "celltype", "karno")]
  times <- c(50, 150, 300)

  pred_exp <- predict_survsvm(mod, newdata = newx, times = times, dist = "exp")
  expect_s3_class(pred_exp, "data.frame")
  expect_equal(nrow(pred_exp), nrow(newx))
  expect_setequal(colnames(pred_exp), paste0("t=", times))
  M <- as.matrix(pred_exp)
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)
  ord <- order(times)
  diffs <- t(apply(M[, ord, drop = FALSE], 1, diff))
  expect_true(all(diffs <= 1e-8))

  pred_wb <- predict_survsvm(mod, newdata = newx, times = times, dist = "weibull", shape = 1.3)
  expect_s3_class(pred_wb, "data.frame")
  expect_equal(nrow(pred_wb), nrow(newx))
  expect_setequal(colnames(pred_wb), paste0("t=", times))
  Mw <- as.matrix(pred_wb)
  expect_true(all(is.finite(Mw)))
  expect_gte(min(Mw), -1e-12)
  expect_lte(max(Mw), 1 + 1e-12)
  diffs_w <- t(apply(Mw[, ord, drop = FALSE], 1, diff))
  expect_true(all(diffs_w <= 1e-8))
})

test_that("predict_survsvm() errors when newdata lacks required predictors", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("survival")

  df <- survival::veteran
  mod <- fit_survsvm(
    survival::Surv(time, status) ~ age + celltype + karno,
    data = df, type = "regression", gamma.mu = 0.1, kernel = "lin_kernel"
  )

  bad_new <- df[1:3, c("age", "karno")] 
  expect_error(
    predict_survsvm(mod, newdata = bad_new, times = c(50, 100)),
    regexp = "not.*found|do not match|undefined|independ|covariat|column"
  )
})

test_that("tune_survsvm() returns a compact grid; refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("survival")

  Surv <- survival::Surv

  df <- survival::veteran
  set.seed(777)

  grid <- list(
    gamma.mu = c(0.05, 0.1),              
    kernel   = c("lin_kernel", "add_kernel")
  )

  res <- tune_survsvm(
    formula = Surv(time, status) ~ age + celltype + karno,
    data = df,
    times = c(100, 300),
    metrics = c("cindex", "ibs"),
    param_grid = grid,
    folds = 2,
    seed = 777,
    refit_best = FALSE
  )

  expect_s3_class(res, "tune_surv")
  expect_true(all(c("gamma.mu", "kernel") %in% names(res)))
  expect_true(all(c("cindex", "ibs") %in% names(res)))
  expect_true(nrow(res) >= 1)

  best_mod <- tune_survsvm(
    formula = Surv(time, status) ~ age + celltype + karno,
    data = df,
    times = c(100, 300),
    metrics = c("cindex", "ibs"),
    param_grid = grid,
    folds = 2,
    seed = 777,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_true(inherits(best_mod, "tune_surv"))
  expect_identical(best_mod$learner, "survsvm")
  expect_identical(attr(best_mod, "engine"), "survivalsvm")
  expect_true(!is.null(attr(best_mod, "tuning_results")))
})
