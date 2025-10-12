
test_that("fit_bnnsurv() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_bnnsurv(
    survival::Surv(time, status) ~ age + karno + diagtime + prior,
    data = df,
    k = 3, num_base_learners = 5
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "bnnsurv")
  expect_identical(mod$engine, "bnnSurvival")
  expect_true(!is.null(mod$model))
})

test_that("predict_bnnsurv() (native grid) returns bounded, monotone survival", {
  skip_on_cran()
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_bnnsurv(
    survival::Surv(time, status) ~ age + karno + diagtime + prior,
    data = df, k = 3, num_base_learners = 5
  )

  newx <- df[1:6, ]
  pred <- predict_bnnsurv(mod, newdata = newx) 

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(newx))
  expect_true(ncol(pred) >= 2)               
  expect_true(all(startsWith(colnames(pred), "t=")))

  M <- as.matrix(pred)
  expect_true(is.numeric(M))
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)

  tg <- as.numeric(sub("^t=", "", colnames(pred)))
  ord <- order(tg)
  M_ord <- M[, ord, drop = FALSE]
  diffs <- t(apply(M_ord, 1, diff))
  expect_true(all(diffs <= 1e-8))
})

test_that("predict_bnnsurv() with requested times is shaped right and monotone", {
  skip_on_cran()
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_bnnsurv(
    survival::Surv(time, status) ~ age + karno + diagtime + prior,
    data = df, k = 3, num_base_learners = 5
  )

  times <- c(100, 25, 200)  
  newx  <- df[1:5, ]
  pred  <- predict_bnnsurv(mod, newdata = newx, times = times)

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(newx))
  expect_equal(ncol(pred), length(times))
  expect_setequal(colnames(pred), paste0("t=", as.character(times)))

  M <- as.matrix(pred)
  expect_true(is.numeric(M))
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)

  ord <- order(times)
  M_ord <- M[, paste0("t=", as.character(times[ord])), drop = FALSE]
  diffs <- t(apply(M_ord, 1, diff))
  expect_true(all(diffs <= 1e-8))
})

test_that("predict_bnnsurv() errors when newdata lacks required predictors", {
  skip_on_cran()
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_bnnsurv(
    survival::Surv(time, status) ~ age + karno + diagtime + prior,
    data = df, k = 2, num_base_learners = 5
  )

  bad_new <- df[1:3, c("age", "karno", "diagtime")] 
  expect_error(
    predict_bnnsurv(mod, newdata = bad_new, times = c(50, 100)),
    regexp = paste0(
      "object .* not found|undefined columns|subset.*select|",
      "do not match"  
    )
  )
})

test_that("tune_bnnsurv() returns a compact grid and refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survival")

  Surv <- survival::Surv

  df <- survival::veteran
  set.seed(2025)

  grid <- expand.grid(
    k = c(2, 3),
    num_base_learners = c(5, 8),  
    sample_fraction = c(0.7, 1),
    stringsAsFactors = FALSE
  )

  res <- tune_bnnsurv(
    formula = Surv(time, status) ~ age + karno + diagtime + prior,
    data = df,
    times = c(50, 100),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("k", "num_base_learners", "sample_fraction") %in% names(res)))
  expect_true(all(c("cindex", "ibs") %in% names(res)))
  expect_true(nrow(res) >= 1)

  best_mod <- tune_bnnsurv(
    formula = Surv(time, status) ~ age + karno + diagtime + prior,
    data = df,
    times = c(50, 100),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_identical(best_mod$learner, "bnnsurv")
  expect_identical(best_mod$engine, "bnnSurvival")
  expect_true(!is.null(best_mod$model))
})
