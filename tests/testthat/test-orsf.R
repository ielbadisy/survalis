
test_that("fit_orsf() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("aorsf")
  skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_orsf(
    survival::Surv(time, status) ~ age + karno,
    data = df,
    n_tree = 100, mtry = 2
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "orsf")
  expect_identical(attr(mod, "engine"), "aorsf")
  expect_true(!is.null(mod$model))
})

test_that("predict_orsf() returns bounded, monotone survival at requested times", {
  skip_on_cran()
  skip_if_not_installed("aorsf")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_orsf(
    survival::Surv(time, status) ~ age + karno,
    data = df,
    n_tree = 120, mtry = 2
  )

  times <- c(150, 30, 300, 200)
  newx  <- df[1:6, c("age", "karno")]

  pred <- predict_orsf(mod, newdata = newx, times = times)

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

test_that("predict_orsf() errors when newdata lacks required predictors", {
  skip_on_cran()
  skip_if_not_installed("aorsf")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_orsf(
    survival::Surv(time, status) ~ age + karno,
    data = df, n_tree = 80, mtry = 2
  )

  bad_new <- df[1:3, c("age")] 
  expect_error(
    predict_orsf(mod, newdata = bad_new, times = c(50, 100)),
    regexp = paste0(
      "object .* not found|undefined columns|subset.*select|",
      "do not match|is missing|required .*column|",
      "training data.*have columns.*not contained in.*new_data"
    )
  )
  
})

test_that("tune_orsf() returns a compact grid and refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("aorsf")
  skip_if_not_installed("survival")

  Surv <- survival::Surv

  df <- survival::veteran
  set.seed(2025)

  grid <- expand.grid(
    n_tree     = c(80, 120),   
    mtry       = c(1, 2),
    min_events = c(5, 10),
    stringsAsFactors = FALSE
  )

  res <- tune_orsf(
    formula = Surv(time, status) ~ age + karno,
    data = df,
    times = c(60, 180),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = FALSE
  )

  expect_s3_class(res, "tuned_surv")
  expect_true(all(c("n_tree", "mtry", "min_events") %in% names(res)))
  expect_true(all(c("cindex", "ibs") %in% names(res)))
  expect_true(nrow(res) >= 1)

  best_mod <- tune_orsf(
    formula = Surv(time, status) ~ age + karno,
    data = df,
    times = c(60, 180),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_identical(best_mod$learner, "orsf")
  expect_identical(attr(best_mod, "engine"), "aorsf")
  expect_true(!is.null(best_mod$model))
})