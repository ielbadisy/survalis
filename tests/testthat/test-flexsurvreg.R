
test_that("fit_flexsurvreg() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_flexsurvreg(
    survival::Surv(time, status) ~ age + celltype + karno,
    data = df,
    dist = "weibull"
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "flexsurvreg")
  expect_identical(attr(mod, "engine"), "flexsurv")
  expect_true(!is.null(mod$model))
})

test_that("predict_flexsurvreg() returns bounded, monotone survival at requested times", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_flexsurvreg(
    survival::Surv(time, status) ~ age + celltype + karno,
    data = df, dist = "weibull"
  )
  times <- c(150, 30, 300, 200)
  newx  <- df[1:6, ]

  pred <- predict_flexsurvreg(mod, newdata = newx, times = times)

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

test_that("predict_flexsurvreg() errors when newdata lacks required predictors", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_flexsurvreg(
    survival::Surv(time, status) ~ age + celltype + karno,
    data = df, dist = "weibull"
  )

  bad_new <- df[1:3, c("age", "karno")]  
  expect_error(
    predict_flexsurvreg(mod, newdata = bad_new, times = c(50, 100)),
    regexp = paste0(
      "object .* not found|undefined columns|subset.*select|",
      "do not match|variable lengths differ|newdata.*missing|terms|",
      'Value of covariate ".*" not supplied in "newdata"'
    )
  )  
})

test_that("tune_flexsurvreg() returns a tuned grid and refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("flexsurv")
  skip_if_not_installed("survival")

  Surv <- survival::Surv

  df <- survival::veteran
  set.seed(2025)

  grid <- c("weibull", "exponential")  

  res <- tune_flexsurvreg(
    formula = Surv(time, status) ~ age + celltype + karno,
    data = df,
    times = c(60, 180, 300),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = FALSE
  )

  expect_s3_class(res, "tuned_surv")
  expect_true("dist" %in% names(res))
  expect_true(all(c("cindex", "ibs") %in% names(res)))
  expect_true(nrow(res) >= 1)

  best_mod <- tune_flexsurvreg(
    formula = Surv(time, status) ~ age + celltype + karno,
    data = df,
    times = c(60, 180, 300),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_identical(best_mod$learner, "flexsurvreg")
  expect_identical(attr(best_mod, "engine"), "flexsurv")
  expect_true(!is.null(best_mod$model))
})