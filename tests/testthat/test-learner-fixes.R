# Regression tests for learner bugs found while benchmarking:
# aalen scaling, single-time predictions, survmetalearner via formula and data,
# xgboost tuning parameters, tuner default grids, and errors hidden by parallel workers.

lung_complete <- function() {
  d <- survival::lung
  d$status <- as.integer(d$status == 2)
  d <- d[, c("time", "status", "age", "sex", "ph.ecog", "ph.karno", "pat.karno", "meal.cal", "wt.loss")]
  stats::na.omit(d)
}
lung_form <- survival::Surv(time, status) ~ age + sex + ph.ecog + ph.karno + pat.karno + meal.cal + wt.loss
vet_form <- survival::Surv(time, status) ~ karno + age + celltype

test_that("aalen predicts varying survival when covariates are on very different scales", {
  skip_on_cran()
  skip_if_not_installed("timereg")
  d <- lung_complete()
  n_train <- 120
  tm <- default_times(d$time, d$status, n = 3L, range = c(0.25, 0.75))
  mod <- fit_aalen(lung_form, d[seq_len(n_train), ])
  pred <- as.matrix(predict_aalen(mod, d[-seq_len(n_train), ], tm))
  expect_gt(min(apply(pred, 2, stats::sd)), 0.01)
  expect_lt(min(pred), 0.99)
  # the standardization is applied at predict time: prediction uses the stored scaling
  expect_true(all(c("center", "scale") %in% names(mod$scaling)))
  expect_true(all(c("meal.cal", "age") %in% names(mod$scaling$center)))
})

test_that("predictions accept a single evaluation time", {
  skip_on_cran()
  d <- survival::veteran
  tr <- d[1:90, ]
  te <- d[91:130, ]
  t1 <- stats::median(d$time)
  for (id in c("flexsurvreg", "bnnsurv", "cforest")) {
    pkg <- c(flexsurvreg = "flexsurv", bnnsurv = "bnnSurvival", cforest = "party")[[id]]
    skip_if_not_installed(pkg)
    fit_fun <- get(paste0("fit_", id))
    pred_fun <- get(paste0("predict_", id))
    mod <- fit_fun(vet_form, tr)
    p1 <- as.matrix(pred_fun(mod, newdata = te, times = t1))
    expect_equal(dim(p1), c(nrow(te), 1L), info = id)
    expect_false(anyNA(p1), info = id)
    # the single-time column equals the corresponding column of a multi-time prediction
    pm <- as.matrix(pred_fun(mod, newdata = te, times = c(t1 / 2, t1, 2 * t1)))
    expect_equal(unname(p1[, 1]), unname(pm[, 2]), tolerance = 1e-6, info = id)
  }
})

test_that(".rows_by_times() orients one-time and many-time results", {
  expect_equal(dim(survalis:::.rows_by_times(c(0.9, 0.8, 0.7))), c(3L, 1L))
  expect_equal(dim(survalis:::.rows_by_times(matrix(1:6, nrow = 2))), c(3L, 2L))
})

test_that("survmetalearner runs from a formula and data", {
  skip_on_cran()
  skip_if_not_installed("nnls")
  skip_if_not_installed("randomForestSRC")
  d <- survival::veteran
  d$trt <- factor(d$trt)
  mod <- fit_survmetalearner(formula = vet_form, data = d, learners = c("coxph", "glmnet"), folds = 3)
  expect_s3_class(mod, "survmetalearner")
  expect_setequal(mod$learners, c("coxph", "glmnet"))
  expect_equal(unname(colSums(mod$weights)), rep(1, ncol(mod$weights)), tolerance = 1e-6)
  # times that differ from the fitted grid are handled by interpolating the weights
  tm <- c(60, 120, 250)
  p <- as.matrix(predict_survmetalearner(mod, d[1:8, ], tm))
  expect_equal(dim(p), c(8L, 3L))
  expect_false(anyNA(p))
  p1 <- as.matrix(predict_survmetalearner(mod, d[1:8, ], 120))
  expect_equal(dim(p1), c(8L, 1L))
  # and it runs through benchmark()
  res <- benchmark(formula = vet_form, data = d, learners = "survmetalearner", times = tm,
                   metrics = c("cindex", "ibs"), tune = FALSE, folds = 3, ncores = 1,
                   suppress_errors = FALSE)
  expect_gt(nrow(res), 0L)
  expect_false(anyNA(res$value))
})

test_that("xgboost learner accepts the parameters the tuner passes", {
  skip_on_cran()
  skip_if_not_installed("xgboost")
  d <- survival::veteran
  mod <- fit_xgboost(vet_form, d, nrounds = 20, max_depth = 2, eta = 0.1)
  expect_s3_class(mod, "mlsurv_model")
  grid <- expand.grid(nrounds = 20, max_depth = c(2, 3), eta = 0.1,
                      aft_loss_distribution = "extreme", aft_loss_distribution_scale = 1,
                      objective = "survival:aft", stringsAsFactors = FALSE)
  tm <- default_times(d$time, d$status, n = 3L, range = c(0.25, 0.75))
  res <- tune_xgboost(vet_form, d, tm, param_grid = grid, metrics = "ibs", folds = 3, seed = 1)
  expect_false(any(as.logical(res$failed)))
  expect_true("ibs" %in% names(res))
})

test_that("tune() reports a clear error when every configuration fails", {
  skip_on_cran()
  skip_if_not_installed("xgboost")
  d <- survival::veteran
  # an unsupported objective makes fit_xgboost() fail for every configuration
  bad_grid <- expand.grid(nrounds = 10, max_depth = 2, eta = 0.1,
                          aft_loss_distribution = "extreme", aft_loss_distribution_scale = 1,
                          objective = "not_an_objective", stringsAsFactors = FALSE)
  expect_error(
    tune(vet_form, d, "xgboost", grid = bad_grid, resampling = cv(v = 3, seed = 1), metric = "ibs"),
    "All 1 tuning configurations of 'xgboost' failed"
  )
})

test_that("tuners for survsvm and survdnn have default grids", {
  expect_true(is.list(formals(tune_survsvm)$param_grid) || is.call(formals(tune_survsvm)$param_grid))
  expect_false(identical(formals(tune_survsvm)$param_grid, quote(expr = )))
  expect_false(identical(formals(tune_survdnn)$param_grid, quote(expr = )))
  # the argument check that used to reject the default grid no longer fires
  err <- tryCatch(
    tune_survdnn(vet_form, survival::veteran[1:3, ], times = c(50, 100), folds = 3, param_grid = list(
      hidden = list(c(4)), lr = 0.01, activation = "relu", epochs = 1, loss = "cox",
      optimizer = "adam", dropout = 0, batch_norm = TRUE)),
    error = function(e) conditionMessage(e)
  )
  expect_false(is.character(err) && grepl("!missing", err, fixed = TRUE))
})

test_that("compare() reports the real error when a metric fails on parallel workers", {
  skip_on_cran()
  skip_on_os("windows")
  d <- survival::veteran
  d$trt <- factor(d$trt)
  tm <- c(23.5, 62, 145.75)
  msg <- function(nc) tryCatch(
    compare(vet_form, d, models = c("coxph", "rsf"), times = tm, resampling = cv(v = 3, seed = 1),
            metrics = c("cindex", "brier"), ncores = nc),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg(1), "single time point")
  expect_match(msg(2), "single time point")
  expect_false(grepl("object 'value' not found", msg(2), fixed = TRUE))
})

test_that("glmnet accepts a dataset with a single predictor", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  d <- survival::veteran
  f1 <- survival::Surv(time, status) ~ karno
  mod <- fit_glmnet(f1, d)
  tm <- default_times(d$time, d$status, n = 3L, range = c(0.25, 0.75))
  p <- as.matrix(predict_glmnet(mod, d[1:30, ], tm))
  expect_equal(dim(p), c(30L, 3L))
  expect_false(anyNA(p))
  expect_gt(max(apply(p, 2, stats::sd)), 0.01)
})
