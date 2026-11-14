
test_that("fit_ranger() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_ranger(
    survival::Surv(time, status) ~ age + karno + celltype,
    data = df,
    num.trees = 100, mtry = 2, min.node.size = 3
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "ranger")
  expect_identical(attr(mod, "engine"), "ranger")
  expect_true(!is.null(mod$model))
})

test_that("predict_ranger() returns bounded, monotone survival at requested times", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_ranger(
    survival::Surv(time, status) ~ age + karno + celltype,
    data = df,
    num.trees = 120, mtry = 2, min.node.size = 3
  )

  times <- c(150, 30, 300, 200)
  newx  <- df[1:6, ]

  pred <- predict_ranger(mod, newdata = newx, times = times)

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

test_that("predict_ranger() errors when newdata lacks required predictors", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_ranger(
    survival::Surv(time, status) ~ age + karno + celltype,
    data = df,
    num.trees = 80, mtry = 2, min.node.size = 3
  )

  bad_new <- df[1:3, c("age", "karno")]  
  expect_error(
    predict_ranger(mod, newdata = bad_new, times = c(50, 100)),
    regexp = paste0(
      "object .* not found|undefined columns|subset.*select|",
      "do not match|newdata.*missing|terms|model\\.matrix|",
      "No such column|not contained in.*data|variable lengths differ|",
      "One or more independent variables not found in data"
    )
  )
  
})

test_that("tune_ranger() returns a compact grid and refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  skip_if_not_installed("survival")

  Surv <- survival::Surv

  df <- survival::veteran
  set.seed(2025)

  grid <- expand.grid(
    num.trees = c(80, 120),  
    mtry = c(1, 2),
    min.node.size = c(3, 5),
    stringsAsFactors = FALSE
  )

  res <- tune_ranger(
    formula = Surv(time, status) ~ age + karno + celltype,
    data = df,
    times = c(60, 180),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("num.trees", "mtry", "min.node.size") %in% names(res)))
  expect_true(all(c("cindex", "ibs") %in% names(res)))
  expect_true(nrow(res) >= 1)

  best_mod <- tune_ranger(
    formula = Surv(time, status) ~ age + karno + celltype,
    data = df,
    times = c(60, 180),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 2025,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_true(inherits(best_mod, "tune_surv"))
  expect_identical(best_mod$learner, "ranger")
  expect_identical(attr(best_mod, "engine"), "ranger")
  expect_true(!is.null(best_mod$model))
})
test_that("fit_ranger() passes case weights to the engine", {
  testthat::skip_if_not_installed("ranger")
  lung <- survival::lung
  m <- fit_ranger(survival::Surv(time, status) ~ age + sex, data = lung,
                  case.weights = rep(1:2, length.out = nrow(lung)), num.trees = 20)
  expect_equal(dim(predict_ranger(m, lung[1:2, ], times = 300)), c(2L, 1L))
})

test_that("ranger and orsf accept transformed terms such as log1p() in the formula", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  skip_if_not_installed("aorsf")
  skip_if_not_installed("survival")

  d <- stats::na.omit(survival::pbc[, c("time", "status", "age", "bili", "albumin")])
  d$status <- as.integer(d$status == 2)
  f <- survival::Surv(time, status) ~ age + log1p(bili) + albumin

  for (eng in c("ranger", "orsf")) {
    spec <- if (eng == "ranger") list(num.trees = 50, seed = 1) else list(n_tree = 20)
    m <- fit(f, data = d, model = eng, spec = spec)
    p <- predict(m, newdata = d[1:5, ], times = c(1000, 2000))
    expect_equal(dim(as.matrix(p)), c(5L, 2L))
    expect_true(all(as.matrix(p) >= 0 & as.matrix(p) <= 1))
  }

  set.seed(1)
  ev <- evaluate(f, d, model = "ranger", spec = list(num.trees = 50),
                 times = c(1000, 2000), metrics = "cindex",
                 resampling = cv(3, seed = 1))
  expect_true(is.finite(ev$summary$mean[ev$summary$metric == "cindex"]))
})
