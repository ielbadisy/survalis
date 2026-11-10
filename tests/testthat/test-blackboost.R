test_that("tune_blackboost() returns a minimal grid with metrics and refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("mboost")
  skip_if_not_installed("partykit")
  skip_if_not_installed("survival")

  df <- survival::veteran
  set.seed(777)

  Surv <- survival::Surv

  grid <- expand.grid(
    mstop = c(25, 50),
    nu = 0.1,
    maxdepth = 2,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  res <- tune_blackboost(
    formula = Surv(time, status) ~ age + karno + celltype,   
    data = df,
    times = c(30, 90, 150),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 777,
    refit_best = FALSE
  )

  expect_s3_class(res, "tuned_surv")
  expect_true(all(c("mstop", "nu", "maxdepth") %in% names(res)))
  expect_true(all(c("cindex", "ibs") %in% names(res)))

  best_mod <- tune_blackboost(
    formula = Surv(time, status) ~ age + karno + celltype,  
    data = df,
    times = c(30, 90, 150),
    param_grid = grid,
    metrics = c("cindex", "ibs"),
    folds = 2,
    seed = 777,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_identical(best_mod$learner, "blackboost")
  expect_true(inherits(best_mod$model, "mboost"))
})

test_that("predict_blackboost() handles one row, no outcome columns and unsorted times", {
  skip_on_cran()
  skip_if_not_installed("mboost")
  skip_if_not_installed("partykit")

  df <- survival::veteran
  mod <- fit_blackboost(Surv(time, status) ~ age + karno + trt, data = df)
  tt <- c(200, 10, 100, 5000)
  b <- as.matrix(predict_blackboost(mod, df, times = tt))
  expect_identical(colnames(b), paste0("t=", tt))
  one <- as.matrix(predict_blackboost(mod, df[1, ], times = tt))
  expect_equal(dim(one), c(1L, 4L))
  expect_equal(unname(one[1, ]), unname(b[1, ]))
  nd <- df[1:5, c("age", "karno", "trt")]
  expect_equal(unname(as.matrix(predict_blackboost(mod, nd, times = tt))), unname(b[1:5, ]))
  # beyond the last event time survival is carried forward per observation
  far <- as.matrix(predict_blackboost(mod, df[1:3, ], times = c(999, 5000)))
  expect_equal(far[, 1], far[, 2])
})
