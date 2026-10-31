
test_that("fit_bnnsurv() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_bnnsurv(
    survival::Surv(time, status) ~ age + karno + diagtime + prior,
    data = df,
    k = 3, num_base_learners = 5
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "bnnsurv")
  expect_identical(mod$engine, "survalis")
  expect_true(!is.null(mod$model))
})

test_that("predict_bnnsurv() (native grid) returns bounded, monotone survival", {
  skip_on_cran()
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
  expect_identical(best_mod$engine, "survalis")
  expect_true(!is.null(best_mod$model))
})

test_that("one learner on all rows with k = n is the Kaplan-Meier estimator", {
  d <- survival::veteran
  mod <- fit_bnnsurv(
    survival::Surv(time, status) ~ age + karno,
    data = d, k = nrow(d), num_base_learners = 1,
    replace = FALSE, sample_fraction = 1
  )
  tt <- c(0.5, 1, 7, 50, 100.3, 500, 999, 2000)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d)
  expect_equal(
    unname(unlist(predict_bnnsurv(mod, d[1, ], times = tt))),
    summary(km, times = tt, extend = TRUE)$surv,
    tolerance = 1e-12
  )
})

test_that(".weighted_km() matches survfit() with case weights", {
  d <- survival::veteran
  set.seed(3)
  w <- stats::runif(nrow(d))
  grid <- sort(unique(d$time))
  S <- survalis:::.weighted_km(d$time, d$status, w,
                               row = rep(1L, nrow(d)), nrow = 1L, grid = grid)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d, weights = w)
  expect_equal(S[, 1], summary(km, times = grid)$surv, tolerance = 1e-12)
})

test_that("tied neighbours share weight, so predictions ignore row order", {
  d <- survival::veteran
  f <- survival::Surv(time, status) ~ karno
  m1 <- fit_bnnsurv(f, d, k = 3, num_base_learners = 1,
                    replace = FALSE, sample_fraction = 1)
  perm <- rev(seq_len(nrow(d)))
  m2 <- fit_bnnsurv(f, d[perm, ], k = 3, num_base_learners = 1,
                    replace = FALSE, sample_fraction = 1)
  expect_equal(predict_bnnsurv(m1, d[1:20, ], times = c(30, 100, 300)),
               predict_bnnsurv(m2, d[1:20, ], times = c(30, 100, 300)),
               tolerance = 1e-12)
})

test_that("predictions read the step function at the requested times", {
  d <- survival::veteran
  mod <- fit_bnnsurv(survival::Surv(time, status) ~ age + karno, data = d,
                     k = 10, num_base_learners = 3)
  full <- as.matrix(predict_bnnsurv(mod, d[1:4, ]))
  grid <- mod$model$grid
  at <- as.matrix(predict_bnnsurv(mod, d[1:4, ], times = c(grid[5] + 0.5, grid[5])))
  expect_equal(unname(at[, 1]), unname(full[, 5]))
  expect_equal(unname(at[, 2]), unname(full[, 5]))
  before <- as.matrix(predict_bnnsurv(mod, d[1:4, ], times = min(grid) / 2))
  expect_true(all(before == 1))
})

test_that("factors are one-hot encoded and new data reuse the training levels", {
  d <- survival::veteran
  d$celltype <- factor(d$celltype)
  f <- survival::Surv(time, status) ~ celltype + karno + age
  set.seed(11)
  mod <- fit_bnnsurv(f, d, k = 5, num_base_learners = 6,
                     num_features_per_base_learner = 1)
  cols <- mod$model$columns
  expect_true(all(c("celltypesmallcell", "celltypeadeno", "celltypelarge") %in% cols))
  # a factor is drawn as one predictor: all of its dummies or none
  for (bl in mod$model$learners) {
    dummies <- grepl("^celltype", cols[bl$cols])
    expect_true(all(dummies) || !any(dummies))
  }
  nd <- d[c(5, 90), ]
  nd$celltype <- factor(as.character(nd$celltype))
  p_sub <- predict_bnnsurv(mod, nd, times = c(30, 100))
  p_full <- predict_bnnsurv(mod, d[c(5, 90), ], times = c(30, 100))
  expect_equal(p_sub, p_full, ignore_attr = TRUE)

  bad <- d[1, ]
  bad$celltype <- factor("unknown")
  expect_error(predict_bnnsurv(mod, bad, times = 30))
})

test_that("kernel weights and the euclidean metric give valid curves", {
  d <- survival::veteran
  mod <- fit_bnnsurv(survival::Surv(time, status) ~ age + karno + prior, data = d,
                     metric = "euclidean", weighting_function = function(x) exp(-x))
  M <- as.matrix(predict_bnnsurv(mod, d[1:10, ], times = c(30, 100, 300)))
  expect_true(all(M >= 0 & M <= 1))
  expect_true(all(t(apply(M, 1, diff)) <= 0))
  expect_error(
    fit_bnnsurv(survival::Surv(time, status) ~ age, data = d,
                weighting_function = function(x) -x) |>
      predict_bnnsurv(d[1:2, ], times = 30),
    "non-negative"
  )
})
