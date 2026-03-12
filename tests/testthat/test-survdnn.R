
test_that("fit_survdnn() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("survdnn")
  skip_if_not_installed("torch")
  skip_if_not_installed("survival")
  skip_if_not(torch::torch_is_installed(), "Torch backend is not installed")

  df <- survival::veteran
  Surv <- survival::Surv

  mod <- fit_survdnn(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    loss = "cox",
    hidden = c(16L, 8L),
    epochs = 5,             
    lr = 1e-3,
    verbose = FALSE
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "survdnn")
  expect_identical(attr(mod, "engine"), "survdnn")
  expect_true(!is.null(mod$model))
})

test_that("predict_survdnn() returns bounded, monotone-ish survival at requested times", {
  skip_on_cran()
  skip_if_not_installed("survdnn")
  skip_if_not_installed("torch")
  skip_if_not_installed("survival")
  skip_if_not(torch::torch_is_installed(), "Torch backend is not installed")

  df  <- survival::veteran
  Surv <- survival::Surv
  mod <- fit_survdnn(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    loss = "cox",
    hidden = c(16L, 8L),
    epochs = 5,
    lr = 1e-3,
    verbose = FALSE
  )

  times <- c(150, 30, 300, 200)
  newx  <- df[1:6, ]

  pred <- predict_survdnn(mod, newdata = newx, times = times)

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(newx))
  expect_equal(ncol(pred), length(times))

  norm_names <- sub("^t=", "", colnames(pred))
  expect_setequal(norm_names, as.character(times))

  M <- as.matrix(pred)
  mode(M) <- "numeric"
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)

  ord <- order(as.numeric(norm_names))
  M_ord <- M[, ord, drop = FALSE]
  diffs <- t(apply(M_ord, 1, diff))
  expect_true(all(diffs <= 1e-6))
})

test_that("predict_survdnn() exposes lp and risk prediction types", {
  skip_on_cran()
  skip_if_not_installed("survdnn")
  skip_if_not_installed("torch")
  skip_if_not_installed("survival")
  skip_if_not(torch::torch_is_installed(), "Torch backend is not installed")

  df  <- survival::veteran
  Surv <- survival::Surv
  mod <- fit_survdnn(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    loss = "cox",
    hidden = c(16L, 8L),
    epochs = 5,
    lr = 1e-3,
    optimizer = "adam",
    dropout = 0.1,
    batch_norm = TRUE,
    verbose = FALSE
  )

  lp <- predict_survdnn(mod, newdata = df[1:4, ], type = "lp")
  risk <- predict_survdnn(mod, newdata = df[1:4, ], times = 90, type = "risk")

  expect_type(lp, "double")
  expect_length(lp, 4)
  expect_type(risk, "double")
  expect_length(risk, 4)
  expect_true(all(is.finite(lp)))
  expect_true(all(is.finite(risk)))
  expect_gte(min(risk), -1e-12)
  expect_lte(max(risk), 1 + 1e-12)
})

test_that("predict_survdnn() errors when newdata lacks required predictors", {
  skip_on_cran()
  skip_if_not_installed("survdnn")
  skip_if_not_installed("torch")
  skip_if_not_installed("survival")
  skip_if_not(torch::torch_is_installed(), "Torch backend is not installed")

  df  <- survival::veteran
  Surv <- survival::Surv
  mod <- fit_survdnn(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    loss = "cox",
    hidden = c(8L),
    epochs = 3,
    lr = 1e-3,
    verbose = FALSE
  )

  bad_new <- df[1:3, c("age", "karno")]  
  expect_error(
    predict_survdnn(mod, newdata = bad_new, times = c(50, 100)),
    regexp = paste0(
      "object .* not found|undefined columns|subset.*select|",
      "do not match|newdata.*missing|terms|model\\.matrix|",
      'Value of covariate ".*" not supplied in "newdata"|',
      "names of independing variables.*do not match|",
      "independent variables not found"
    )
  )
})

test_that("tune_survdnn() returns a compact grid and refit_best yields a fitted model", {
  skip_on_cran()
  skip_if_not_installed("survdnn")
  skip_if_not_installed("torch")
  skip_if_not_installed("survival")
  skip_if_not(torch::torch_is_installed(), "Torch backend is not installed")

  Surv <- survival::Surv
  df   <- survival::veteran
  set.seed(2025)

  grid <- list(
    hidden     = list(c(8L), c(16L, 8L)), 
    lr         = c(1e-3),
    activation = c("relu"),
    epochs     = c(3),                      
    loss       = c("cox")
  )

  res <- tune_survdnn(
    formula = Surv(time, status) ~ age + karno + celltype,
    data = df,
    times = c(90),
    metrics = c("cindex", "ibs"),
    param_grid = grid,
    folds = 2,
    seed = 2025,
    refit_best = FALSE
  )

  expect_s3_class(res, "tuned_surv")
  expect_true(all(c("hidden", "lr", "activation", "epochs", "loss") %in% names(res)))
  expect_true(all(c("cindex", "ibs") %in% names(res)))
  expect_true(nrow(res) >= 1)

  best_mod <- tune_survdnn(
    formula = Surv(time, status) ~ age + karno + celltype,
    data = df,
    times = c(90),
    metrics = c("cindex", "ibs"),
    param_grid = grid,
    folds = 2,
    seed = 2025,
    refit_best = TRUE
  )

  expect_s3_class(best_mod, "mlsurv_model")
  expect_identical(best_mod$learner, "survdnn")
  expect_identical(attr(best_mod, "engine"), "survdnn")
  expect_true(!is.null(attr(best_mod, "tuning_results")))
})
