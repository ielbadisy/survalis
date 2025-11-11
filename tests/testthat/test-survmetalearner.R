
test_that("fit_survmetalearner() learns nonnegative, normalized weights", {
  set.seed(123)
  n <- 80
  x <- rnorm(n)
  time   <- rexp(n, rate = exp(0.2 + 0.4 * x)) * 300
  status <- rbinom(n, 1, 0.7)
  df <- data.frame(x = x)
  times <- c(50, 150, 300)

  mk_surv <- function(lambda, tvec) exp(-outer(rep(1, n), tvec, `*`) * lambda)
  S1 <- mk_surv(lambda = exp(-0.1 + 0.3 * x) / 300, tvec = times)
  S2 <- mk_surv(lambda = exp( 0.2 - 0.2 * x) / 300, tvec = times)
  colnames(S1) <- paste0("t=", times)
  colnames(S2) <- paste0("t=", times)

  base_preds  <- list(L1 = S1, L2 = S2)
  base_models <- list(L1 = list(), L2 = list())

  meta <- fit_survmetalearner(
    base_preds = base_preds,
    time = time,
    status = status,
    times = times,
    base_models = base_models,
    formula = survival::Surv(time, status) ~ x,
    data = df
  )

  expect_s3_class(meta, "mlsurv_model")
  expect_true(inherits(meta, "survmetalearner"))
  expect_identical(meta$learners, c("L1", "L2"))
  expect_true(is.matrix(meta$weights))
  expect_identical(rownames(meta$weights), c("L1", "L2"))
  expect_setequal(colnames(meta$weights), paste0("t=", times))

  W <- meta$weights
  expect_true(all(W >= -1e-12))
  expect_true(all(abs(colSums(W) - 1) < 1e-6))
})

test_that("predict_survmetalearner() stacks base predictions correctly", {
  assign("predict_L1", function(object, newdata, times, ...) {
    lam <- exp(-0.1 + 0.3 * newdata$x) / 300
    out <- exp(-outer(lam, times, `*`))
    colnames(out) <- paste0("t=", times)
    as.data.frame(out)
  }, envir = .GlobalEnv)
  assign("predict_L2", function(object, newdata, times, ...) {
    lam <- exp(0.2 - 0.2 * newdata$x) / 300
    out <- exp(-outer(lam, times, `*`))
    colnames(out) <- paste0("t=", times)
    as.data.frame(out)
  }, envir = .GlobalEnv)
  on.exit({
    rm("predict_L1", envir = .GlobalEnv)
    rm("predict_L2", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(1)
  n <- 60
  x <- rnorm(n)
  time   <- rexp(n, rate = exp(0.2 + 0.4 * x)) * 300
  status <- rbinom(n, 1, 0.7)
  df <- data.frame(x = x)
  times <- c(40, 120, 240)

  mk_surv <- function(lambda, tvec) exp(-outer(rep(1, n), tvec, `*`) * lambda)
  S1 <- mk_surv(lambda = exp(-0.1 + 0.3 * x) / 300, tvec = times)
  S2 <- mk_surv(lambda = exp( 0.2 - 0.2 * x) / 300, tvec = times)
  colnames(S1) <- paste0("t=", times)
  colnames(S2) <- paste0("t=", times)

  base_preds  <- list(L1 = S1, L2 = S2)
  base_models <- list(L1 = list(), L2 = list())

  meta <- fit_survmetalearner(
    base_preds = base_preds,
    time = time,
    status = status,
    times = times,
    base_models = base_models,
    formula = survival::Surv(time, status) ~ x,
    data = df
  )

  newx <- data.frame(x = c(-1, 0, 1))
  pred <- predict_survmetalearner(meta, newdata = newx, times = times)

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(newx))
  expect_equal(ncol(pred), length(times))
  expect_setequal(colnames(pred), paste0("t=", times))
  M <- as.matrix(pred)
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)

  ord <- order(times)
  diffs <- t(apply(M[, ord, drop = FALSE], 1, diff))
  expect_true(all(diffs <= 1e-8))

  B1 <- as.matrix(get("predict_L1", envir = .GlobalEnv)(NULL, newx, times))
  B2 <- as.matrix(get("predict_L2", envir = .GlobalEnv)(NULL, newx, times))
  W  <- meta$weights[, paste0("t=", times), drop = FALSE] 

  stacked_check <- matrix(NA_real_, nrow = nrow(newx), ncol = length(times))
  for (j in seq_along(times)) {
    stacked_check[, j] <- cbind(B1[, j], B2[, j]) %*% W[, j]
  }
  expect_true(max(abs(stacked_check - M)) < 1e-6)
})

test_that("plot_survmetalearner_weights() returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  set.seed(2)
  n <- 20
  x <- rnorm(n)
  df <- data.frame(x = x)
  times <- c(50, 100)

  # simple base survival predictions (n x |times|)
  S1 <- exp(-outer(rep(1, n), times, `*`) * 0.002)
  S2 <- exp(-outer(rep(1, n), times, `*`) * 0.004)
  colnames(S1) <- paste0("t=", times)
  colnames(S2) <- paste0("t=", times)

  # fit a minimal meta-learner to produce `meta`
  meta <- fit_survmetalearner(
    base_preds = list(L1 = S1, L2 = S2),
    time = rexp(n, 1/200),
    status = rbinom(n, 1, 0.6),
    times = times,
    base_models = list(L1 = list(), L2 = list()),
    formula = survival::Surv(time, status) ~ x,
    data = df
  )

  # Silence any ggplot lifecycle warnings during this test only
  old <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old), add = TRUE)

  p <- suppressWarnings(plot_survmetalearner_weights(meta))
  expect_s3_class(p, "ggplot")
})


test_that("cv_survmetalearner() runs a small CV and returns a summary", {
  skip_on_cran()
  skip_if_not_installed("rsample")
  skip_if_not_installed("survival")

  assign("predict_L1", function(object, newdata, times, ...) {
    lam <- exp(-0.1 + 0.3 * newdata$x) / 300
    out <- exp(-outer(lam, times, `*`))
    colnames(out) <- paste0("t=", times)
    as.data.frame(out)
  }, envir = .GlobalEnv)
  assign("predict_L2", function(object, newdata, times, ...) {
    lam <- exp(0.2 - 0.2 * newdata$x) / 300
    out <- exp(-outer(lam, times, `*`))
    colnames(out) <- paste0("t=", times)
    as.data.frame(out)
  }, envir = .GlobalEnv)
  on.exit({
    rm("predict_L1", envir = .GlobalEnv)
    rm("predict_L2", envir = .GlobalEnv)
  }, add = TRUE)

  set.seed(7)
  n <- 60
  x <- rnorm(n)
  df <- data.frame(
    x = x,
    time = rexp(n, rate = exp(0.2 + 0.3 * x)) * 300,
    status = rbinom(n, 1, 0.7)
  )
  times <- c(100, 200)

  mk_surv <- function(lambda, tvec) exp(-outer(rep(1, n), tvec, `*`) * lambda)
  S1 <- mk_surv(lambda = exp(-0.1 + 0.3 * x) / 300, tvec = times)
  S2 <- mk_surv(lambda = exp( 0.2 - 0.2 * x) / 300, tvec = times)
  colnames(S1) <- paste0("t=", times)
  colnames(S2) <- paste0("t=", times)

  base_models <- list(L1 = list(), L2 = list())
  base_preds  <- list(L1 = S1, L2 = S2)

  res <- cv_survmetalearner(
    formula = survival::Surv(time, status) ~ x,
    data = df,
    times = times,
    base_models = base_models,
    base_preds = base_preds,
    folds = 2,
    metrics = c("cindex"),
    seed = 99,
    verbose = FALSE
  )

  expect_true(inherits(res, "cv_survmetalearner_result"))
  expect_true(is.data.frame(res$cv_results))
  expect_true(is.data.frame(res$summary))
  expect_true(all(c("mean", "sd") %in% names(res$summary)))
  expect_true(inherits(res$model, "survmetalearner"))
})



test_that("single-base meta equals base", {
  veteran <- survival::veteran
  veteran$status <- as.integer(veteran$status == 1)
  form  <- Surv(time, status) ~ age + karno + trt + celltype + prior
  times <- as.integer(quantile(veteran$time, c(0.25, 0.5, 0.75)))

  mod   <- fit_coxph(form, data = veteran)
  bp    <- predict_coxph(mod, veteran, times)

  meta  <- fit_survmetalearner(list(coxph = bp), veteran$time, veteran$status,
                               times, list(coxph = mod), form, veteran)

  pred  <- predict_survmetalearner(meta, veteran, times)
  expect_true(max(abs(as.matrix(pred) - as.matrix(bp))) < 1e-10)
  out   <- score_survmodel(meta, times, c("cindex","ibs","iae","ise"))
  expect_true(all(c("cindex","ibs","iae","ise") %in% out$metric))
  out_b <- score_survmodel(meta, times[2], "brier")
  expect_equal(out_b$metric, "brier")
})

test_that("multi-base meta shapes, weights, and CV", {
  veteran <- survival::veteran
  veteran$status <- as.integer(veteran$status == 1)
  form  <- Surv(time, status) ~ age + karno + trt + celltype + prior
  times <- as.integer(quantile(veteran$time, c(0.25, 0.5, 0.75)))

  mod   <- fit_coxph(form, data = veteran)
  bp1   <- predict_coxph(mod, veteran, times)
  bp2   <- predict_coxph(mod, veteran, times)

  meta  <- fit_survmetalearner(
            base_preds  = list(coxph = bp1, coxph_alt = bp2),
            time        = veteran$time, status = veteran$status, times = times,
            base_models = list(coxph = mod, coxph_alt = mod),
            formula     = form, data = veteran)

  expect_true(all(abs(colSums(meta$weights) - 1) < 1e-10))
  expect_true(all(meta$weights >= 0))

  pred  <- predict_survmetalearner(meta, veteran[1:5, ], times)
  expect_true(is.matrix(as.matrix(pred)))
  expect_equal(nrow(pred), 5)
  expect_equal(ncol(pred), length(times))

  expect_error(
    predict_survmetalearner(meta, veteran, c(times, max(times)+1)),
    "subset of the training times"
  )

  cvres <- cv_survmetalearner(form, veteran, times,
                              base_models = list(coxph = mod, coxph_alt = mod),
                              base_preds  = list(coxph = bp1, coxph_alt = bp2),
                              folds = 3, metrics = c("cindex","ibs"),
                              seed = 7, verbose = FALSE)
  expect_true(all(c("model","cv_results","summary") %in% names(cvres)))
  expect_true(all(c("cindex","ibs") %in% unique(cvres$cv_results$metric)))
})
