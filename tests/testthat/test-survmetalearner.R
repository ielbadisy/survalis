
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

test_that("fit_survmetalearner() ignores subjects censored before t*", {
  set.seed(1)
  n <- 200
  time <- rexp(n) * 100
  status <- rbinom(n, 1, 0.5)
  df <- data.frame(x = rnorm(n))
  times <- c(30, 80)
  known <- sapply(times, function(t) time > t | status == 1)
  truth <- sapply(times, function(t) as.numeric(time > t))
  # a: exact where status at t* is known, predicts "alive" for early censored
  S_a <- truth; S_a[!known] <- 1
  # b: matches the naive target I(T > t*) everywhere, noisy on a few known subjects
  S_b <- truth
  noisy <- which(known)[1:5]
  S_b[noisy] <- 1 - S_b[noisy]
  m <- fit_survmetalearner(base_preds = list(a = S_a, b = S_b), time = time,
                           status = status, times = times,
                           base_models = list(a = list(), b = list()),
                           formula = survival::Surv(time, status) ~ x, data = df)
  expect_true(all(is.finite(m$weights)))
  expect_equal(unname(colSums(m$weights)), rep(1, length(times)))
  expect_gt(min(m$weights["a", ]), 0.99)
})
