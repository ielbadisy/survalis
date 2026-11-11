sim_svm <- function(n = 150, censor = TRUE, seed = 1) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)))
  lt <- 1 + 0.8 * d$x1 - 0.5 * d$x2^2 + (d$g == "b") + rnorm(n, 0, 0.5)
  ct <- if (censor) log(rexp(n, 0.1)) else rep(Inf, n)
  d$time <- exp(pmin(lt, ct))
  d$status <- as.numeric(lt <= ct)
  d
}

svm_parts <- function(mod, data) {
  enc <- mod$model$encoder
  mf <- stats::model.frame(mod$formula, data)
  Z <- survalis:::.ssvm_encode(enc, mf)
  K <- survalis:::.ssvm_kernel(enc, Z, Z)
  list(K = K, y = log(data$time), f = survalis:::.ssvm_predict_f(mod$model, data))
}

test_that("fit_survsvm() returns a proper mlsurv_model wrapper", {
  mod <- fit_survsvm(Surv(time, status) ~ age + celltype + karno,
                     data = veteran, gamma.mu = 0.1, kernel = "lin_kernel")
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "survsvm")
  expect_identical(attr(mod, "engine"), "survalis")
  expect_true(length(mod$model$beta) > 0)
})

test_that("without censoring the solution is the LS-SVM closed form", {
  d <- sim_svm(censor = FALSE)
  n <- nrow(d)
  for (k in c("lin_kernel", "add_kernel", "rbf_kernel", "poly_kernel")) {
    mod <- fit_survsvm(Surv(time, status) ~ x1 + x2 + g, d, kernel = k, gamma.mu = 2)
    p <- svm_parts(mod, d)
    M <- rbind(c(0, rep(1, n)), cbind(1, p$K + diag(n) / 2))
    s <- solve(M, c(0, p$y))
    expect_equal(p$f, drop(p$K %*% s[-1]) + s[1], tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("with censoring the solution satisfies the optimality conditions", {
  d <- sim_svm()
  for (g in c(0.01, 1, 100)) {
    mod <- fit_survsvm(Surv(time, status) ~ x1 + x2 + g, d, kernel = "add_kernel", gamma.mu = g)
    p <- svm_parts(mod, d)
    r <- p$y - p$f
    active <- d$status == 1 | r > 0
    # beta_i = gamma * r_i on the active set, 0 elsewhere; sum(beta) = 0
    beta <- numeric(nrow(d))
    beta[active] <- g * r[active]
    expect_equal(drop(p$K %*% beta) + mod$model$b, p$f, tolerance = 1e-7, ignore_attr = TRUE)
    expect_equal(sum(beta), 0, tolerance = 1e-7)
    expect_equal(sum(mod$model$beta), 0, tolerance = 1e-7)
  }
})

test_that("linear kernel with a large error weight recovers least squares", {
  d <- sim_svm(censor = FALSE)
  mod <- fit_survsvm(Surv(time, status) ~ x1 + x2 + g, d, kernel = "lin_kernel", gamma.mu = 1e6)
  ols <- stats::fitted(stats::lm(log(time) ~ x1 + x2 + g, d))
  expect_equal(survalis:::.ssvm_predict_f(mod$model, d), ols, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("kernels are positive semi-definite", {
  d <- sim_svm()
  for (k in c("lin_kernel", "add_kernel", "rbf_kernel", "poly_kernel")) {
    mod <- fit_survsvm(Surv(time, status) ~ x1 + x2 + g, d, kernel = k)
    ev <- eigen(svm_parts(mod, d)$K, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(min(ev), -1e-8)
  }
})

test_that("add_kernel tolerates constant predictors and ordered factors", {
  d <- sim_svm()
  d$const <- 3
  d$o <- factor(sample(c("lo", "mid", "hi"), nrow(d), TRUE),
                levels = c("lo", "mid", "hi"), ordered = TRUE)
  mod <- fit_survsvm(Surv(time, status) ~ x1 + const + o + g, d)
  pr <- predict_survsvm(mod, d[1:3, ], times = c(1, 3))
  expect_true(all(is.finite(as.matrix(pr))))
})

test_that("predict_survsvm() returns calibrated, well-formed probabilities", {
  d <- sim_svm(n = 300)
  mod <- fit_survsvm(Surv(time, status) ~ x1 + x2 + g, d)
  nd <- sim_svm(n = 3000, censor = FALSE, seed = 2)
  times <- stats::quantile(nd$time, c(0.25, 0.5, 0.75))
  for (ds in c("km", "lognormal", "exp", "weibull")) {
    M <- as.matrix(predict_survsvm(mod, nd, times = times, dist = ds, shape = 1.3))
    expect_equal(dim(M), c(nrow(nd), 3L))
    expect_true(all(is.finite(M) & M >= 0 & M <= 1))
    expect_true(all(apply(M, 1, diff) <= 1e-12))
  }
  # the average predicted curve matches the Kaplan-Meier curve of the sample
  M <- as.matrix(predict_survsvm(mod, nd, times = times))
  km <- summary(survival::survfit(Surv(time, status) ~ 1, d), times = times)$surv
  expect_lt(max(abs(colMeans(M) - km)), 0.05)
})

test_that("fitting does not consume the random number generator", {
  d <- sim_svm()
  set.seed(9); a <- runif(1)
  set.seed(9); invisible(fit_survsvm(Surv(time, status) ~ x1 + g, d)); b <- runif(1)
  expect_identical(a, b)
})

test_that("bad inputs give clear errors", {
  d <- sim_svm()
  expect_error(fit_survsvm(Surv(time, status) ~ x1, d, type = "vanbelle1"), "regression")
  expect_error(fit_survsvm(Surv(time, status) ~ x1, d, gamma.mu = -1), "positive")
  d0 <- d; d0$time[1] <- 0
  expect_error(fit_survsvm(Surv(time, status) ~ x1, d0), "positive")
  mod <- fit_survsvm(Surv(time, status) ~ x1 + g, d)
  expect_error(predict_survsvm(mod, d[1:3, "x1", drop = FALSE], times = 1))
  nd <- d[1:2, ]; nd$x1[1] <- NA
  expect_error(predict_survsvm(mod, nd, times = 1), "missing")
})

test_that("tune_survsvm() returns ranked results and refits the best model", {
  res <- tune_survsvm(Surv(time, status) ~ age + celltype + karno, data = veteran,
                      times = c(100, 300), param_grid = list(gamma.mu = c(0.01, 0.1),
                      kernel = "add_kernel"), folds = 3)
  expect_s3_class(res, "tune_surv")
  expect_equal(nrow(res), 2L)
  best <- tune_survsvm(Surv(time, status) ~ age + celltype + karno, data = veteran,
                       times = c(100, 300), param_grid = list(gamma.mu = 0.1,
                       kernel = "lin_kernel"), folds = 3, refit_best = TRUE)
  expect_s3_class(best, "mlsurv_model")
  expect_false(is.null(attr(best, "tuning_results")))
})

test_that("fit_survsvm() stores the outcome column names", {
  s <- fit_survsvm(survival::Surv(time, status) ~ age + karno, survival::veteran)
  expect_identical(c(s$time, s$status), c("time", "status"))
})
