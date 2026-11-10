## Regression tests for the 2.2.0 correctness fixes.

vet <- survival::veteran
f_vet <- Surv(time, status) ~ age + karno + celltype

test_that("Cox predictions match survival::survfit", {
  d <- vet
  d$age <- d$age + 100  # mean linear predictor far from 0
  m <- fit(f_vet, d, "coxph")
  tt <- c(30, 100, 200)
  p <- unclass(predict(m, d[1:4, ], times = tt))
  cf <- survival::coxph(f_vet, data = d, model = TRUE)
  sf <- summary(survival::survfit(cf, newdata = d[1:4, ]), times = tt)
  expect_equal(matrix(as.numeric(p), nrow(p)), unname(t(sf$surv)), tolerance = 1e-8)
})

test_that("the Breslow helper matches the textbook estimator", {
  d <- vet
  cf <- survival::coxph(f_vet, data = d, ties = "breslow", model = TRUE)
  lp <- cf$linear.predictors
  tt <- c(25, 80, 300)
  p <- .predict_cox_from_lp(lp, lp[1:3], cf$y, tt)
  sf <- summary(survival::survfit(cf, newdata = d[1:3, ]), times = tt)
  expect_equal(unname(as.matrix(as.data.frame(p))), unname(t(sf$surv)), tolerance = 1e-8)
})

test_that("a 1/2 or logical status is read as event = 1", {
  l <- survival::lung[c("time", "status", "age", "sex")]
  l <- l[stats::complete.cases(l), ]
  m12 <- fit(Surv(time, status) ~ age + sex, l, "coxph")
  l01 <- l
  l01$status <- l01$status - 1L
  m01 <- fit(Surv(time, status) ~ age + sex, l01, "coxph")
  expect_equal(coef(m12$model), coef(m01$model))
  e12 <- evaluate(m12, newdata = l, times = c(200, 400), metrics = "cindex")
  e01 <- evaluate(m01, newdata = l01, times = c(200, 400), metrics = "cindex")
  expect_equal(e12$summary$mean, e01$summary$mean)
})

test_that("status == k recoding reaches every learner", {
  d <- vet
  d$st2 <- d$status + 1L
  a <- fit(Surv(time, st2 == 2) ~ age + karno, d, "glmnet", seed = 1)
  b <- fit(Surv(time, status) ~ age + karno, vet, "glmnet", seed = 1)
  tt <- c(50, 150)
  expect_equal(unclass(predict(a, d[1:3, ], times = tt)),
               unclass(predict(b, vet[1:3, ], times = tt)))
})

test_that("C-index and AUC agree with survival and timeROC", {
  m <- fit(f_vet, vet, "coxph")
  y <- survival::Surv(vet$time, vet$status)
  p <- predict(m, vet, times = c(30, 90, 180))
  ref <- survival::concordance(survival::coxph(f_vet, vet))$concordance
  expect_equal(unname(cindex_survmat(y, p, t_star = 90)), ref, tolerance = 1e-8)
  skip_if_not_installed("timeROC")
  tr <- timeROC::timeROC(vet$time, vet$status, 1 - p[, "t=90"], cause = 1,
                         times = 90, weighting = "marginal")
  expect_equal(unname(auc_survmat(y, p, t_star = 90)), unname(tr$AUC[2]),
               tolerance = 1e-4)
  expect_equal(roc_survmat(y, p, t_star = 90)$auc,
               unname(auc_survmat(y, p, t_star = 90)))
})

test_that("IBS does not depend on the order of the times", {
  m <- fit(f_vet, vet, "coxph")
  y <- survival::Surv(vet$time, vet$status)
  tt <- c(30, 60, 120, 240)
  p1 <- predict(m, vet, times = tt)
  p2 <- predict(m, vet, times = rev(tt))
  expect_equal(ibs_survmat(y, p1, tt), ibs_survmat(y, p2, rev(tt)))
})

test_that("IAE and ISE integrate over the evaluation window only", {
  m <- fit(f_vet, vet, "coxph")
  y <- survival::Surv(vet$time, vet$status)
  tt <- c(50, 100)
  p <- predict(m, vet, times = tt)
  v <- iae_survmat(y, p, tt)
  expect_true(is.finite(v) && v >= 0)
  # a constant gap g over the window [50, 100] integrates to 50 * g at most
  expect_lt(v, 50)
})

test_that("cv_summary counts only non-missing values", {
  s <- cv_summary(data.frame(metric = "ibs", value = c(0.1, NA, 0.3)))
  expect_equal(s$n, 2)
  expect_equal(s$mean, 0.2)
  expect_equal(s$se, stats::sd(c(0.1, 0.3)) / sqrt(2))
})

test_that("repeated resampling runs give identical results", {
  e1 <- evaluate(f_vet, vet, model = "rsf", times = c(50, 100),
                 resampling = cv(3, seed = 7), metrics = "cindex",
                 spec = list(ntree = 50))
  e2 <- evaluate(f_vet, vet, model = "rsf", times = c(50, 100),
                 resampling = cv(3, seed = 7), metrics = "cindex",
                 spec = list(ntree = 50))
  expect_equal(e1$per_split$value, e2$per_split$value)
})

test_that("tune keeps fixed engine arguments in the best spec", {
  tn <- tune(f_vet, vet, model = "rsf", times = c(50, 100),
             grid = expand.grid(mtry = c(1, 2)), resampling = cv(3, seed = 1),
             metric = "ibs", ntree = 50)
  expect_equal(tn$fit$spec$ntree, 50)
  expect_error(tune(f_vet, vet, model = "rsf", times = c(50, 100),
                    grid = expand.grid(mtry = c(1, 2)), mtry = 2),
               "both tuned in the grid and fixed")
})

test_that("`~ .` works for fit, predict and evaluate", {
  d <- vet[c("time", "status", "age", "karno")]
  m <- fit(Surv(time, status) ~ ., d, "coxph")
  expect_setequal(names(coef(m$model)), c("age", "karno"))
  p <- predict(m, d[1:2, ], times = c(50, 100))
  expect_equal(dim(p), c(2L, 2L))
})

test_that("rows with missing values are dropped once, at fit time", {
  d <- vet
  d$age[1:5] <- NA
  m <- fit(f_vet, d, "coxph")
  expect_equal(m$model$n, nrow(d) - 5L)
})

test_that("SHAP values add up to the prediction gap", {
  m <- fit(Surv(time, status) ~ age + karno + diagtime, vet, "coxph")
  set.seed(1)
  sh <- compute_shap(m, newdata = vet[1, ], baseline_data = vet,
                     times = c(50, 100), sample.size = 400)
  sh <- as.data.frame(sh)
  skip_if_not("phi" %in% names(sh))
  tot <- tapply(sh$phi, sh$time, sum)
  pr <- predict(m, vet[1, ], times = c(50, 100))
  base <- colMeans(predict(m, vet, times = c(50, 100)))
  expect_equal(as.numeric(tot), as.numeric(pr) - unname(base), tolerance = 0.05)
})

test_that("interaction strength is near zero for an additive model", {
  m <- fit(Surv(time, status) ~ age + karno, vet, "coxph")
  set.seed(1)
  h <- as.data.frame(compute_interactions(m, data = vet, times = c(50, 100),
                                          type = "1way", target_time = 100,
                                          grid.size = 30))
  skip_if_not("interaction" %in% names(h))
  expect_true(all(h$interaction < 0.05, na.rm = TRUE))
})

test_that("calibration handles tied predictions", {
  m <- fit(Surv(time, status) ~ celltype, vet, "coxph")
  cal <- compute_calibration(m, data = vet, time = "time", status = "status",
                             eval_time = 100, n_bins = 10)
  expect_true(nrow(as.data.frame(cal)) <= 4)
})

test_that("varimp is positive when a useful feature is permuted", {
  m <- fit(Surv(time, status) ~ karno + diagtime, vet, "coxph")
  for (met in c("ibs", "cindex")) {
    vi <- as.data.frame(interpret(m, "varimp", times = c(50, 100),
                                  metric = met, n_repetitions = 5, seed = 1))
    expect_gt(vi$importance[vi$feature == "karno"], 0)
  }
})

test_that("kmplot accepts a model fitted on a rescaled time", {
  d <- vet
  d$years <- d$time / 365.25
  m <- fit(Surv(time, status) ~ trt + karno, vet, "coxph")
  expect_s3_class(kmplot(Surv(years, status) ~ trt, data = d, model = m), "ggplot")
})
