test_that("fit_coxboost(cv = TRUE) uses the cross-validated number of steps", {
  skip_on_cran()
  skip_if_not_installed("CoxBoost")
  d <- survival::veteran; f <- Surv(time, status) ~ age + karno + diagtime
  expect_equal(fit_coxboost(f, d)$model$stepno, 100)
  set.seed(1); m <- fit_coxboost(f, d, cv = TRUE, K = 5, stepno = 60)
  x <- as.matrix(d[, c("age", "karno", "diagtime")])
  set.seed(1)
  opt <- CoxBoost::cv.CoxBoost(d$time, d$status, x, maxstepno = 60, K = 5,
                               penalty = 9 * sum(d$status))$optimal.step
  expect_equal(m$model$stepno, opt)
  expect_equal(dim(as.matrix(predict(m$model, newdata = x[1:3, ], type = "lp"))), c(1L, 3L))
})

test_that("coxboost takes unpen.index by name and fixed folds", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("CoxBoost")
  d <- veteran
  x <- model.matrix(~ age + karno + diagtime + celltype, d)[, -1]
  set.seed(1); fid <- sample(rep(1:3, length.out = nrow(d)))
  fo <- split(seq_len(nrow(d)), fid)
  ref_cv <- CoxBoost::cv.CoxBoost(d$time, d$status, x, maxstepno = 40, K = 3, folds = fo,
                                  penalty = 9 * sum(d$status), unpen.index = 2)
  ref <- CoxBoost::CoxBoost(d$time, d$status, x, stepno = ref_cv$optimal.step,
                            penalty = 9 * sum(d$status), unpen.index = 2)
  f <- Surv(time, status) ~ age + karno + diagtime + celltype
  m <- fit(f, d, "coxboost", spec = list(stepno = 40, folds = fid, unpen.index = "karno"))
  expect_identical(m$cv$optimal.step, ref_cv$optimal.step)
  expect_equal(m$cv$mean.logplik, ref_cv$mean.logplik)
  expect_equal(coef(m$model), coef(ref))
  # a list of folds and positions give the same fit; a factor name selects all its dummies
  m2 <- fit(f, d, "coxboost", spec = list(stepno = 40, folds = fo, unpen.index = 2))
  expect_equal(coef(m2$model), coef(ref))
  m3 <- fit(f, d, "coxboost", spec = list(stepno = 10, unpen.index = "celltype"))
  expect_identical(m3$model$unpen.index, 4:6)
  expect_error(fit(f, d, "coxboost", spec = list(unpen.index = "nope")), "unpen.index")
  expect_error(fit(f, d, "coxboost", spec = list(folds = 1:3)), "one fold label")
  expect_null(fit(f, d, "coxboost", spec = list(stepno = 5))$cv)
})
