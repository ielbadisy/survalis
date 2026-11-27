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
