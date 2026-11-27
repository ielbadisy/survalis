test_that("predict_aalen() validates inputs without mocking", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("timereg")
  testthat::skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_aalen(survival::Surv(time, status) ~ trt + karno + age, data = df, max.time = 600)

  expect_error(predict_aalen(list(model = mod$model), newdata = df, times = 1:5))

  bad <- mod
  bad$learner <- "coxph"

  expect_error(predict_aalen(bad, newdata = df, times = 1:5))

  bad_new <- subset(df, select = -age)
  
  expect_error(predict_aalen(mod, newdata = bad_new, times = c(50, 100)))

  expect_error(predict_aalen(mod, newdata = df[1:5, ], times = c("a", "b")))
})

test_that("aalen S(t) does not depend on the other requested times", {
  r <- with(survival::rotterdam, data.frame(time = dtime / 365.25, status = death, age, size, grade,
                                            lnodes = log1p(nodes), lpgr = log1p(pgr), ler = log1p(er)))
  m <- fit(Surv(time, status) ~ age + size + grade + lnodes + lpgr + ler, r, model = "aalen")
  a <- predict(m, r[1:20, ], times = c(2, 5, 8))
  b <- predict(m, r[1:20, ], times = 8)
  expect_equal(unclass(a)[, 3], unclass(b)[, 1], tolerance = 1e-10, ignore_attr = TRUE)
  # S(t) = min over model jump times s <= t of exp(-x'B(s))
  nd <- .aalen_scale(.aalen_relevel(r[1:20, ], m$xlevels), m$scaling)
  jt <- m$model$cum[, 1]
  raw <- timereg::predict.aalen(m$model, newdata = nd, times = c(jt[jt <= 8], 8), se = FALSE)$S0
  expect_equal(unclass(b)[, 1], apply(raw, 1, min), tolerance = 1e-10, ignore_attr = TRUE)
})
