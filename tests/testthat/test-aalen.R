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
