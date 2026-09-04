test_that(".is_surv_lhs detects a Surv() outcome syntactically", {
  expect_true(.is_surv_lhs(Surv(time, status) ~ trt))
  expect_true(.is_surv_lhs(survival::Surv(time, status) ~ trt))
  expect_false(.is_surv_lhs(mpg ~ wt))
  expect_false(.is_surv_lhs(~wt))
})

test_that("plot() on a Surv() formula draws a Kaplan-Meier curve", {
  testthat::skip_on_cran()
  p <- plot(Surv(time, status) ~ trt, data = veteran)
  expect_s3_class(p, "ggplot")
})

test_that("plot() on a Surv() formula requires data", {
  expect_error(plot(Surv(time, status) ~ trt), "requires a data frame")
})

test_that("plot() on an ordinary formula still uses base graphics, unaffected", {
  testthat::skip_on_cran()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot(mpg ~ wt, data = mtcars)
  expect_null(out)
})
