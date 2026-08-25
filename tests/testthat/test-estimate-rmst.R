test_that("estimate_rmst: marginal mode (trt_col = NULL) returns a standardized RMST with a CI", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
  times <- seq(10, 300, by = 10)

  est <- estimate_rmst(mod, data = veteran, tau = 200, times = times, R = 15, seed = 1)

  expect_s3_class(est, "rmst_estimate")
  expect_identical(est$mode, "marginal")
  expect_null(est$trt_col)
  expect_null(est$levels)
  expect_true(is.finite(est$estimate))
  expect_gt(est$estimate, 0)
  expect_length(est$boot, 15)
  expect_lte(est$lo, est$hi)
})

test_that("estimate_rmst: contrast mode (trt_col set) g-computes a two-arm contrast", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
  times <- seq(10, 300, by = 10)

  rc <- estimate_rmst(mod, data = veteran, tau = 200, times = times,
                      trt_col = "trt", R = 15, seed = 1)

  expect_s3_class(rc, "rmst_estimate")
  expect_identical(rc$mode, "contrast")
  expect_true(is.finite(rc$estimate))
  expect_length(rc$boot, 15)
  expect_lte(rc$lo, rc$hi)
  expect_equal(rc$levels, c(1, 2))
})

test_that("estimate_rmst validates inputs", {
  mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
  times <- seq(10, 300, by = 10)

  expect_error(
    estimate_rmst(mod, data = veteran, tau = 200, times = times, trt_col = "not_a_col"),
    "trt_col"
  )
  expect_error(
    estimate_rmst(mod, data = veteran, tau = -1, times = times),
    "tau"
  )
})

test_that("plot_rmst dispatches marginal vs. contrast by trt_col, like plot_survmat's group arg", {
  testthat::skip_on_cran()

  mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
  times <- seq(10, 300, by = 10)

  p_marg <- plot_rmst(mod, data = veteran, tau_seq = c(100, 200), times = times, R = 5, seed = 1)
  expect_s3_class(p_marg, "ggplot")
  expect_equal(nrow(attr(p_marg, "curve")), 2)

  p_contrast <- plot_rmst(mod, data = veteran, tau_seq = c(100, 200), times = times,
                                trt_col = "trt", R = 5, seed = 1)
  expect_s3_class(p_contrast, "ggplot")
  expect_equal(nrow(attr(p_contrast, "curve")), 2)
})
