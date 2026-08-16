test_that("plot_survcurve() builds an unstratified KM curve without a risk table", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("ggplot2")

  p <- plot_survcurve(
    survival::Surv(time, status) ~ 1,
    data = veteran,
    risk.table = FALSE
  )

  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})

test_that("plot_survcurve() builds a stratified KM curve with risk table and p-value", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")

  p <- plot_survcurve(
    survival::Surv(time, status) ~ trt,
    data = veteran,
    conf.int = TRUE,
    risk.table = TRUE,
    pval = TRUE
  )

  expect_s3_class(p, "patchwork")
})

test_that(".survcurve_data() produces one row at time 0 with surv = 1 per stratum", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  fit <- survival::survfit(survival::Surv(time, status) ~ trt, data = veteran)
  dt <- .survcurve_data(fit)

  expect_true(is.data.frame(dt))
  expect_true(all(c("time", "surv", "strata") %in% names(dt)))

  starts <- dt[dt$time == 0, ]
  expect_equal(nrow(starts), length(unique(dt$strata)))
  expect_true(all(starts$surv == 1))
})

test_that(".survcurve_pvalue() matches survival::survdiff() by hand", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  form <- survival::Surv(time, status) ~ trt
  pv <- .survcurve_pvalue(form, veteran)

  sd <- survival::survdiff(form, data = veteran)
  expected <- stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)

  expect_equal(pv, expected, tolerance = 1e-10)
})
