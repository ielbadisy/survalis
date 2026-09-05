f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")
}

test_that("estimate() returns a survalis_estimate with a bootstrap CI", {
  skip_verb()
  e <- estimate(f, veteran, "coxph", estimand = "rmst_diff", treatment = "trt",
                tau = 300, R = 30, refit = FALSE, seed = 1)
  expect_s3_class(e, "survalis_estimate")
  expect_equal(e$estimand, "rmst_diff")
  expect_identical(e$levels, sort(unique(veteran$trt)))
  expect_true(e$conf.low <= e$estimate && e$estimate <= e$conf.high)
  expect_equal(nrow(e$curves), length(e$curves$time))
  expect_output(print(e), "RMST difference at tau = 300")
})

test_that("estimate() supports the four estimands", {
  skip_verb()
  for (q in c("rmst_diff", "rmst_ratio", "survival_diff", "survival_ratio")) {
    e <- estimate(f, veteran, "coxph", estimand = q, treatment = "trt",
                  tau = 200, R = 15, refit = FALSE, seed = 1)
    expect_true(is.finite(e$estimate), info = q)
  }
})

test_that("estimate() validates treatment, tau and grid coverage", {
  skip_verb()
  expect_error(estimate(f, veteran, "coxph", treatment = "nope", tau = 100),
               "must name a column")
  vv <- veteran; vv$grp <- vv$age  # continuous -> not binary
  expect_error(estimate(Surv(time, status) ~ karno, vv, "coxph", treatment = "grp", tau = 100),
               "exactly two distinct values")
  expect_error(estimate(f, veteran, "coxph", treatment = "trt", tau = 100,
                        times = c(10, 20, 30)),
               "must reach `tau`")
})

test_that("estimate(refit = TRUE) refits per replicate", {
  skip_verb()
  e <- estimate(f, veteran, "coxph", treatment = "trt", tau = 200, R = 20,
                refit = TRUE, seed = 1)
  expect_true(e$refit)
  expect_length(e$boot, 20L)
  expect_s3_class(plot(e), "ggplot")
})
