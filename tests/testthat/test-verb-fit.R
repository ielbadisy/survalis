f <- Surv(time, status) ~ age + karno + celltype

skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")
}

test_that("fit() validates the outcome and the model id", {
  skip_verb()
  expect_error(fit(time ~ age, veteran, "coxph"), "Surv\\(\\) outcome")
  expect_error(fit(f, veteran, "not_a_learner"), "Unknown model")
})

test_that("fit() returns a survalis_fit that keeps the mlsurv_model contract", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  expect_s3_class(m, "survalis_fit")
  expect_s3_class(m, "mlsurv_model")
  expect_identical(m$model_id, "coxph")
  expect_identical(m$time, "time")
  expect_identical(m$status, "status")
  expect_true(is.list(m$spec))
  expect_output(print(m), "survalis_fit.*coxph")
})

test_that("fit() merges ... into spec once, with a warning", {
  skip_verb()
  expect_warning(fit(f, veteran, "ranger", num.trees = 50), "pass engine arguments")
})

test_that("fit(seed=) makes stochastic learners reproducible", {
  skip_verb()
  m1 <- fit(f, veteran, "ranger", spec = list(num.trees = 100), seed = 1)
  m2 <- fit(f, veteran, "ranger", spec = list(num.trees = 100), seed = 1)
  p1 <- predict(m1, veteran[1:10, ], times = c(90, 180))
  p2 <- predict(m2, veteran[1:10, ], times = c(90, 180))
  expect_equal(unclass(as.matrix(p1)), unclass(as.matrix(p2)), ignore_attr = TRUE)
})
