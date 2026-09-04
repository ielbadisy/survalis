f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")
}

test_that("evaluate() on a formula returns a survalis_eval with a summary", {
  skip_verb()
  e <- evaluate(f, veteran, "coxph", times = c(90, 180, 365), resampling = cv(3, seed = 1))
  expect_s3_class(e, "survalis_eval")
  expect_setequal(e$summary$metric, c("cindex", "ibs"))
  expect_equal(e$times, c(90, 180, 365))
  expect_equal(nrow(e$per_split), 3L * 2L)
  expect_true(all(e$summary$mean >= 0))
  expect_output(print(e), "survalis_eval.*coxph")
})

test_that("evaluate() on a fit reuses the stored recipe", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  e1 <- evaluate(m, times = c(90, 180), resampling = cv(3, seed = 42))
  e2 <- evaluate(f, veteran, "coxph", times = c(90, 180), resampling = cv(3, seed = 42))
  expect_equal(e1$per_split$value, e2$per_split$value)
})

test_that("evaluate() honors the resampling spec (holdout, repeats)", {
  skip_verb()
  eh <- evaluate(f, veteran, "coxph", times = 180, resampling = holdout(0.7, seed = 1),
                 metrics = "cindex")
  expect_equal(nrow(eh$per_split), 1L)

  er <- evaluate(f, veteran, "coxph", times = c(90, 180), resampling = cv(3, repeats = 2, seed = 1))
  expect_equal(length(unique(er$per_split$id)), 6L)
})

test_that("evaluate() default grid is recorded", {
  skip_verb()
  e <- evaluate(f, veteran, "coxph", resampling = cv(3, seed = 1))
  expect_equal(e$times, default_times(veteran$time, veteran$status))
})
