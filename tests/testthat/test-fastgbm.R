test_that("fit_fastgbm() returns an mlsurv_model and predicts a valid survmat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fastgbm")

  mod <- fit_fastgbm(Surv(time, status) ~ age + karno + celltype, veteran,
                     ntrees = 40, pexp_bins = 8)
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "fastgbm")
  expect_identical(attr(mod, "engine"), "fastgbm")

  p <- predict_fastgbm(mod, veteran[1:6, ], times = c(60, 120, 240))
  expect_equal(dim(p), c(6L, 3L))
  expect_true(all(p >= 0 & p <= 1))
  expect_identical(names(p), c("t=60", "t=120", "t=240"))
  # monotone non-increasing over time
  expect_true(all(p[, 1] >= p[, 3]))
})

test_that("fastgbm plugs into the verb interface", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fastgbm")

  expect_true("fastgbm" %in% list_survlearners()$learner)
  m <- fit(Surv(time, status) ~ age + karno, veteran, "fastgbm",
           spec = list(ntrees = 40, pexp_bins = 8))
  expect_s3_class(m, "survalis_fit")
  S <- predict(m, veteran[1:4, ], times = c(90, 180))
  expect_s3_class(S, "survmat")
})

test_that("predict_fastgbm() rejects a non-pexp objective", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fastgbm")
  mod <- fit_fastgbm(Surv(time, status) ~ age, veteran, objective = "cox", ntrees = 20)
  expect_error(predict_fastgbm(mod, veteran[1:2, ], times = 90), "objective = 'pexp'")
})
