test_that("fit_nbsurv() returns an mlsurv_model and predicts a valid survmat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("nbsurv")

  mod <- fit_nbsurv(Surv(time, status) ~ age + karno + celltype, veteran)
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "nbsurv")
  expect_identical(attr(mod, "engine"), "nbsurv")

  p <- predict_nbsurv(mod, veteran[1:6, ], times = c(60, 120, 240))
  expect_equal(dim(p), c(6L, 3L))
  expect_true(all(p >= 0 & p <= 1))
  expect_identical(names(p), c("t=60", "t=120", "t=240"))
  # monotone non-increasing over time
  expect_true(all(p[, 1] >= p[, 3]))
})

test_that("nbsurv plugs into the verb interface", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("nbsurv")

  expect_true("nbsurv" %in% list_survlearners()$learner)
  m <- fit(Surv(time, status) ~ age + karno, veteran, "nbsurv")
  expect_s3_class(m, "survalis_fit")
  S <- predict(m, veteran[1:4, ], times = c(90, 180))
  expect_s3_class(S, "survmat")
})
