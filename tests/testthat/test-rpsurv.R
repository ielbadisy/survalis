test_that("fit_rpsurv() returns an mlsurv_model and predicts a valid survmat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rpsurv")

  mod <- fit_rpsurv(Surv(time, status) ~ age + karno + celltype, veteran, df = 3)
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "rpsurv")
  expect_identical(attr(mod, "engine"), "rpsurv")

  p <- predict_rpsurv(mod, veteran[1:6, ], times = c(60, 120, 240))
  expect_equal(dim(p), c(6L, 3L))
  expect_true(all(p >= 0 & p <= 1))
  expect_identical(names(p), c("t=60", "t=120", "t=240"))
  # monotone non-increasing over time
  expect_true(all(p[, 1] >= p[, 3]))
})

test_that("rpsurv plugs into the verb interface", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rpsurv")

  expect_true("rpsurv" %in% list_survlearners()$learner)
  m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "rpsurv",
           spec = list(df = 3))
  expect_s3_class(m, "survalis_fit")
  S <- predict(m, veteran[1:4, ], times = c(90, 180))
  expect_s3_class(S, "survmat")
})

test_that("predict_rpsurv() handles newdata with unseen-but-known factor subsets", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rpsurv")

  mod <- fit_rpsurv(Surv(time, status) ~ age + celltype, veteran, df = 3)
  # subset that may not contain every celltype level
  sub <- veteran[veteran$celltype == "squamous", ][1:3, ]
  p <- predict_rpsurv(mod, sub, times = c(50, 150))
  expect_equal(dim(p), c(3L, 2L))
  expect_true(all(p >= 0 & p <= 1))
})
