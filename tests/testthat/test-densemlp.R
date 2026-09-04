test_that("fit_densemlp() returns an mlsurv_model and predicts a valid survmat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("densemlp")

  mod <- fit_densemlp(Surv(time, status) ~ age + karno + celltype, veteran,
                      n_bins = 8, epochs = 15, hidden_units = c(16, 8))
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "densemlp")
  expect_identical(attr(mod, "engine"), "densemlp")

  p <- predict_densemlp(mod, veteran[1:6, ], times = c(60, 120, 240))
  expect_equal(dim(p), c(6L, 3L))
  expect_true(all(p >= 0 & p <= 1))
  expect_identical(names(p), c("t=60", "t=120", "t=240"))
})

test_that("densemlp interpolates onto arbitrary times and plugs into verbs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("densemlp")

  expect_true("densemlp" %in% list_survlearners()$learner)
  m <- fit(Surv(time, status) ~ age + karno, veteran, "densemlp",
           spec = list(n_bins = 8, epochs = 15))
  expect_s3_class(m, "survalis_fit")
  S <- predict(m, veteran[1:4, ], times = c(45, 137, 303))
  expect_s3_class(S, "survmat")
  expect_equal(attr(S, "times"), c(45, 137, 303))
})

test_that("predict_densemlp() rejects a cox-loss fit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("densemlp")
  mod <- fit_densemlp(Surv(time, status) ~ age, veteran, loss = "cox", epochs = 10)
  expect_error(predict_densemlp(mod, veteran[1:2, ], times = 90), "loss = 'brier'")
})
