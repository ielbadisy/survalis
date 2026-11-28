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

  expect_true("fastgbm" %in% listlearners()$learner)
  m <- fit(Surv(time, status) ~ age + karno, veteran, "fastgbm",
           spec = list(ntrees = 40, pexp_bins = 8))
  expect_s3_class(m, "survalis_fit")
  S <- predict(m, veteran[1:4, ], times = c(90, 180))
  expect_s3_class(S, "survmat")
})

test_that("predict_fastgbm() works for the cox and aft objectives", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fastgbm")
  tt <- c(30, 90, 180)
  for (o in c("cox", "aft")) {
    mod <- fit_fastgbm(Surv(time, status) ~ age + karno, veteran, objective = o, ntrees = 20)
    S <- as.matrix(predict_fastgbm(mod, veteran[1:4, ], times = tt))
    expect_equal(dim(S), c(4L, 3L))
    # each column equals the engine's single-time prediction for those subjects
    for (j in seq_along(tt)) {
      ref <- predict(mod$model, newdata = veteran[1:4, ], type = "survival", times = tt[j])
      expect_equal(unname(S[, j]), as.vector(ref), tolerance = 1e-8)
    }
    expect_true(all(diff(t(S)) <= 0))
  }
})

test_that("predict_fastgbm() works without outcome columns and for one row", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fastgbm")

  mod <- fit_fastgbm(Surv(time, status) ~ age + karno + trt, veteran, ntrees = 30)
  tt <- c(200, 10, 100)
  b <- as.matrix(predict_fastgbm(mod, veteran, times = tt))
  expect_identical(colnames(b), paste0("t=", tt))
  nd <- veteran[1:5, c("age", "karno", "trt")]
  expect_equal(unname(as.matrix(predict_fastgbm(mod, nd, times = tt))), unname(b[1:5, ]))
  one <- as.matrix(predict_fastgbm(mod, nd[1, , drop = FALSE], times = tt))
  expect_equal(unname(one[1, ]), unname(b[1, ]))
})

test_that("tune() searches a default fastgbm grid", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("fastgbm")
  expect_true("fastgbm" %in% listlearners(has_tune = TRUE)$learner)
  g <- eval(formals(tune_fastgbm)$param_grid)
  expect_true(all(names(g) %in% names(formals(fastgbm::fastgbm))))
  small <- expand.grid(ntrees = c(20, 40), learning_rate = 0.1, max_depth = 2, min_node_size = 10)
  tu <- tune(Surv(time, status) ~ ., veteran, "fastgbm", grid = small, times = c(60, 120),
             resampling = cv(3, seed = 1), objective = "cox")
  expect_s3_class(tu, "survalis_tune")
  expect_equal(nrow(tu$results), 2L)
  expect_true(all(diff(tu$results$ibs) >= 0))
  expect_identical(tu$fit$spec$objective, "cox")
  expect_s3_class(predict(tu, veteran[1:3, ], times = c(60, 120)), "survmat")
})
