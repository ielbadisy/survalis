test_that("fit_bart() fits BART survival and returns expected structure", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("BART")
  testthat::skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_bart(
    survival::Surv(time, status) ~ age + karno + celltype,
    data = df,
    K = 3
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "bart")
  expect_identical(attr(mod, "engine"), "BART")
  expect_true(!is.null(mod$model))
  expect_true(!is.null(mod$eval_times))
  expect_type(mod$eval_times, "double")
  expect_true(length(mod$eval_times) >= 2)
})

test_that("predict_bart() returns well-formed survival probabilities and a proper survival matrix", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("BART")
  testthat::skip_if_not_installed("survival")

  set.seed(123)
  df <- survival::veteran
  idx <- sample(seq_len(nrow(df)), floor(0.7 * nrow(df)))
  train <- df[idx, ]
  test  <- df[-idx, ]

  mod <- fit_bart(
    survival::Surv(time, status) ~ age + karno + celltype,
    data = train,
    K = 3
  )

  req_times <- c(10, 30, 60, 61.2)
  pred <- predict_bart(mod, newdata = test, times = req_times)

  expect_true(is.data.frame(pred))
  expect_equal(nrow(pred), nrow(test))
  expect_equal(ncol(pred), length(req_times))
  expect_identical(colnames(pred), paste0("t=", req_times))

  M <- as.matrix(pred)
  expect_true(all(is.finite(M)))
  expect_true(all(M >= 0 & M <= 1))

  ord <- order(req_times)
  Ms  <- M[, ord, drop = FALSE]
  expect_true(all(Ms[, 1] >= Ms[, ncol(Ms)] - 1e-6))

  et <- sort(unique(mod$eval_times))
  t_below <- min(et) - 5
  t_above <- max(et) + 5
  t_mid   <- et[ceiling(length(et) / 2)]

  pred_align <- predict_bart(
    mod, newdata = test,
    times = c(t_below, min(et), t_above, max(et), t_mid, t_mid)
  )

  A <- as.matrix(pred_align)
  expect_true(max(abs(A[, 1] - A[, 2]), na.rm = TRUE) < 1e-6)
  expect_true(max(abs(A[, 3] - A[, 4]), na.rm = TRUE) < 1e-6)
  expect_true(max(abs(A[, 5] - A[, 6]), na.rm = TRUE) < 1e-6)

  pred_one <- predict_bart(mod, newdata = test[1:7, ], times = t_mid)
  expect_true(is.data.frame(pred_one))
  expect_equal(nrow(pred_one), 7)
  expect_equal(ncol(pred_one), 1)
  expect_identical(colnames(pred_one), paste0("t=", t_mid))
  expect_true(all(pred_one[[1]] >= 0 & pred_one[[1]] <= 1))
})

test_that("predict_bart() warns on learner mismatch and errors on bad inputs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("BART")
  testthat::skip_if_not_installed("survival")

  df <- survival::veteran
  mod <- fit_bart(
    survival::Surv(time, status) ~ age + karno + celltype,
    data = df,
    K = 3
  )

  bad <- mod; bad$learner <- "coxph"
  expect_warning(
    predict_bart(bad, newdata = df[1:5, ], times = c(10, 30)),
    regexp = "may not come from fit_bart"
  )

  bad_new <- subset(df, select = -age)
  expect_error(predict_bart(mod, newdata = bad_new, times = c(10, 30)))

  expect_error(predict_bart(mod, newdata = df[1:5, ], times = c("a", "b")))
})
