test_that(".finalize_survmat() clamps and enforces monotonic survival over time", {
  finalize <- getFromNamespace(".finalize_survmat", "survalis")

  raw <- matrix(
    c(1.1, 0.8, 0.85,
      0.9, -0.1, 0.4),
    nrow = 2,
    byrow = TRUE
  )

  out <- finalize(raw, times = c(100, 50, 200))
  M <- as.matrix(out)

  expect_identical(colnames(out), c("t=100", "t=50", "t=200"))
  expect_true(all(M >= 0 & M <= 1))

  ord <- order(c(100, 50, 200))
  dec_ok <- apply(M[, ord, drop = FALSE], 1, function(r) all(diff(r) <= 1e-8))
  expect_true(all(dec_ok))
})

test_that("score_survmodel() standardizes prediction matrices before scoring", {
  local_predict_dummy <- function(object, newdata, times) {
    data.frame(`t=100` = rep(1.2, nrow(newdata)), `t=200` = rep(0.4, nrow(newdata)))
  }

  assign("predict_dummy", local_predict_dummy, envir = .GlobalEnv)
  on.exit(rm("predict_dummy", envir = .GlobalEnv), add = TRUE)

  df <- survival::veteran[, c("time", "status", "age")]
  mod <- structure(
    list(
      learner = "dummy",
      formula = survival::Surv(time, status) ~ age,
      data = df
    ),
    class = "mlsurv_model"
  )

  res <- score_survmodel(mod, times = c(100, 200), metrics = c("cindex", "auc", "ibs"))

  expect_s3_class(res, "data.frame")
  expect_setequal(res$metric, c("cindex", "auc", "ibs"))
  expect_true(all(is.finite(res$value)))
})

test_that("plot_survmat() returns ggplot objects for individual and grouped curves", {
  S <- data.frame(
    `t=1` = c(0.95, 0.90, 0.92),
    `t=2` = c(0.80, 0.70, 0.78),
    `t=3` = c(0.60, 0.45, 0.55),
    check.names = FALSE
  )

  expect_s3_class(plot_survmat(S), "ggplot")
  expect_s3_class(
    plot_survmat(S, group = c("A", "B", "A"), show_individual = TRUE),
    "ggplot"
  )

  expect_error(
    plot_survmat(S, group = c("A", "B")),
    "length equal to nrow"
  )
})

test_that("survmat_to_rmst() interpolates correctly for off-grid tau (regression test for prior under-integration bug)", {
  times <- 0:5
  S <- matrix(c(1, 0.9, 0.8, 0.7, 0.6, 0.5), nrow = 1)
  colnames(S) <- paste0("t=", times)

  # On-grid tau: exact trapezoidal sum over whole intervals [0,1],[1,2].
  expect_equal(survmat_to_rmst(S, times, tau = 2), 1.8, tolerance = 1e-8)

  # Off-grid tau: must interpolate the fractional interval [3,3.5], not
  # silently truncate to the tau = 3 value.
  rmst_tau3   <- survmat_to_rmst(S, times, tau = 3)
  rmst_tau3.5 <- survmat_to_rmst(S, times, tau = 3.5)
  expect_equal(rmst_tau3, 2.55, tolerance = 1e-8)
  expect_equal(rmst_tau3.5, 2.875, tolerance = 1e-8)
  expect_gt(rmst_tau3.5, rmst_tau3)

  # tau beyond the last grid point holds at the full definite integral
  # (no extrapolation past the observed curve).
  expect_equal(survmat_to_rmst(S, times, tau = 10), 3.75, tolerance = 1e-8)
})
