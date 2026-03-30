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

  expect_s3_class(res, "tbl_df")
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
