test_that("predict_xgboost() returns LP when times=NULL and survival when times provided", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("survival")

  set.seed(2025)

  df <- survival::veteran
  df$celltype <- droplevels(df$celltype)

  mod <- fit_xgboost(
    Surv(time, status) ~ age + karno + celltype,
    data    = df,
    nrounds = 20 
  )

  newx <- df[1:5, c("age", "karno", "celltype", "time", "status")]

  lp <- predict_xgboost(mod, newdata = newx, times = NULL)
  expect_type(lp, "double")
  expect_length(lp, nrow(newx))
  expect_true(all(is.finite(lp)))

  times <- c(100, 200, 300)
  psurv <- predict_xgboost(mod, newdata = newx, times = times)

  expect_s3_class(psurv, "data.frame")
  expect_equal(nrow(psurv), nrow(newx))
  expect_equal(ncol(psurv), length(times))
  expect_setequal(colnames(psurv), paste0("t=", times))

  M <- as.matrix(psurv)
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)

  ord <- order(times)
  diffs <- t(apply(M[, ord, drop = FALSE], 1, diff))
  expect_true(all(diffs <= 1e-8))
})
