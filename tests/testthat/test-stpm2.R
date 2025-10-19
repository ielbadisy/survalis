
test_that("fit_stpm2() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("rstpm2")
  skip_if_not_installed("survival")

  df <- survival::veteran
  Surv <- survival::Surv

  mod <- fit_stpm2(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    df = 4
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "stpm2")
  expect_identical(mod$engine, "rstpm2")
  expect_true(!is.null(mod$model))
})

test_that("predict_stpm2() returns bounded, monotone survival at requested times", {
  skip_on_cran()
  skip_if_not_installed("rstpm2")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  Surv <- survival::Surv
  mod <- fit_stpm2(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    df = 3
  )

  times <- c(150, 30, 300, 200)
  newx  <- df[1:6, ]

  pred <- predict_stpm2(mod, newdata = newx, times = times)

  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), nrow(newx))
  expect_equal(ncol(pred), length(times))
  expect_setequal(colnames(pred), paste0("t=", as.character(times)))

  M <- as.matrix(pred)
  expect_true(is.numeric(M))
  expect_true(all(is.finite(M)))
  expect_gte(min(M), -1e-12)
  expect_lte(max(M), 1 + 1e-12)

  ord <- order(times)
  M_ord <- M[, paste0("t=", as.character(times[ord])), drop = FALSE]
  diffs <- t(apply(M_ord, 1, diff))
  expect_true(all(diffs <= 1e-8))
})

test_that("predict_stpm2() errors when newdata lacks required predictors", {
  skip_on_cran()
  skip_if_not_installed("rstpm2")
  skip_if_not_installed("survival")

  df  <- survival::veteran
  Surv <- survival::Surv
  mod <- fit_stpm2(
    Surv(time, status) ~ age + karno + celltype,
    data = df,
    df = 3
  )

  bad_new <- df[1:3, c("age", "karno")]  
  expect_error(
    predict_stpm2(mod, newdata = bad_new, times = c(50, 100)),
    regexp = paste0(
      "object .* not found|undefined columns|subset.*select|",
      "do not match|newdata.*missing|terms|model\\.matrix|",
      'Value of covariate ".*" not supplied in "newdata"|',
      "names of independing variables.*do not match"
    )
  )
})
