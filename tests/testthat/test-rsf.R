test_that("fit_rsf() returns a proper mlsurv_model wrapper", {
  skip_on_cran()
  skip_if_not_installed("randomForestSRC")
  skip_if_not_installed("survival")

  df <- survival::veteran
  Surv <- survival::Surv  

  mod <- fit_rsf(
    Surv(time, status) ~ age + celltype + karno,  
    data = df,
    ntree = 200, mtry = 2, nodesize = 5
  )

  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "rsf")
  expect_identical(attr(mod, "engine"), "randomForestSRC")
  expect_true(!is.null(mod$model))
})

test_that("fit_rsf() accepts a survival::Surv() outcome and records time/status", {
  skip_on_cran()
  skip_if_not_installed("randomForestSRC")
  mod <- fit_rsf(survival::Surv(time, status) ~ age + karno, data = survival::veteran,
                 ntree = 50)
  expect_identical(mod$time, "time")
  expect_identical(mod$status, "status")
  p <- predict_rsf(mod, survival::veteran[1:3, ], times = c(50, 100))
  expect_equal(dim(p), c(3L, 2L))
})
