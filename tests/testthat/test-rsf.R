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
