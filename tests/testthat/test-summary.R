
test_that("summary.mlsurv_model() returns the expected structure and values", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("cli")

  df <- survival::veteran
  fml <- survival::Surv(time, status) ~ age + trt + celltype

  mod <- list(
    learner = "coxph",
    formula = fml,
    data = df
  )
  class(mod) <- "mlsurv_model"
  attr(mod, "engine") <- "survival"

  res <- summary(mod) 
  expect_s3_class(res, "summary.mlsurv_model")

  expect_identical(res$learner, "coxph")
  expect_identical(res$engine, "survival")
  expect_identical(res$formula, fml)

  ds <- res$data_summary
  expect_type(ds, "list")
  expect_named(ds, c("observations", "predictors", "time_range", "event_rate"), ignore.order = TRUE)

  expect_identical(ds$observations, nrow(df))

  expect_true("age" %in% ds$predictors)
  expect_true("trt" %in% ds$predictors)
  expect_true(any(grepl("^celltype", ds$predictors)))

  expect_length(ds$time_range, 2)
  expect_true(is.numeric(ds$time_range))
  expect_identical(ds$time_range, range(df$time, na.rm = TRUE))

  mf <- model.frame(fml, df)
  y  <- model.response(mf)  
  exp_event_rate <- mean(y[, "status"], na.rm = TRUE)
  expect_equal(ds$event_rate, exp_event_rate, tolerance = 1e-12)

  expect_no_error(summary(mod))
})

test_that("summary.mlsurv_model() errors on wrong class", {
  bad <- list(formula = y ~ x, data = data.frame())

  meth <- getS3method("summary", "mlsurv_model")

  expect_error(meth(bad), "inherits")   
})


test_that("summary() S3 dispatch works for class mlsurv_model", {
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("cli")

  df <- survival::veteran
  fml <- survival::Surv(time, status) ~ age + trt + celltype
  mod <- list(learner = "coxph", formula = fml, data = df)
  class(mod) <- "mlsurv_model"
  attr(mod, "engine") <- "survival"

  s <- summary(mod)
  expect_s3_class(s, "summary.mlsurv_model")
})