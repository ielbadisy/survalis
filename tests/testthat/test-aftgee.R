test_that("fit_aftgee() fits and returns expected structure", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("aftgee")
  testthat::skip_if_not_installed("survival")

  df <- survival::veteran

  mod <- fit_aftgee(survival::Surv(time, status) ~ trt + karno + age,
                    data = df, corstr = "independence")

  expect_true(is.list(mod))
  expect_s3_class(mod, "mlsurv_model")
  expect_identical(mod$learner, "aftgee")
  expect_identical(attr(mod, "engine"), "aftgee")
  expect_true(!is.null(mod$formula))
  expect_true(!is.null(mod$data))
  expect_true(!is.null(mod$model))
  expect_true(any(grepl("aftgee", tolower(class(mod$model)))))
})

test_that("predict_aftgee() returns survival matrix with correct shape, names, bounds", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("aftgee")
  testthat::skip_if_not_installed("survival")

  set.seed(123)
  df <- survival::veteran
  idx <- sample(seq_len(nrow(df)), floor(0.7 * nrow(df)))
  train <- df[idx, ]
  test  <- df[-idx, ]

  mod <- fit_aftgee(survival::Surv(time, status) ~ trt + karno + age,
                    data = train)

  times <- c(20, 60, 120, 200)
  pred <- predict_aftgee(mod, newdata = test, times = times)

  expect_true(is.data.frame(pred))
  expect_equal(nrow(pred), nrow(test))
  expect_equal(ncol(pred), length(times))
  expect_identical(colnames(pred), paste0("t=", times))

  M <- as.matrix(pred)
  expect_true(all(is.finite(M)))
  expect_true(all(M >= 0 & M <= 1))

  expect_true(all(M[, 1] >= M[, ncol(M)] - 1e-6))
})

test_that("predict_aftgee() uses default times when times = NULL", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("aftgee")
  testthat::skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_aftgee(survival::Surv(time, status) ~ trt + karno + age, data = df)

  pred <- predict_aftgee(mod, newdata = df[1:10, ], times = NULL)

  expect_true(is.data.frame(pred))
  expect_equal(ncol(pred), 20)
  expect_true(all(grepl("^t=", colnames(pred))))
})

test_that("predict_aftgee() warns on learner tag mismatch and errors on bad inputs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("aftgee")
  testthat::skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_aftgee(survival::Surv(time, status) ~ trt + karno + age, data = df)

  bad <- mod; bad$learner <- "coxph"
  expect_warning(predict_aftgee(bad, newdata = df[1:5, ], times = c(20, 60)),
                 regexp = "may not come from fit_aftgee")

  bad_new <- subset(df, select = -age)
  expect_error(predict_aftgee(mod, newdata = bad_new, times = c(20, 60)))

  
  expect_error(predict_aftgee(mod, newdata = df[1:5, ], times = c("a", "b")))
})

test_that("predict_aftgee() uses model matrix compatible with coef length", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("aftgee")
  testthat::skip_if_not_installed("survival")

  df  <- survival::veteran
  mod <- fit_aftgee(survival::Surv(time, status) ~ trt + karno + age, data = df)

  beta <- mod$model$coef.res
  X    <- model.matrix(delete.response(terms(mod$formula)), df)

  expect_equal(ncol(X), length(beta), tolerance = 0)
})

test_that("predict_aftgee() uses the Kaplan-Meier baseline of the residuals", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("aftgee")
  df  <- survival::veteran
  mod <- fit_aftgee(survival::Surv(time, status) ~ karno + age, data = df)
  tt  <- c(30, 90, 180)
  S   <- as.matrix(predict_aftgee(mod, newdata = df, times = tt))

  ## S(t | x) = KM of the residuals, evaluated at log(t) - x'beta
  lp <- drop(model.matrix(~ karno + age, df) %*% mod$model$coef.res)
  km <- survival::survfit(survival::Surv(log(df$time) - lp, df$status) ~ 1)
  ref <- vapply(tt, function(t) {
    stats::stepfun(km$time, c(1, km$surv))(log(t) - lp)
  }, numeric(nrow(df)))
  expect_equal(unname(S), unname(ref), tolerance = 1e-12)

  ## on the training data the average curve tracks the marginal KM
  km0 <- summary(survival::survfit(survival::Surv(time, status) ~ 1, df),
                 times = tt)$surv
  expect_lt(max(abs(colMeans(S) - km0)), 0.05)
})

test_that("fit_aftgee() works with missing covariates and default times use the training data", {
  testthat::skip_if_not_installed("aftgee")
  lung <- survival::lung
  m <- fit_aftgee(survival::Surv(time, status) ~ age + ph.ecog, data = lung)
  S <- as.matrix(predict_aftgee(m, lung[1:3, ], times = c(200, 400)))
  expect_true(all(S >= 0 & S <= 1))
  expect_equal(dim(predict_aftgee(m, lung[1:3, c("age", "ph.ecog")])), c(3L, 20L))
})
