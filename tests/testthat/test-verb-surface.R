f <- Surv(time, status) ~ age + karno

expo <- list(
  name = "exponential",
  fit = function(formula, data, ...) survival::survreg(formula, data, dist = "exponential"),
  predict = function(object, newdata, times) {
    rate <- 1 / predict(object, newdata, type = "response")
    exp(-outer(rate, times))
  }
)

test_that("the exported surface is the seven verbs plus helpers and dwscore", {
  ex <- sort(getNamespaceExports("survalis"))
  expect_setequal(ex, c(
    "fit", "evaluate", "compare", "tune", "interpret", "estimate",
    "cv", "holdout", "group_cv", "bootstrap",
    "list_survlearners", "list_metrics", "list_interpretability_methods",
    "dwscore", "Surv"
  ))
})

test_that("list_survlearners() returns ids and tunable flags", {
  ls <- list_survlearners()
  expect_named(ls, c("learner", "tunable"))
  expect_true("coxph" %in% ls$learner)
  expect_true(all(list_survlearners(has_tune = TRUE)$tunable))
})

test_that("a custom learner works with fit() and predict()", {
  m <- fit(f, veteran, model = expo)
  expect_s3_class(m, "survalis_fit")
  expect_identical(m$model_id, "exponential")
  S <- predict(m, veteran[1:4, ], times = c(30, 90))
  expect_equal(dim(S), c(4L, 2L))
  expect_true(all(S >= 0 & S <= 1))
  expect_length(predict(m, veteran[1:4, ], times = c(30, 90), type = "rmst"), 4L)
})

test_that("a custom learner with a malformed predict() errors clearly", {
  bad <- expo
  bad$predict <- function(object, newdata, times) matrix(0.5, nrow(newdata), 1)
  m <- fit(f, veteran, model = bad)
  expect_error(predict(m, veteran[1:3, ], times = c(30, 90)), "2 columns")
  expect_error(fit(f, veteran, model = list(fit = 1, predict = 2)), "custom learner")
  expect_error(fit(f, veteran, model = "nope"), "Unknown model")
})

test_that("a custom learner works with evaluate(), compare() and estimate()", {
  ev <- evaluate(f, veteran, model = expo, times = c(30, 90, 180),
                 resampling = cv(3, seed = 1))
  expect_s3_class(ev, "survalis_eval")
  expect_identical(ev$model, "exponential")

  m <- fit(f, veteran, model = expo)
  ev2 <- evaluate(m, times = c(30, 90, 180), resampling = cv(3, seed = 1))
  expect_equal(ev$summary$mean, ev2$summary$mean)

  cmp <- compare(f, veteran, models = list("coxph", expo), times = c(30, 90, 180),
                 resampling = cv(3, seed = 1))
  expect_setequal(cmp$models, c("coxph", "exponential"))

  d <- veteran
  d$trt <- as.integer(d$trt == 2)
  est <- estimate(Surv(time, status) ~ trt + age + karno, d, model = expo,
                  treatment = "trt", tau = 180, times = c(30, 90, 180), R = 5, seed = 1)
  expect_s3_class(est, "survalis_estimate")
})

test_that("compare() rejects duplicate labels", {
  expect_error(compare(f, veteran, models = list(expo, expo), times = 90), "unique")
})

test_that("tune() runs a supplied grid for learners without a built-in tuner", {
  shape <- list(
    name = "weibull",
    fit = function(formula, data, dist = "weibull", ...) survival::survreg(formula, data, dist = dist),
    predict = function(object, newdata, times) {
      lp <- predict(object, newdata, type = "lp")
      if (object$dist == "exponential") {
        return(exp(-outer(exp(-lp), times)))
      }
      scale <- object$scale
      exp(-outer(exp(-lp / scale), times^(1 / scale)))
    }
  )
  tn <- tune(f, veteran, model = shape,
             grid = data.frame(dist = c("weibull", "exponential")),
             times = c(30, 90, 180), resampling = cv(3, seed = 1), metric = "ibs")
  expect_s3_class(tn, "survalis_tune")
  expect_equal(nrow(tn$results), 2L)
  expect_true(tn$fit$spec$dist %in% c("weibull", "exponential"))
  expect_false(is.unsorted(tn$results$ibs))
  expect_equal(dim(predict(tn, veteran[1:2, ], times = c(30, 90))), c(2L, 2L))
})

test_that("evaluate() scores a fit on new data and a supplied matrix", {
  idx <- 1:100
  m <- fit(f, veteran[idx, ], "coxph")
  ev <- evaluate(m, newdata = veteran[-idx, ], times = c(30, 90, 180))
  expect_s3_class(ev, "survalis_eval")
  expect_null(ev$resampling)
  expect_output(print(ev), "scored once on 37 subjects")

  S <- predict(m, veteran[-idx, ], times = c(30, 90, 180))
  y <- Surv(veteran$time[-idx], veteran$status[-idx])
  ev2 <- evaluate(S, y = y)
  expect_equal(ev2$summary$mean, ev$summary$mean)

  P <- unclass(as.matrix(S))
  attributes(P) <- list(dim = dim(P))
  ev3 <- evaluate(P, y = y, times = c(30, 90, 180))
  expect_equal(ev3$summary$mean, ev$summary$mean)
  expect_error(evaluate(P, y = y), "times")
  expect_error(evaluate(P[-1, ], y = y, times = c(30, 90, 180)), "rows")
})

test_that("interpret() draws ROC and decision curves", {
  m <- fit(f, veteran, "coxph")
  r <- interpret(m, "roc", eval_time = 90)
  expect_s3_class(r, "survalis_interpret")
  expect_s3_class(plot(r), "ggplot")
  d <- interpret(m, "dca", eval_time = 90, thresholds = seq(0.1, 0.9, 0.1))
  expect_s3_class(plot(d), "ggplot")
  expect_error(interpret(m, "roc", times = c(30, 90)), "single time")
})
