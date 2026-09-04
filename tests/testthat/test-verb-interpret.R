f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")
  testthat::skip_if_not_installed("pec")
}

mk <- function() fit(f, veteran, "ranger", spec = list(num.trees = 80), seed = 1)

test_that("interpret() dispatches pdp / ice / ale with a feature", {
  skip_verb()
  m <- mk()
  p <- interpret(m, "pdp", feature = "karno", times = c(90, 180))
  expect_s3_class(p, "survalis_interpret")
  expect_equal(p$feature, "karno")
  expect_s3_class(plot(p), "ggplot")

  i <- interpret(m, "ice", feature = "karno", times = c(90, 180))
  expect_s3_class(i, "survalis_interpret")

  a <- interpret(m, "ale", feature = "age", times = c(90, 180), newdata = veteran[1:40, ])
  expect_s3_class(a, "survalis_interpret")
})

test_that("interpret() dispatches surrogate/tree_surrogate/interaction/counterfactual/calibration", {
  skip_verb()
  m <- mk()
  s <- interpret(m, "surrogate", newdata = veteran[1, , drop = FALSE],
                 times = c(90, 180), target_time = 180, k = 5)
  expect_s3_class(s, "survalis_interpret")

  ts <- interpret(m, "tree_surrogate", times = c(90, 180))
  expect_s3_class(ts, "survalis_interpret")

  ia <- interpret(m, "interaction", times = c(90, 180), target_time = 180, type = "1way")
  expect_s3_class(ia, "survalis_interpret")

  cf <- interpret(m, "counterfactual", newdata = veteran[1, , drop = FALSE],
                  times = c(90, 180), target_time = 180)
  expect_s3_class(cf, "survalis_interpret")

  cal <- interpret(m, "calibration", eval_time = 90, n_boot = 5)
  expect_s3_class(cal, "survalis_interpret")
})

test_that("interpret() dispatches varimp / permute alias and shap", {
  skip_verb()
  m <- mk()
  v <- interpret(m, "varimp", times = c(90, 180))
  expect_s3_class(v, "survalis_interpret")
  expect_s3_class(as.data.frame(v), "data.frame")
  expect_identical(interpret(m, "permute", times = c(90, 180))$plot_fn, v$plot_fn)

  s <- interpret(m, "shap", newdata = veteran[5, ], times = c(90, 180),
                 sample.size = 5)
  expect_s3_class(s, "survalis_interpret")
})

test_that("interpret() shap_method reaches compute_shap()'s own aggregation method", {
  skip_verb()
  m <- mk()
  agg <- interpret(m, "shap", newdata = veteran[5, ], times = c(90, 180),
                   sample.size = 5, aggregate = TRUE, shap_method = "integral")
  expect_s3_class(agg, "survalis_interpret")
  expect_identical(attr(as.data.frame(agg), "shap_method"), "integral")
})

test_that("interpret() enforces method-specific required args", {
  skip_verb()
  m <- mk()
  expect_error(interpret(m, "pdp", times = 90), "requires `feature`")
  expect_error(interpret(m, "counterfactual", newdata = veteran[1:3, ], times = c(90, 180)),
               "requires `target_time`")
  expect_error(interpret(m, "calibration", times = c(90, 180)), "single time")
  expect_error(interpret(m, "nope"), "should be one of")
})

test_that("interpret() accepts a survalis_tune by using its refit fit", {
  skip_verb()
  g <- expand.grid(num.trees = 80, mtry = 2, min.node.size = c(5, 15))
  t <- tune(f, veteran, "ranger", grid = g, times = c(90, 180),
            resampling = cv(2, seed = 1), metric = "ibs")
  expect_s3_class(interpret(t, "varimp", times = c(90, 180)), "survalis_interpret")
})
