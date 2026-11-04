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

test_that("as.data.frame() tidies a survcluster result in long format", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("unsurv")
  m <- fit(f, veteran, "coxph")
  tt <- c(30, 90, 180)
  r <- interpret(m, "survcluster", times = tt, K = 2, seed = 1)
  d <- as.data.frame(r)
  expect_s3_class(d, "data.frame")
  expect_named(d, c("id", "cluster", "medoid", "time", "surv"))
  n <- nrow(r$result$survmat)
  expect_equal(nrow(d), n * length(tt))
  expect_setequal(unique(d$time), tt)
  expect_equal(sort(unique(d$cluster)), 1:2)
  expect_equal(sum(d$medoid), 2 * length(tt))
  # one cluster per observation, matching the clustering fit
  per_id <- unique(d[c("id", "cluster")])
  expect_equal(nrow(per_id), n)
  expect_equal(per_id$cluster[order(per_id$id)], r$result$fit$clusters)
  expect_equal(d$surv[d$id == 3 & d$time == 90],
               unname(r$result$survmat[3, 2]))
  expect_true(all(d$surv >= 0 & d$surv <= 1))
})

test_that("as.data.frame() tidies roc, dca and tree_surrogate results", {
  testthat::skip_on_cran()
  m <- fit(f, veteran, "coxph")
  roc <- as.data.frame(interpret(m, "roc", eval_time = 90))
  expect_named(roc, c("eval_time", "threshold", "sensitivity", "specificity", "auc"))
  expect_true(all(roc$eval_time == 90))
  dca <- as.data.frame(interpret(m, "dca", eval_time = 90))
  expect_named(dca, c("eval_time", "threshold", "model", "treat_all", "treat_none"))
  testthat::skip_if_not_installed("partykit")
  ts <- as.data.frame(interpret(m, "tree_surrogate", times = c(90, 180)))
  expect_named(ts, c("time", "feature", "splits", "r_squared"))
  expect_setequal(unique(ts$time), c(90, 180))
})

test_that("as.data.frame() tidies pdp, ice, ale, calibration and varimp results", {
  testthat::skip_on_cran()
  m <- fit(f, veteran, "coxph")
  tt <- c(90, 180)
  p <- as.data.frame(interpret(m, "pdp", feature = "karno", times = tt))
  expect_named(p, c("feature", "feature_value", "id", "time", "surv", "type"))
  expect_true(all(p$feature == "karno"))
  expect_setequal(unique(p$type), c("ice", "pdp"))
  expect_true(all(is.na(p$id[p$type == "pdp"])))
  pi <- interpret(m, "pdp", feature = "karno", times = tt)
  integ <- as.data.frame(pi, what = "integrated")
  expect_named(integ, c("feature", "feature_value", "integrated"))
  expect_equal(nrow(integ), nrow(pi$result$pdp_integrated))

  i <- as.data.frame(interpret(m, "ice", feature = "celltype", times = tt))
  expect_named(i, c("feature", "feature_value", "id", "time", "surv", "type"))
  expect_true(all(i$type == "ice"))

  a <- interpret(m, "ale", feature = "age", times = tt)
  ad <- as.data.frame(a)
  expect_named(ad, c("feature", "feature_value", "time", "ale"))
  expect_named(as.data.frame(a, what = "integrated"),
               c("feature", "feature_value", "integrated"))
  expect_equal(nrow(ad), nrow(a$result$ale) * length(tt))
  expect_setequal(unique(ad$time), tt)
  expect_equal(ad$ale[ad$time == 180], a$result$ale[["t=180"]])

  cal <- as.data.frame(interpret(m, "calibration", eval_time = 90, n_boot = 3))
  expect_named(cal, c("eval_time", "bin", "mean_pred_surv", "observed_surv",
                      "lower_ci", "upper_ci"))
  expect_true(all(cal$eval_time == 90))

  vi <- as.data.frame(interpret(m, "varimp", times = tt))
  expect_identical(class(vi), "data.frame")
})

test_that("every interpret() method prints and has a tidy data frame form", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("unsurv")
  testthat::skip_if_not_installed("partykit")
  m <- fit(f, veteran, "coxph")
  tt <- c(90, 180)
  one <- veteran[1, , drop = FALSE]
  args <- list(
    pdp = list(feature = "karno", times = tt),
    ice = list(feature = "karno", times = tt),
    ale = list(feature = "age", times = tt),
    shap = list(times = tt),
    shap_matrix = list(newdata = veteran[1:4, ], times = tt),
    varimp = list(times = tt),
    permute = list(times = tt),
    interaction = list(times = tt, target_time = 180, type = "1way"),
    counterfactual = list(newdata = one, times = tt, target_time = 90),
    surrogate = list(newdata = one, times = tt, target_time = 90),
    tree_surrogate = list(times = tt),
    calibration = list(eval_time = 90, n_boot = 3),
    roc = list(eval_time = 90),
    dca = list(eval_time = 90),
    survcluster = list(times = tt, K = 2, seed = 1)
  )
  expect_setequal(names(args), list_interpretability_methods()$method)
  for (meth in names(args)) {
    r <- do.call(interpret, c(list(m, meth), args[[meth]]))
    d <- as.data.frame(r)
    expect_identical(class(d), "data.frame", info = meth)
    expect_gt(nrow(d), 0)
    expect_false(any(grepl("[.=]", names(d))), info = meth)
    expect_false(anyDuplicated(names(d)) > 0, info = meth)
    expect_output(print(r), "\\$result:")
  }
  r <- interpret(m, "shap", times = tt)
  expect_error(as.data.frame(r, what = "integrated"), "applies to")
  r$method <- "unknown"
  r$result <- list(a = 1:2, b = 1:3)
  expect_error(as.data.frame(r), "No data frame form")
})
