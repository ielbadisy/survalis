test_that("core interpretability methods return expected structures", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("gower")
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")

  df <- survival::veteran[, c("time", "status", "age", "karno", "celltype", "trt")]
  mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype + trt, data = df)
  times <- c(50, 100, 150)

  shap_td <- compute_shap(mod, df[1, , drop = FALSE], df, times = times, sample.size = 5)
  expect_s3_class(shap_td, "data.frame")
  expect_true(all(c("feature", "phi", "time") %in% names(shap_td)))

  shap_ag <- compute_shap(
    mod, df[1, , drop = FALSE], df,
    times = times, sample.size = 5, aggregate = TRUE
  )
  expect_s3_class(shap_ag, "data.frame")
  expect_true(all(c("feature", "phi") %in% names(shap_ag)))

  pdp <- compute_pdp(mod, df, feature = "age", times = times, method = "pdp+ice", grid.size = 5)
  expect_type(pdp, "list")
  expect_true(all(c("results", "pdp_integrated") %in% names(pdp)))
  expect_s3_class(pdp$results, "data.frame")

  ale <- compute_ale(mod, df, feature = "age", times = times, grid.size = 5)
  expect_type(ale, "list")
  expect_true(all(c("ale", "integrated") %in% names(ale)))
  expect_s3_class(ale$ale, "data.frame")

  surrogate <- compute_surrogate(
    mod, df[1, , drop = FALSE], df,
    times = times, target_time = 100, k = 3, penalized = FALSE
  )
  expect_s3_class(surrogate, "data.frame")
  expect_true(all(c("feature", "feature_value", "effect", "target_time") %in% names(surrogate)))

  tree_surrogate <- compute_tree_surrogate(mod, df, times = c(100, 150), minsplit = 20, cp = 0.05)
  expect_s3_class(tree_surrogate, "tree_surrogate")
  expect_true(all(c("times", "results", "dynamic") %in% names(tree_surrogate)))

  counterfactual <- compute_counterfactual(
    mod, df[1, , drop = FALSE],
    times = times, target_time = 100,
    features_to_change = c("age", "karno", "trt"),
    grid.size = 10
  )
  expect_s3_class(counterfactual, "data.frame")
  expect_true(all(
    c("feature", "original_value", "suggested_value", "survival_gain", "change_cost", "penalized_gain") %in%
      names(counterfactual)
  ))

  interp_methods <- list_interpretability_methods()
  cf_row <- interp_methods[interp_methods$compute == "compute_counterfactual", , drop = FALSE]
  expect_identical(cf_row$plot, "plot_counterfactual")
  expect_true(cf_row$has_plot)
})

test_that("interaction, variable importance, and calibration helpers run", {
  skip_on_cran()
  skip_if_not_installed("survival")

  df <- survival::veteran[, c("time", "status", "age", "karno", "celltype", "trt")]
  mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype + trt, data = df)
  times <- c(50, 100, 150)
  features <- c("age", "karno", "trt")

  inter_1way <- compute_interactions(mod, df, times = times, target_time = 100, features = features, type = "1way", grid.size = 5)
  expect_s3_class(inter_1way, "data.frame")
  expect_true(all(c("feature", "interaction") %in% names(inter_1way)))

  inter_heat <- compute_interactions(mod, df, times = times, target_time = 100, features = features, type = "heatmap", grid.size = 5)
  expect_s3_class(inter_heat, "data.frame")
  expect_true(all(c("feature1", "feature2", "interaction") %in% names(inter_heat)))

  inter_time <- compute_interactions(mod, df, times = times, features = features, type = "time", grid.size = 5)
  expect_s3_class(inter_time, "data.frame")
  expect_true(all(c("feature", "time", "interaction") %in% names(inter_time)))

  varimp <- compute_varimp(mod, times = 100, metric = "brier", n_repetitions = 2, seed = 1)
  expect_s3_class(varimp, "data.frame")
  expect_true(all(c("feature", "importance", "importance_05", "importance_95", "scaled_importance") %in% names(varimp)))

  calibration <- compute_calibration(
    mod, data = df, time = "time", status = "status",
    eval_time = 100, n_bins = 5, n_boot = 5, seed = 1
  )
  expect_type(calibration, "list")
  expect_true("calibration_table" %in% names(calibration))
  expect_s3_class(calibration$calibration_table, "data.frame")
})

test_that("interpretability helpers satisfy core numeric invariants", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("gower")
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")

  df <- survival::veteran[, c("time", "status", "age", "karno", "celltype", "trt")]
  mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype + trt, data = df)
  times <- c(50, 100, 150)
  features <- c("age", "karno", "trt")

  shap_ag <- compute_shap(
    mod, df[1, , drop = FALSE], df,
    times = times, sample.size = 5, aggregate = TRUE, method = "meanabs"
  )
  expect_true(all(shap_ag$phi >= 0))

  ale <- compute_ale(mod, df, feature = "age", times = times, grid.size = 5)
  ale_mat <- as.matrix(ale$ale[, -1, drop = FALSE])
  expect_true(all(is.finite(ale_mat)))
  expect_true(all(abs(colMeans(ale_mat)) < 1e-8))

  surrogate_pen <- compute_surrogate(
    mod, df[1, , drop = FALSE], df,
    times = times, target_time = 100, k = 2, penalized = TRUE
  )
  expect_lte(nrow(surrogate_pen), 2)
  expect_true(all(diff(abs(surrogate_pen$effect)) <= 1e-8))

  tree_surrogate <- compute_tree_surrogate(mod, df, times = c(100, 150), minsplit = 20, cp = 0.05)
  expect_length(tree_surrogate$results, 2)
  expect_true(tree_surrogate$dynamic)
  expect_true(all(vapply(tree_surrogate$results, function(x) is.finite(x$r_squared), logical(1))))

  inter_heat <- compute_interactions(mod, df, times = times, target_time = 100, features = features, type = "heatmap", grid.size = 5)
  diag_rows <- inter_heat[inter_heat$feature1 == inter_heat$feature2, ]
  expect_true(all(diag_rows$interaction == 0))
  for (i in seq_len(nrow(inter_heat))) {
    rev_val <- inter_heat$interaction[
      inter_heat$feature1 == inter_heat$feature2[i] &
        inter_heat$feature2 == inter_heat$feature1[i]
    ]
    expect_equal(inter_heat$interaction[i], rev_val, tolerance = 1e-12)
  }

  inter_time <- compute_interactions(mod, df, times = times, features = features, type = "time", grid.size = 5)
  expect_setequal(unique(inter_time$time), times)
  expect_true(all(inter_time$interaction >= 0))
  expect_true(all(inter_time$interaction <= 1))

  calibration <- compute_calibration(
    mod, data = df, time = "time", status = "status",
    eval_time = 100, n_bins = 5, n_boot = 5, seed = 1
  )
  ct <- calibration$calibration_table
  expect_true(all(ct$mean_pred_surv >= 0 & ct$mean_pred_surv <= 1))
  expect_true(all(ct$observed_surv >= 0 & ct$observed_surv <= 1))
})

test_that("interpretability helpers error cleanly on unsupported inputs", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("gower")
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")

  df <- survival::veteran[, c("time", "status", "age", "karno", "celltype", "trt")]
  mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype + trt, data = df)
  times <- c(50, 100, 150)

  expect_error(
    compute_counterfactual(mod, df[1, , drop = FALSE], times = times, target_time = 100, features_to_change = "not_a_feature"),
    "No valid features to change"
  )

  expect_error(
    compute_interactions(mod, df, times = times, type = "1way"),
    "target_time must be specified"
  )

  expect_error(
    compute_ale(mod, df, feature = "celltype", times = times, grid.size = 5),
    "only defined for numeric continuous features"
  )

  expect_error(
    compute_surrogate(mod, df[1, , drop = FALSE], df, times = times, target_time = 100, dist.fun = "euclidean"),
    "Specify kernel.width"
  )
})

test_that("interpretability plotting helpers return plot objects", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("gower")
  skip_if_not_installed("rpart")
  skip_if_not_installed("partykit")

  df <- survival::veteran[, c("time", "status", "age", "karno", "celltype", "trt")]
  mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype + trt, data = df)
  times <- c(50, 100, 150)

  shap_td <- compute_shap(mod, df[1, , drop = FALSE], df, times = times, sample.size = 5)
  shap_ag <- compute_shap(mod, df[1, , drop = FALSE], df, times = times, sample.size = 5, aggregate = TRUE)
  pdp <- compute_pdp(mod, df, feature = "age", times = times, method = "pdp+ice", grid.size = 5)
  ale <- compute_ale(mod, df, feature = "age", times = times, grid.size = 5)
  surrogate <- compute_surrogate(
    mod, df[1, , drop = FALSE], df,
    times = times, target_time = 100, k = 3, penalized = FALSE
  )
  counterfactual <- compute_counterfactual(
    mod, df[1, , drop = FALSE],
    times = times, target_time = 100,
    features_to_change = c("age", "karno", "trt"),
    grid.size = 10
  )
  tree_surrogate <- compute_tree_surrogate(mod, df, times = c(100, 150), minsplit = 20, cp = 0.05)
  inter_time <- compute_interactions(mod, df, times = times, features = c("age", "karno", "trt"), type = "time", grid.size = 5)
  varimp <- compute_varimp(mod, times = 100, metric = "brier", n_repetitions = 2, seed = 1)
  calibration <- compute_calibration(
    mod, data = df, time = "time", status = "status",
    eval_time = 100, n_bins = 5, n_boot = 5, seed = 1
  )

  expect_s3_class(plot_shap(shap_td), "ggplot")
  expect_s3_class(plot_shap(shap_ag), "ggplot")
  expect_s3_class(plot_pdp(pdp, feature = "age", which = "per_time"), "ggplot")
  expect_s3_class(plot_pdp(pdp, feature = "age", which = "integrated"), "ggplot")
  expect_s3_class(plot_ale(ale, feature = "age", which = "per_time"), "ggplot")
  expect_s3_class(plot_ale(ale, feature = "age", which = "integrated"), "ggplot")
  expect_s3_class(plot_counterfactual(counterfactual), "ggplot")
  expect_s3_class(plot_surrogate(surrogate), "ggplot")
  expect_s3_class(plot_tree_surrogate(tree_surrogate, type = "importance"), "ggplot")
  expect_s3_class(plot_interactions(inter_time, type = "time"), "ggplot")
  expect_s3_class(plot_varimp(varimp), "ggplot")
  expect_s3_class(plot_calibration(calibration), "ggplot")
})
