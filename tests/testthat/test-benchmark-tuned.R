test_that("benchmark_tuned_survlearners() returns nested CV results and final refits", {
  skip_on_cran()
  skip_if_not_installed("ranger")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("survival")

  Surv <- survival::Surv
  df <- survival::veteran

  res <- benchmark_tuned_survlearners(
    formula = Surv(time, status) ~ age + karno + celltype,
    data = df,
    learners = c("ranger", "glmnet"),
    times = c(90, 180),
    metrics = c("cindex", "ibs"),
    outer_folds = 2,
    inner_folds = 2,
    seed = 2026,
    refit_final = TRUE,
    learner_args = list(
      ranger = list(tune = list(
        param_grid = expand.grid(
          num.trees = c(50, 75),
          mtry = c(1, 2),
          min.node.size = c(3),
          stringsAsFactors = FALSE
        )
      )),
      glmnet = list(tune = list(
        param_grid = list(alpha = c(0, 1))
      ))
    )
  )

  expect_s3_class(res, "nested_surv_benchmark")
  expect_true(is.data.frame(res$outer_results))
  expect_true(is.data.frame(res$outer_summary))
  expect_true(is.data.frame(res$selected_params))

  expect_true(all(c("learner", "id", "outer_fold", "metric", "value") %in% names(res$outer_results)))
  expect_true(all(c("learner", "metric", "mean", "sd", "n", "se", "lower", "upper") %in% names(res$outer_summary)))
  expect_true(all(c("learner", "id", "outer_fold") %in% names(res$selected_params)))

  expect_setequal(unique(res$outer_results$learner), c("ranger", "glmnet"))
  expect_setequal(unique(res$selected_params$learner), c("ranger", "glmnet"))
  expect_setequal(unique(res$outer_results$metric), c("cindex", "ibs"))

  expect_equal(length(res$final_models), 2)
  expect_s3_class(res$final_models$ranger, "mlsurv_model")
  expect_s3_class(res$final_models$glmnet, "mlsurv_model")
  expect_true(!is.null(attr(res$final_models$ranger, "tuning_results")))
  expect_true(!is.null(attr(res$final_models$glmnet, "tuning_results")))
})

test_that("benchmark_tuned_survlearners() supports list-valued tuning parameters", {
  skip_on_cran()
  skip_if_not_installed("survdnn")
  skip_if_not_installed("torch")
  skip_if_not_installed("survival")
  skip_if_not(torch::torch_is_installed(), "Torch backend is not installed")

  Surv <- survival::Surv
  df <- survival::veteran

  res <- benchmark_tuned_survlearners(
    formula = Surv(time, status) ~ age + karno + celltype,
    data = df,
    learners = "survdnn",
    times = 90,
    metrics = c("cindex", "ibs"),
    outer_folds = 2,
    inner_folds = 2,
    seed = 2026,
    learner_args = list(
      survdnn = list(
        tune = list(
          param_grid = list(
            hidden = list(c(8L), c(16L, 8L)),
            lr = c(1e-3),
            activation = c("relu"),
            epochs = c(2L),
            loss = c("cox"),
            optimizer = c("adam"),
            dropout = c(0.1),
            batch_norm = c(TRUE)
          )
        ),
        fit = list(verbose = FALSE)
      )
    )
  )

  expect_s3_class(res, "nested_surv_benchmark")
  expect_true(is.list(res$selected_params$hidden))
  expect_equal(length(res$selected_params$hidden), 2)
  expect_setequal(unique(res$outer_results$learner), "survdnn")
  expect_setequal(unique(res$outer_results$metric), c("cindex", "ibs"))
})
