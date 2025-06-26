fit_ranger <- function(formula, data, ...) {
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  model <- ranger::ranger(
    formula = formula,
    data = data,
    ...,
    respect.unordered.factors = "order",
    case.weights = NULL
  )

  structure(list(
    model = model,
    learner = "ranger",
    formula = formula,
    data = data
  ), class = "mlsurv_model", engine = "ranger")
}


predict_ranger <- function(model, newdata, times) {
  requireNamespace("ranger")

  pred_obj <- predict(model$model, data = newdata)
  pred <- pred_obj$survival
  model_times <- model$model$unique.death.times

  if (is.null(pred)) stop("Prediction returned NULL survival matrix.")
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)  # ensure matrix for single row

  # preallocate output
  surv_mat <- matrix(NA_real_, nrow = nrow(pred), ncol = length(times))
  for (i in seq_len(nrow(pred))) {
    y <- pred[i, ]
    surv_mat[i, ] <- stats::approx(
      x = model_times, y = y, xout = times,
      method = "linear", rule = 2, ties = "ordered"
    )$y
  }

  colnames(surv_mat) <- paste0("t=", times)
  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(surv_mat)))
  return(as.data.frame(surv_mat))
}


# example test
library(survival)
data(veteran, package = "survival")

mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred <- predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))

pred

#----------------- add tuner

tune_ranger <- function(formula, data, times,
                        param_grid = expand.grid(
                          num.trees = c(100, 300),
                          mtry = c(1, 2, 3),
                          min.node.size = c(3, 5)
                        ),
                        metrics = c("cindex", "ibs"),
                        folds = 5, seed = 123, ...) {

  purrr::pmap_dfr(param_grid, function(num.trees, mtry, min.node.size) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_ranger,
      pred_fun = predict_ranger,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(num.trees, mtry, min.node.size) |>
      bind_cols(tidyr::pivot_wider(summary[, c("metric", "mean")], names_from = metric, values_from = mean))
  }) |> arrange(dplyr::across(any_of(metrics[1])))
}


res_ranger <- tune_ranger(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    num.trees = c(100, 300),
    mtry = c(1, 2),
    min.node.size = c(3, 5)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3
)

print(res_ranger)

