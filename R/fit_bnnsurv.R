
fit_bnnsurv <- function(formula, data,
                        k = 2,
                        num_base_learners = 50,
                        num_features_per_base_learner = NULL,
                        metric = "mahalanobis",
                        weighting_function = function(x) x * 0 + 1,
                        replace = TRUE,
                        sample_fraction = NULL,
                        ...) {
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  model <- bnnSurvival::bnnSurvival(
    formula = formula,
    data = data,
    k = k,
    num_base_learners = num_base_learners,
    num_features_per_base_learner = num_features_per_base_learner,
    metric = metric,
    weighting_function = weighting_function,
    replace = replace,
    sample_fraction = sample_fraction,
    ...
  )

  structure(list(
    model = model,
    learner = "bnnsurv",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "bnnsurv")
}


predict_bnnsurv <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "bnnsurv")
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  # predict full survival curves
  pred <- predict(object$model, newdata)

  # extract prediction grid
  pred_times <- bnnSurvival::timepoints(pred)
  pred_matrix <- bnnSurvival::predictions(pred)

  surv_probs <- t(apply(pred_matrix, 1, function(row) {
    approx(x = pred_times, y = row, xout = times, method = "linear", rule = 2)$y
  }))

  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))

  return(as.data.frame(surv_probs))
}


library(survival)
library(bnnSurvival)

mod_bnnsurv <- fit_bnnsurv(Surv(time, status) ~ age + karno + diagtime + prior, data = veteran)

pred_bnnsurv <- predict_bnnsurv(mod_bnnsurv, newdata = veteran[1:5, ], times = c(1, 100, 200, 300))
print(round(pred_bnnsurv, 3))


#-------- add tuner()

#' Tune bnnSurvival learner
#'
#' @param formula Survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data Training dataset
#' @param times Time points at which to evaluate predictions
#' @param param_grid A data.frame of hyperparameter combinations
#' @param metrics Metrics to evaluate (e.g., "cindex", "ibs")
#' @param folds Number of CV folds
#' @param seed Random seed
#' @param ... Additional args passed to fit_bnnsurv()
#'
#' @return Tibble summarizing average performance per parameter setting
#' @export
tune_bnnsurv <- function(formula, data, times,
                         param_grid = expand.grid(
                           k = c(2, 5),
                           num_base_learners = c(30, 50),
                           sample_fraction = c(0.5, 1)
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5,
                         seed = 123, ...) {

  purrr::pmap_dfr(param_grid, function(k, num_base_learners, sample_fraction) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_bnnsurv,
      pred_fun = predict_bnnsurv,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      k = k,
      num_base_learners = num_base_learners,
      sample_fraction = sample_fraction,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      k = k,
      num_base_learners = num_base_learners,
      sample_fraction = sample_fraction
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  }) |> arrange(dplyr::across(any_of(metrics[1])))
}



res_bnnsurv <- tune_bnnsurv(
  formula = Surv(time, status) ~ age + karno + diagtime + prior,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    k = c(2, 20),
    num_base_learners = c(10, 50),
    sample_fraction = c(0.5, 0.75)
  ),
  folds = 3,
  metrics = c("cindex", "ibs")
)

print(res_bnnsurv)
