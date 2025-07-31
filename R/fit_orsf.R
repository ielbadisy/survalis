#' Fit a Survival Forest using aorsf
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data.frame containing the variables
#' @param ... Additional arguments passed to `aorsf::orsf()`
#'
#' @return A list of class 'mlsurv_model' with attributes for unified prediction
#' @export
fit_orsf <- function(formula, data, ...) {
  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  model <- aorsf::orsf(
    formula = formula,
    data = data,
    na_action = "omit",
    ...
  )

  structure(list(
    model = model,
    learner = "aorsf",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "aorsf")
}


#' Predict survival probabilities using ORSF model
#'
#' @param object An object from `fit_orsf()`
#' @param newdata New data for prediction
#' @param times Time points at which to predict survival probabilities
#' @param ... Additional arguments passed to `aorsf::predict()`
#'
#' @return A data.frame: rows = observations, columns = time points
#' @export
predict_orsf <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "aorsf")
  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  pred <- predict(
    object$model,
    new_data = newdata,
    pred_horizon = times,
    pred_type = "surv",
    na_action = "pass",
    boundary_checks = FALSE,
    ...
  )

  if (is.vector(pred)) {
    pred <- matrix(pred, nrow = 1)
  }

  colnames(pred) <- paste0("t=", times)
  rownames(pred) <- paste0("ID_", seq_len(nrow(newdata)))

  as.data.frame(pred)
}

library(aorsf)
library(survival)

mod_orsf <- fit_orsf(Surv(time, status) ~ age + ph.ecog, data = lung)
pred_orsf <- predict_orsf(mod_orsf, newdata = lung[1:5, ], times = c(100, 200, 300))
print(round(pred_orsf, 3))

cv_results_orsf <- cv_survlearner(
  formula = Surv(time, status) ~ age + ph.ecog,
  data = lung,
  fit_fun = fit_orsf,
  pred_fun = predict_orsf,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42,
  verbose = FALSE
)

cv_summary(cv_results_orsf)
cv_plot(cv_results_orsf)


#----------- add tuner

tune_orsf <- function(formula, data, times,
                      param_grid = expand.grid(
                        n_tree = c(100, 300),
                        mtry = c(2, 3),
                        min_events = c(5, 10)
                      ),
                      metrics = c("cindex", "ibs"),
                      folds = 5,
                      seed = 123,
                      ...) {

  purrr::pmap_dfr(param_grid, function(n_tree, mtry, min_events) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_orsf,
      pred_fun = predict_orsf,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      n_tree = n_tree,
      mtry = mtry,
      n_split = min_events,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      n_tree = n_tree,
      mtry = mtry,
      min_events = min_events
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  }) |> dplyr::arrange(dplyr::across(any_of(metrics[1])))
}



res_orsf <- tune_orsf(
  formula = Surv(time, status) ~ age + ph.ecog,
  data = lung,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    n_tree = c(100, 300),
    mtry = c(1, 2),             # o larger than number of covariates
    min_events = c(5, 10)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3
)

print(res_orsf)
