
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


#' Tune ORSF (Oblique Random Survival Forest) model using cross-validation
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data.frame with survival outcome and predictors
#' @param times Time points for survival prediction
#' @param param_grid A data.frame with tuning grid: n_tree, mtry, min_events
#' @param metrics Performance metrics to evaluate (e.g., "cindex", "ibs")
#' @param folds Number of folds for cross-validation
#' @param seed Random seed for reproducibility
#' @param refit_best Logical. If TRUE, returns best fitted model. Else, returns CV scores.
#' @param ... Additional arguments passed to `fit_orsf()`
#'
#' @return A tibble of scores (if `refit_best = FALSE`) or best model (if `refit_best = TRUE`)
#' @export
tune_orsf <- function(formula, data, times,
                      param_grid = expand.grid(
                        n_tree = c(100, 300),
                        mtry = c(2, 3),
                        min_events = c(5, 10)
                      ),
                      metrics = c("cindex", "ibs"),
                      folds = 5,
                      seed = 123,
                      refit_best = FALSE,
                      ...) {

  results <- purrr::pmap_dfr(param_grid, function(n_tree, mtry, min_events) {
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
  })

  results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_orsf(
      formula = formula,
      data = data,
      n_tree = best_row$n_tree,
      mtry = best_row$mtry,
      n_split = best_row$min_events,
      ...
    )
    return(best_model)
  }
}



# Grid search only
res_orsf <- tune_orsf(
  formula = Surv(time, status) ~ age + ph.ecog,
  data = lung,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    n_tree = c(100, 300),
    mtry = c(1, 2),
    min_events = c(5, 10)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3
)

print(res_orsf)

# Refit best
mod_orsf_best <- tune_orsf(
  formula = Surv(time, status) ~ age + ph.ecog,
  data = lung,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    n_tree = c(100, 300),
    mtry = c(1, 2),
    min_events = c(5, 10)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

summary(mod_orsf_best)
