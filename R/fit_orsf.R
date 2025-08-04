
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
    learner = "orsf", # Oblique Random Survival Forests

    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "aorsf")
}


predict_orsf <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "orsf") {
    warning("Object passed to predict_orsf() may not come from fit_orsf().")
  }

  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  survmat <- predict(
    object$model,
    new_data = newdata,
    pred_horizon = times,
    pred_type = "surv",
    na_action = "pass",
    boundary_checks = FALSE,
    ...
  )

  if (is.vector(survmat)) {
    survmat <- matrix(survmat, nrow = 1)
  }

  colnames(survmat) <- paste0("t=", times)

  as.data.frame(survmat)
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
