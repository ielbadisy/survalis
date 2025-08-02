
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

tune_bnnsurv <- function(formula, data, times,
                         param_grid = expand.grid(
                           k = c(2, 3),
                           num_base_learners = c(30, 50),
                           sample_fraction = c(0.5, 1),
                           stringsAsFactors = FALSE
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5,
                         seed = 123,
                         refit_best = FALSE,
                         ...) {
  results <- purrr::pmap_dfr(param_grid, function(k, num_base_learners, sample_fraction) {
    cv_result <- tryCatch({
      cv_survlearner(
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
    }, error = function(e) NULL)

    if (is.null(cv_result)) {
      return(tibble::tibble(
        k = k,
        num_base_learners = num_base_learners,
        sample_fraction = sample_fraction,
        failed = TRUE
      ))
    }

    summary <- cv_summary(cv_result)

    tibble::tibble(
      k = k,
      num_base_learners = num_base_learners,
      sample_fraction = sample_fraction,
      failed = FALSE
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  results <- results |> dplyr::filter(!failed)

  if (nrow(results) == 0) {
    warning("All tuning combinations failed.")
    return(tibble::tibble())
  }

  results <- dplyr::arrange(results, dplyr::desc(.data[[metrics[1]]]))

  if (refit_best) {
    best <- results[1, ]
    model <- fit_bnnsurv(
      formula = formula,
      data = data,
      k = best$k,
      num_base_learners = best$num_base_learners,
      sample_fraction = best$sample_fraction,
      ...
    )
    return(model)
  }

  return(results)
}


param_grid = expand.grid(
  k = c(2),
  num_base_learners = c(30),
  sample_fraction = c(0.5, 1)
)

# Return performance table
res_bnnsurv <- tune_bnnsurv(
  formula = Surv(time, status) ~ age + karno + diagtime + prior,
  data = veteran,
  times = c(100, 200),
  param_grid = param_grid,
  refit_best = FALSE
)
print(res_bnnsurv)

# Return best model (mlsurv_model class)
mod_bnnsurv <- tune_bnnsurv(
  formula = Surv(time, status) ~ age + karno + diagtime + prior,
  data = veteran,
  times = c(100, 200),
  param_grid = param_grid,
  refit_best = TRUE
)

summary(mod_bnnsurv)
predict_bnnsurv(mod_bnnsurv, newdata = veteran[1:5, ], times = c(100, 200))
