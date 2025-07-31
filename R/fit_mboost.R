
fit_mboost <- function(formula, data, weights = NULL, mstop = 100, nu = 0.1,
                       minsplit = 10, minbucket = 4, maxdepth = 2, ...) {
  stopifnot(requireNamespace("mboost", quietly = TRUE))
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  data <- as.data.frame(data)
  family <- mboost::CoxPH()

  control <- mboost::boost_control(mstop = mstop, nu = nu, ...)
  tree_control <- partykit::ctree_control(minsplit = minsplit, minbucket = minbucket, maxdepth = maxdepth)

  # conditionally include weights
  model <- if (!is.null(weights)) {
    mboost::blackboost(
      formula = formula,
      data = data,
      weights = weights,
      family = family,
      control = control,
      tree_controls = tree_control
    )
  } else {
    mboost::blackboost(
      formula = formula,
      data = data,
      family = family,
      control = control,
      tree_controls = tree_control
    )
  }

  structure(
    list(
      model = model,
      learner = "mboost",
      formula = formula,
      data = data,
      time = all.vars(formula)[[2]],
      status = all.vars(formula)[[3]]
    ),
    class = "mlsurv_model",
    engine = "mboost"
  )
}



survival_curve_to_prob <- function(eval_time, event_times, survival_prob) {

    if (event_times[1] != -Inf) {
    event_times <- c(-Inf, event_times)
    survival_prob <- rbind(1, survival_prob)
  }
  if (event_times[length(event_times)] != Inf) {
    event_times <- c(event_times, Inf)
    survival_prob <- rbind(survival_prob, 0)
  }

  index <- findInterval(eval_time, event_times)
  survival_prob[index, , drop = FALSE]
}


predict_mboost <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "mboost")
  stopifnot(requireNamespace("mboost", quietly = TRUE))

  surv_fit <- mboost::survFit(object$model, newdata = newdata)

  sf_time <- surv_fit$time
  sf_surv <- surv_fit$surv  # rows = time, cols = obs

  # if max(eval times) > max(sf_time), pad with last survival value
  if (max(times) > max(sf_time)) {
    last_time <- max(sf_time)
    extra_times <- sort(setdiff(times, sf_time))
    extra_times <- extra_times[extra_times > last_time]

    if (length(extra_times) > 0) {
      pad_vals <- matrix(
        rep(tail(sf_surv, 1), each = length(extra_times)),
        nrow = length(extra_times), byrow = TRUE
      )
      sf_time <- c(sf_time, extra_times)
      sf_surv <- rbind(sf_surv, pad_vals)
    }
  }

  # add boundaries to support interpolation
  sf_time_ext <- c(-Inf, sf_time)
  sf_surv_ext <- rbind(1, sf_surv)

  # interpolate at requested times (constant survival interpolation)
  idx <- findInterval(times, sf_time_ext)
  surv_probs <- sf_surv_ext[idx, , drop = FALSE]
  surv_probs <- t(surv_probs)

  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- rownames(newdata)

  return(as.data.frame(surv_probs))
}




library(survival)
data("veteran", package = "survival")

mod_mboost <- fit_mboost(Surv(time, status) ~ age + karno + celltype, data = veteran)

predict_mboost(mod_mboost, newdata = veteran[1:5, ], times = c(5, 10, 40))




# Load data and required packages
library(survival)
data("veteran", package = "survival")

# Fit model with CV
cv_results_mboost <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  fit_fun = fit_mboost,
  pred_fun = predict_mboost,
  times = c(5, 10, 40),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42
)

# Summary and visualization
print(cv_results_mboost)

#----------- add tuner


tune_mboost <- function(formula, data, times,
                        param_grid = expand.grid(
                          mstop = c(50, 100, 200),
                          nu = c(0.05, 0.1),
                          maxdepth = c(2, 3)
                        ),
                        metrics = c("cindex", "ibs"),
                        folds = 5,
                        seed = 123,
                        refit_best = FALSE,
                        ...) {

  results <- purrr::pmap_dfr(param_grid, function(mstop, nu, maxdepth) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_mboost,
      pred_fun = predict_mboost,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      mstop = mstop,
      nu = nu,
      maxdepth = maxdepth,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      mstop = mstop,
      nu = nu,
      maxdepth = maxdepth
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
    best_model <- fit_mboost(
      formula = formula,
      data = data,
      mstop = best_row$mstop,
      nu = best_row$nu,
      maxdepth = best_row$maxdepth,
      ...
    )
    return(best_model)
  }
}

# Grid search only
res_mboost <- tune_mboost(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(5, 10, 40),
  param_grid = expand.grid(
    mstop = c(100, 200),
    nu = c(0.05, 0.1),
    maxdepth = c(2, 3)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3
)
print(res_mboost)

# Refit best model
mod_mboost_best <- tune_mboost(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(5, 10, 40),
  param_grid = expand.grid(
    mstop = c(100, 200),
    nu = c(0.05, 0.1),
    maxdepth = c(2, 3)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)
summary(mod_mboost_best)

