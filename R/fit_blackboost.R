fit_blackboost <- function(formula, data, weights = NULL, mstop = 100, nu = 0.1,
                       minsplit = 10, minbucket = 4, maxdepth = 2, ...) {

  stopifnot(requireNamespace("mboost", quietly = TRUE))
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  data <- as.data.frame(data)
  family <- mboost::CoxPH()

  control <- mboost::boost_control(mstop = mstop, nu = nu, ...)
  tree_control <- partykit::ctree_control(
    minsplit = minsplit,
    minbucket = minbucket,
    maxdepth = maxdepth
  )

  model <- suppressMessages(suppressWarnings({
    if (!is.null(weights)) {
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
  }))

  structure(
    list(
      model = model,
      learner = "blackboost",
      formula = formula,
      data = data,
      time = all.vars(formula)[[2]],
      status = all.vars(formula)[[3]]
    ),
    class = "mlsurv_model",
    engine = "mboost"
  )
}

predict_blackboost <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "blackboost") {
    warning("Object passed to predict_blackboost() may not come from fit_blackboost().")
  }

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
  survmat <- sf_surv_ext[idx, , drop = FALSE]
  survmat <- t(survmat)

  colnames(survmat) <- paste0("t=", times)
  rownames(survmat) <- rownames(newdata)

  as.data.frame(survmat)
}




library(survival)
data("veteran", package = "survival")

mod <- fit_blackboost(Surv(time, status) ~ age + karno + celltype, data = veteran)

predict_blackboost(mod, newdata = veteran[1:5, ], times = c(5, 10, 40))




# Load data and required packages
library(survival)
data("veteran", package = "survival")

# Fit model with CV
cv_results_blackboost <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  fit_fun = fit_blackboost,
  pred_fun = predict_blackboost,
  times = c(5, 10, 40),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42
)

# Summary and visualization
print(cv_results_blackboost)

#----------- add tuner


tune_blackboost <- function(formula, data, times,
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
      fit_fun = fit_blackboost,
      pred_fun = predict_blackboost,
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
    best_model <- fit_blackboost(
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
res_blackboost <- tune_blackboost(
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
print(res_blackboost)

# Refit best model
mod_blackboost_best <- tune_blackboost(
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
summary(mod_blackboost_best)

