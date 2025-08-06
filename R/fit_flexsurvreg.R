fit_flexsurvreg <- function(formula, data, dist = "weibull", ...) {

  stopifnot(requireNamespace("flexsurv", quietly = TRUE))

  model <- flexsurv::flexsurvreg(formula = formula, data = data, dist = dist, ...)

  structure(list(
    model = model,
    learner = "flexsurvreg",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "flexsurv")
}

predict_flexsurvreg <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "flexsurvreg") {
    warning("Object passed to predict_flexsurvreg() may not come from fit_flexsurvreg().")
  }


  stopifnot(requireNamespace("flexsurv", quietly = TRUE))


  pred <- predict(
    object$model,
    newdata = newdata,
    type = "survival",
    times = times
  )

  # Unpack .pred (a list-column of tibbles)
  .pred <- pred$.pred
  n <- nrow(newdata)
  k <- length(times)
  survmat <- matrix(NA_real_, nrow = n, ncol = k)
  colnames(survmat) <- paste0("t=", times)

  for (i in seq_len(n)) {
    survmat[i, ] <- .pred[[i]]$.pred_survival
    }

  as.data.frame(survmat)
}

library(survival)
library(flexsurv)

data(veteran, package = "survival")

# Fit the model
mod_flex <- fit_flexsurvreg(Surv(time, status) ~ age + celltype + karno,
                         data = veteran,
                         dist = "weibull")

# Predict survival probabilities at given time points
times <- c(100, 200, 300)
predicted_surv <- predict_flexsurvreg(mod_flex, newdata = veteran[1:5, ], times = times)

predicted_surv



#------------ add tuner
tune_flexsurvreg <- function(formula, data, times,
                          param_grid = c("weibull", "exponential", "lognormal"),
                          metrics = c("cindex", "ibs"),
                          folds = 5, seed = 123,
                          refit_best = FALSE, ...) {

  results <- purrr::map_dfr(param_grid, function(d) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_flexsurvreg,
      pred_fun = predict_flexsurvreg,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      dist = d,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(dist = d) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(
          summary[, c("metric", "mean")],
          names_from = metric,
          values_from = mean
        )
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
    best_model <- fit_flexsurvreg(
      formula = formula,
      data = data,
      dist = best_row$dist,
      ...
    )
    return(best_model)  # will be of class mlsurv_model
  }
}


library(survival)

veteran <- survival::veteran


param_grid = c("weibull", "exponential", "lognormal")


# Get best dist only
res <- tune_flexsurvreg(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  param_grid = param_grid,
  times = c(100, 200, 300),
  refit_best = FALSE
  )

print(res)

# Fit and get best model (ready for predict(), summary())
mod <- tune_flexsurvreg(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  param_grid = param_grid,
  times = c(100, 200, 300),
  refit_best = TRUE
)

summary(mod)
predict_flexsurvreg(mod, newdata = veteran[1:5, ], times = 1:100)




cv_results_xgb <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  fit_fun = fit_flexsurvreg,
  pred_fun = predict_flexsurvreg,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs", "iae", "ise"),
  folds = 5,
  verbose = TRUE
)



