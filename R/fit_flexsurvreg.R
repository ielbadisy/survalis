#' Fit a Parametric Survival Regression Model Using flexsurvreg
#'
#' This function fits a fully parametric survival regression model using the
#' \pkg{flexsurv} package. It supports a variety of parametric distributions
#' (e.g., Weibull, exponential, log-normal) and returns a model object compatible
#' with the `mlsurv_model` interface for downstream survival predictions.
#'
#' @param formula A \code{\link[survival]{Surv}}-based survival formula of the form
#'   \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing the variables in the model.
#' @param dist Character string specifying the parametric distribution to use
#'   (default = `"weibull"`). See \code{\link[flexsurv]{flexsurvreg}} for available
#'   options.
#' @param ... Additional arguments passed to \code{\link[flexsurv]{flexsurvreg}}.
#'
#' @return An object of class \code{mlsurv_model} containing:
#'   \itemize{
#'     \item \code{model} - the fitted \code{flexsurvreg} model
#'     \item \code{learner} - string identifier `"flexsurvreg"`
#'     \item \code{formula} - the model formula
#'     \item \code{data} - the training data
#'     \item \code{time} - the name of the time-to-event column
#'     \item \code{status} - the name of the event indicator column
#'   }
#'
#' @examples
#' 
#' mod_flex <- fit_flexsurvreg(Surv(time, status) ~ age + celltype + karno,
#'                             data = veteran,
#'                             dist = "weibull")
#' @export

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

#' Predict Survival Probabilities from a flexsurvreg Model
#'
#' This function generates survival probability predictions at specified time
#' points from a model fitted with \code{\link{fit_flexsurvreg}}.
#'
#' @param object A fitted \code{mlsurv_model} object returned by
#'   \code{\link{fit_flexsurvreg}}.
#' @param newdata A data frame of new observations for which to predict survival
#'   probabilities.
#' @param times Numeric vector of time points at which to estimate survival
#'   probabilities.
#' @param ... Additional arguments passed to \code{\link[flexsurv]{predict.flexsurvreg}}.
#'
#' @return A data frame of survival probabilities with one row per observation in
#'   \code{newdata} and one column per time point. Column names are of the form
#'   `"t=100"`, `"t=200"`, etc.
#'
#' @examples
#' mod_flex <- fit_flexsurvreg(Surv(time, status) ~ age + celltype + karno,
#'                             data = veteran,
#'                             dist = "weibull")
#' predict_flexsurvreg(mod_flex, newdata = veteran[1:5, ],
#'                     times = c(100, 200, 300))
#' @export

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


#' Tune Parametric Survival Models with flexsurvreg
#'
#' Performs hyperparameter tuning for \code{\link{fit_flexsurvreg}} over a grid
#' of parametric distributions using cross-validation. Returns either a summary
#' table of performance metrics or a refitted best model.
#'
#' @param formula A \code{\link[survival]{Surv}}-based survival formula.
#' @param data A data frame containing the variables in the model.
#' @param times Numeric vector of time points at which to evaluate performance.
#' @param param_grid Character vector of distributions to test
#'   (default = \code{c("weibull", "exponential", "lognormal")}).
#' @param metrics Character vector of performance metrics to compute
#'   (default = \code{c("cindex", "ibs")}).
#' @param folds Integer, number of cross-validation folds (default = 5).
#' @param seed Random seed for reproducibility.
#' @param refit_best Logical; if \code{TRUE}, refits the best-performing model
#'   on the full dataset and returns it. If \code{FALSE}, returns the tuning results.
#' @param ... Additional arguments passed to \code{\link{fit_flexsurvreg}}.
#'
#' @return If \code{refit_best = FALSE}, returns a data frame of tuning results
#'   with one row per distribution and columns for each metric.  
#'   If \code{refit_best = TRUE}, returns the best \code{mlsurv_model} object.
#'
#' @examples
#'
#' # Distribution tuning without refitting
#' res <- tune_flexsurvreg(Surv(time, status) ~ age + karno + celltype,
#'                         data = veteran,
#'                         param_grid = c("weibull", "exponential", "lognormal"),
#'                         times = c(100, 200, 300))
#' print(res)
#'
#' # Refit the best distribution
#' best_mod <- tune_flexsurvreg(Surv(time, status) ~ age + karno + celltype,
#'                              data = veteran,
#'                              param_grid = c("weibull", "exponential", "lognormal"),
#'                              times = c(100, 200, 300),
#'                              refit_best = TRUE)
#' summary(best_mod)
#' @export

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

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))
  } else {
    results <- results |> dplyr::arrange(!!rlang::sym(metrics[1]))
  }

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
