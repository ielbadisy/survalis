#' Fit a Survival Tree Model using `rpart`
#'
#' Fits a survival tree using the \pkg{rpart} package with an exponential splitting rule.
#' This learner is compatible with the `mlsurv_model` interface and can be used with
#' the `cv_survlearner()` and `tune_*()` functions in the \pkg{survalis} framework.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#' @param data A `data.frame` containing the variables in the model.
#' @param minsplit Minimum number of observations that must exist in a node in order for a split to be attempted.
#'   Default is `20`.
#' @param minbucket Minimum number of observations in any terminal node. Default is `round(minsplit / 3)`.
#' @param cp Complexity parameter for pruning. Default is `0.01`.
#' @param maxcompete Number of competitor splits retained in the output. Default is `4`.
#' @param maxsurrogate Number of surrogate splits retained in the output. Default is `5`.
#' @param usesurrogate How surrogates are used in the splitting process. Default is `2`.
#' @param xval Number of cross-validations to perform in `rpart`. Default is `10`.
#' @param surrogatestyle Controls selection of surrogate splits. Default is `0`.
#' @param maxdepth Maximum depth of any node of the final tree. Default is `30`.
#'
#' @return An object of class `"mlsurv_model"` containing:
#' \itemize{
#'   \item `model` -- fitted `rpart` survival tree object
#'   \item `learner` -- `"rpart"`
#'   \item `formula`, `data`, `time`, and `status`
#' }
#'
#' @details
#' This function fixes the `method = "exp"` argument to fit an exponential
#' splitting survival tree, which models the hazard function assuming
#' exponential survival within terminal nodes.
#'
#' @references
#' Therneau TM, Atkinson EJ. (2019). *An Introduction to Recursive Partitioning Using the RPART Routines*. Mayo Clinic.
#'
#' @examples
#' mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' pred_rpart <- predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300))
#' head(pred_rpart)
#'
#' @export

fit_rpart <- function(formula, data,
                      minsplit = 20,
                      minbucket = round(minsplit / 3),
                      cp = 0.01,
                      maxcompete = 4,
                      maxsurrogate = 5,
                      usesurrogate = 2,
                      xval = 10,
                      surrogatestyle = 0,
                      maxdepth = 30) {
  stopifnot(requireNamespace("rpart", quietly = TRUE))

  method <- "exp"  # fixed to survival tree with exponential splitting

  model <- rpart::rpart(
    formula,
    data = data,
    method = method,
    control = list(
      minsplit = minsplit, minbucket = minbucket, cp = cp,
      maxcompete = maxcompete, maxsurrogate = maxsurrogate,
      usesurrogate = usesurrogate, xval = xval,
      surrogatestyle = surrogatestyle, maxdepth = maxdepth
    ),
    model = TRUE, x = TRUE, y = TRUE
  )

  structure(list(
    model = model,
    learner = "rpart",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2]
  ), class = "mlsurv_model", engine = "rpart")
}

#' Predict Survival Probabilities from an `rpart` Survival Tree
#'
#' Predicts survival probabilities at specified time points from a fitted
#' `rpart` survival tree model, assuming an exponential survival distribution
#' within each terminal node.
#'
#' @param object An `"mlsurv_model"` object created by [fit_rpart()].
#' @param newdata A `data.frame` containing the predictor variables for prediction.
#' @param times Numeric vector of time points at which to compute survival probabilities.
#' @param ... Additional arguments passed to internal prediction functions.
#'
#' @return A `data.frame` with one row per observation in `newdata` and one column
#' per requested time point, containing predicted survival probabilities.
#'
#' @details
#' Predictions are based on converting the predicted mean survival times from the
#' survival tree into survival probabilities under an exponential assumption:
#' \deqn{S(t) = \exp(-t / \hat{\mu})}
#' where \eqn{\hat{\mu}} is the predicted mean survival time from the terminal node.
#'
#' @examples
#'
#' mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export

predict_rpart <- function(object, newdata, times, ...) {

  if (!is.null(object$learner) && object$learner != "rpart") {
    warning("Object passed to predict_rpart() may not come from fit_rpart().")
  }

  stopifnot(object$learner == "rpart")
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  model_party <- partykit::as.party(object$model)

  # predict expected survival times from the tree
  predicted_times <- stats::predict(model_party, newdata = newdata, type = "response")

  # convert predicted survival times to survival probabilities using exponential assumption
  survmat <- sapply(times, function(t) {
    exp(-t / predicted_times)
  })

  .finalize_survmat(survmat, times = times)
}

#' Tune a Survival Tree Model (`rpart`) via Cross-Validation
#'
#' Performs hyperparameter tuning for the `rpart` survival tree model using
#' cross-validation, returning either a table of results or the best fitted model.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#' @param data A `data.frame` containing the variables in the model.
#' @param times Numeric vector of time points for evaluation.
#' @param param_grid A `data.frame` or `expand.grid()` result specifying hyperparameter combinations.
#'   Must include `minsplit`, `cp`, and `maxdepth`.
#' @param metrics Character vector of evaluation metrics to compute (e.g., `"cindex"`, `"ibs"`).
#' @param folds Number of cross-validation folds. Default is `5`.
#' @param seed Random seed for reproducibility.
#' @param refit_best Logical; if `TRUE`, the best hyperparameter combination is refitted on the full dataset.
#' @param ... Additional arguments passed to [fit_rpart()].
#'
#' @return If `refit_best = FALSE`, returns a tibble summarizing mean CV performance for each parameter combination.
#' If `refit_best = TRUE`, returns an `"mlsurv_model"` object for the best parameters, with a `"tuning_results"` attribute.
#'
#' @examples
#'
#' # Example tuning without refit
#' res_rpart <- tune_rpart(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   param_grid = expand.grid(
#'     minsplit = c(10, 20),
#'     cp = c(0.001, 0.01),
#'     maxdepth = c(10, 30)
#'   ),
#'   metrics = c("cindex", "ibs"),
#'   folds = 3,
#'   seed = 42,
#'   refit_best = FALSE
#' )
#'
#' print(res_rpart)
#'
#' # Example tuning with refit
#' mod_rpart_best <- tune_rpart(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   param_grid = expand.grid(
#'     minsplit = c(10, 20),
#'     cp = c(0.001, 0.01),
#'     maxdepth = c(10, 30)
#'   ),
#'   metrics = c("cindex", "ibs"),
#'   folds = 3,
#'   seed = 42,
#'   refit_best = TRUE
#' )
#'
#' @export

tune_rpart <- function(formula, data, times,
                       param_grid = expand.grid(
                         minsplit = c(10, 20),
                         cp = c(0.001, 0.01),
                         maxdepth = c(10, 30)
                       ),
                       metrics = c("cindex", "ibs"),
                       folds = 5,
                       seed = 123,
                       refit_best = TRUE,
                       ...) {
  res <- purrr::pmap_dfr(param_grid, function(minsplit, cp, maxdepth) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_rpart,
      pred_fun = predict_rpart,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      minsplit = minsplit,
      cp = cp,
      maxdepth = maxdepth,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      minsplit = minsplit,
      cp = cp,
      maxdepth = maxdepth
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    res <- dplyr::arrange(res, dplyr::desc(!!sym(metrics[1])))
  } else {
    res <- dplyr::arrange(res, !!sym(metrics[1]))
  }

  if (refit_best) {
    best_params <- res[1, ]
    model <- fit_rpart(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      minsplit = best_params$minsplit,
      cp = best_params$cp,
      maxdepth = best_params$maxdepth,
      ...
    )
    attr(model, "tuning_results") <- res
    class(model) <- c("tune_surv", class(model))
    return(model)
  }

  return(res)
}
