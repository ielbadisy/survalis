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
#' splitting survival tree (exponential deviance splits). Predictions are the
#' Kaplan-Meier curves of the training subjects in each terminal node.
#'
#' @references
#' Therneau TM, Atkinson EJ. (2019). *An Introduction to Recursive Partitioning Using the RPART Routines*. Mayo Clinic.
#'
#' @examples
#' mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' pred_rpart <- predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300))
#' head(pred_rpart)
#'
#' @keywords internal
#' @noRd

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
#' `rpart` survival tree model, using the Kaplan-Meier curve of the training
#' observations in each terminal node.
#'
#' @param object An `"mlsurv_model"` object created by `fit_rpart()`.
#' @param newdata A `data.frame` containing the predictor variables for prediction.
#' @param times Numeric vector of time points at which to compute survival probabilities.
#' @param ... Additional arguments passed to internal prediction functions.
#'
#' @return A `data.frame` with one row per observation in `newdata` and one column
#' per requested time point, containing predicted survival probabilities.
#'
#' @details
#' Each new observation is dropped down the tree to its terminal node. The
#' predicted survival curve is the Kaplan-Meier estimate computed from the
#' training observations in that node, read as a step function with
#' \eqn{S(t) = 1} before the first event time in the node.
#'
#' @examples
#'
#' mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @keywords internal
#' @noRd

predict_rpart <- function(object, newdata, times, ...) {

  if (!is.null(object$learner) && object$learner != "rpart") {
    warning("Object passed to predict_rpart() may not come from fit_rpart().")
  }

  stopifnot(object$learner == "rpart")
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  model_party <- partykit::as.party(object$model)
  leaf_km <- .rpart_leaf_km(model_party)

  # map each new observation to its terminal node
  newdata <- as.data.frame(newdata)
  nodes <- as.character(stats::predict(model_party, newdata = newdata, type = "node"))

  # Kaplan-Meier curve of the training observations in each terminal node,
  # read as a right-continuous step function with S(t) = 1 before the first
  # event time
  survmat <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(times))
  for (nd in unique(nodes)) {
    km <- leaf_km[[nd]]
    idx <- which(nodes == nd)
    row <- .step_surv(matrix(km$surv, nrow = 1), km$time, times)
    survmat[idx, ] <- matrix(row, nrow = length(idx), ncol = length(times),
                             byrow = TRUE)
  }

  .finalize_survmat(survmat, times = times)
}

# Per-terminal-node Kaplan-Meier curves from the training data of a party tree
.rpart_leaf_km <- function(model_party) {
  fitted_nodes <- model_party$fitted[["(fitted)"]]
  y <- model_party$fitted[["(response)"]]
  y <- as.matrix(y)
  split_idx <- split(seq_along(fitted_nodes), fitted_nodes)
  lapply(split_idx, function(i) {
    km <- survival::survfit(survival::Surv(y[i, 1], y[i, 2]) ~ 1)
    list(time = km$time, surv = km$surv)
  })
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
#' @param ncores Integer number of CPU cores passed to \code{cv_survlearner}
#'   for fold evaluation (default \code{1}).
#' @param refit_best Logical; if `TRUE`, the best hyperparameter combination is refitted on the full dataset.
#' @param ... Additional arguments passed to `fit_rpart()`.
#'
#' @return If `refit_best = FALSE`, returns a tibble summarizing mean CV performance for each parameter combination.
#' If `refit_best = TRUE`, returns an `"mlsurv_model"` object for the best parameters, with a `"tuning_results"` attribute.
#'
#' @examples
#'
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
#' @keywords internal
#' @noRd

tune_rpart <- function(formula, data, times,
                       param_grid = expand.grid(
                         minsplit = c(10, 20),
                         cp = c(0.001, 0.01),
                         maxdepth = c(10, 30)
                       ),
                       metrics = c("cindex", "ibs"),
                       folds = 5,
                       seed = 123,
                       ncores = 1,
                       refit_best = TRUE,
                       ...) {
  res <- .pmap_rbind_dt(param_grid, function(minsplit, cp, maxdepth) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_rpart,
      pred_fun = predict_rpart,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ncores = ncores,
      minsplit = minsplit,
      cp = cp,
      maxdepth = maxdepth,
      ...
    )

    summary <- cv_summary(cv_results)
    .wide_metric_row(list(minsplit = minsplit, cp = cp, maxdepth = maxdepth), summary)
  }, learner = "rpart")

  maximize <- .metric_maximize(metrics[1])
  res <- .arrange_by_metric_dt(res, metrics[1], maximize)

  if (refit_best) {
    best_params <- res[1, ]
    model <- fit_rpart(
      formula = formula,
      data = .complete_cases_df(data, all.vars(formula)),
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
