#' Fit a Componentwise Gradient Boosted Cox Model (blackboost)
#'
#' Fits a survival boosting model using \pkg{mboost}'s \code{blackboost} with
#' a Cox proportional hazards loss (\code{mboost::CoxPH()}) and shallow tree
#' base-learners (\pkg{partykit}). Returns an `mlsurv_model` compatible with the
#' `survalis` pipeline.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A `data.frame` containing all variables referenced in \code{formula}.
#' @param weights Optional case weights passed to \code{mboost::blackboost()}.
#' @param mstop Integer; number of boosting iterations (stopping iteration). Default \code{100}.
#' @param nu Learning rate/shrinkage in \code{(0,1]}; default \code{0.1}.
#' @param minbucket Minimum bucket size per tree node.
#' @param maxdepth  Maximum tree depth.
#' @param minsplit, minbucket, maxdepth Tree control parameters passed via
#'   \code{partykit::ctree_control()} (defaults: \code{10}, \code{4}, \code{2}).
#' @param ... Additional arguments forwarded to \code{mboost::boost_control()} and/or
#'   \code{mboost::blackboost()}.
#'
#' @details
#' The base-learner is a conditional inference tree controlled by
#' \code{partykit::ctree_control()}. The loss is the partial likelihood for Cox PH
#' via \code{mboost::CoxPH()}. Use \code{mstop} and \code{nu} to control complexity.
#'
#' @return An object of class `mlsurv_model` with elements:
#' \describe{
#'   \item{model}{Fitted \code{mboost} model.}
#'   \item{learner}{\code{"blackboost"}.}
#'   \item{formula, data, time, status}{Original inputs/metadata.}
#' }
#'
#' @references
#' Bühlmann P, Hothorn T (2007). Boosting algorithms: regularization, prediction
#' and model fitting. \emph{Statistical Science}.
#'
#' @examples
#' \donttest{
#'   mod <- fit_blackboost(
#'     Surv(time, status) ~ age + karno + celltype,
#'     data = veteran
#'   )
#'   head(predict_blackboost(mod, newdata = veteran[1:5, ], times = c(100, 200)))
#' }
#' @keywords internal
#' @noRd

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

  model <- suppressMessages({
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
  })

  structure(
    list(
      model = model,
      learner = "blackboost",
      formula = formula,
      data = data,
      time = all.vars(formula[[2]])[[1]],
      status = all.vars(formula[[2]])[[2]]
    ),
    class = "mlsurv_model",
    engine = "mboost"
  )
}

#' Predict Survival Probabilities from a blackboost Model
#'
#' Generates survival probabilities at requested times from a fitted
#' \pkg{mboost} Cox boosting model produced by `fit_blackboost()`.
#'
#' @param object An `mlsurv_model` returned by `fit_blackboost()`.
#' @param newdata A `data.frame` of new observations for prediction.
#' @param times Numeric vector of time points at which to compute survival
#'   probabilities.
#' @param ... Additional arguments forwarded to \code{mboost::survFit()}.
#'
#' @details
#' Predictions use \code{mboost::survFit()} to obtain a step function
#' \(S(t)\) per observation, evaluated at \code{times} as a step function
#' (survival 1 before the first event time, last value carried forward beyond
#' the last one).
#'
#' @return A base `data.frame` with one row per observation in \code{newdata} and
#' columns named \code{"t=<time>"} (character), containing survival probabilities
#' in \[0, 1].
#'
#' @examples
#' mod <- fit_blackboost(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' predict_blackboost(mod, newdata = veteran[1:5, ], times = c(5, 10, 40))
#' @keywords internal
#' @noRd

predict_blackboost <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "blackboost") {
    warning("Object passed to predict_blackboost() may not come from fit_blackboost().")
  }

  stopifnot(requireNamespace("mboost", quietly = TRUE))

  newdata <- as.data.frame(newdata)
  n_new <- nrow(newdata)
  # mboost::survFit() fails on a single row: duplicate it and keep the first result.
  nd <- if (n_new == 1L) newdata[c(1L, 1L), , drop = FALSE] else newdata
  surv_fit <- mboost::survFit(object$model, newdata = nd)

  # surv_fit$surv: rows = time, cols = obs. Step lookup with S = 1 before the first
  # grid time and right-constant beyond the last one.
  S <- t(as.matrix(surv_fit$surv))
  survmat <- .step_surv(S, surv_fit$time, times)[seq_len(n_new), , drop = FALSE]
  .finalize_survmat(survmat, times = times)
}

#' Tune blackboost Hyperparameters (Cross-Validation)
#'
#' Cross-validates \pkg{mboost} Cox boosting models over a hyperparameter grid
#' and selects the best configuration according to the primary metric. Optionally
#' refits the best model on the full dataset.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A `data.frame` containing variables referenced in \code{formula}.
#' @param times Numeric vector of evaluation time points (same scale as the training time).
#' @param param_grid A `data.frame` (e.g., from \code{expand.grid()}) with columns
#'   such as \code{mstop}, \code{nu}, and \code{maxdepth}.
#' @param metrics Character vector of metrics to compute (e.g., \code{"cindex"}, \code{"ibs"}).
#'   The first entry is the primary selection metric.
#' @param folds Integer; number of CV folds. Default \code{5}.
#' @param seed Integer random seed for reproducibility. Default \code{123}.
#' @param ncores Integer number of CPU cores passed to \code{cv_survlearner}
#'   for fold evaluation (default \code{1}).
#' @param refit_best Logical; if \code{TRUE}, refits the best configuration on the full
#'   data and returns it. Default \code{FALSE}.
#' @param ... Additional arguments forwarded to `fit_blackboost()`.
#'
#' @return If \code{refit_best = FALSE}, a `data.frame` (class `"tuned_surv"`) of grid
#' results with metric columns, sorted by the primary metric. If \code{refit_best = TRUE},
#' a fitted `mlsurv_model` from `fit_blackboost()` using the selected hyperparameters.
#'
#' @details
#' Internally calls \code{cv_survlearner()} with \code{fit_blackboost()}/\code{predict_blackboost()}
#' so tuning mirrors the production path. Typical grids vary \code{mstop}, \code{nu}, and
#' tree \code{maxdepth}.
#'
#' @examples
#' \donttest{
#' res_blackboost <- tune_blackboost(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(5, 10, 40),
#'   param_grid = expand.grid(
#'     mstop = c(100, 200),
#'     nu = c(0.05, 0.1),
#'     maxdepth = c(2, 3)
#'   ),
#'   metrics = c("cindex", "ibs"),
#'   folds = 3
#' )
#' print(res_blackboost)
#'
#' mod_blackboost_best <- tune_blackboost(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(5, 10, 40),
#'   param_grid = expand.grid(
#'     mstop = c(100, 200),
#'     nu = c(0.05, 0.1),
#'     maxdepth = c(2, 3)
#'   ),
#'   metrics = c("cindex", "ibs"),
#'   folds = 3,
#'   refit_best = TRUE
#' )
#' summary(mod_blackboost_best)
#' }
#' @keywords internal
#' @noRd

tune_blackboost <- function(formula, data, times,
                        param_grid = expand.grid(
                          mstop = c(50, 100, 200),
                          nu = c(0.05, 0.1),
                          maxdepth = c(2, 3)
                        ),
                        metrics = c("cindex", "ibs"),
                        folds = 5,
                        seed = 123,
                        ncores = 1,
                        refit_best = FALSE,
                        ...) {

  results <- .pmap_rbind_dt(param_grid, function(mstop, nu, maxdepth) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_blackboost,
      pred_fun = predict_blackboost,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ncores = ncores,
      mstop = mstop,
      nu = nu,
      maxdepth = maxdepth,
      ...
    )

    summary <- cv_summary(cv_results)
    .wide_metric_row(list(mstop = mstop, nu = nu, maxdepth = maxdepth), summary)
  }, learner = "blackboost")

  maximize <- .metric_maximize(metrics[1])
  results <- .arrange_by_metric_dt(results, metrics[1], maximize)

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
