#' Fit a fastgbm Survival Model
#'
#' Fits a histogram-based gradient boosting survival model with the
#' \pkg{fastgbm} package and returns an object compatible with the
#' `mlsurv_model` interface. The piecewise-exponential objective
#' (`objective = "pexp"`) is used by default because it models the hazard
#' jointly over covariates and time, so predicted survival curves can be
#' evaluated at arbitrary time points.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing the formula variables.
#' @param objective fastgbm survival objective: \code{"pexp"} (default),
#'   \code{"cox"} (Breslow baseline), or \code{"aft"} (normal error). All three give
#'   survival curves at arbitrary times.
#' @param ntrees Number of boosting rounds (default \code{200}).
#' @param learning_rate Shrinkage (default \code{0.1}).
#' @param max_depth Maximum tree depth (default \code{5}).
#' @param pexp_bins Number of piecewise-exponential time intervals (default
#'   \code{10}).
#' @param ... Additional arguments passed to \code{fastgbm::fastgbm()} (for
#'   example \code{min_node_size}, \code{subsample}, \code{colsample},
#'   \code{lambda}, \code{early_stopping} + \code{validation}, \code{seed}).
#'
#' @return An object of class \code{"mlsurv_model"} with elements \code{model},
#'   \code{learner} (\code{"fastgbm"}), \code{formula}, \code{data}, \code{time},
#'   \code{status}, and \code{objective}.
#'
#' @seealso \code{predict_fastgbm}, \code{fastgbm},
#'   \code{fit} for the verb interface
#'
#' @examplesIf requireNamespace("fastgbm", quietly = TRUE)
#' mod <- fit_fastgbm(Surv(time, status) ~ age + karno + celltype, veteran,
#'                    ntrees = 50, pexp_bins = 8)
#' summary(mod)
#'
#' @keywords internal
#' @noRd
fit_fastgbm <- function(formula, data, objective = c("pexp", "cox", "aft"),
                        ntrees = 200, learning_rate = 0.1, max_depth = 5,
                        pexp_bins = 10, ...) {
  stopifnot(requireNamespace("fastgbm", quietly = TRUE))
  objective <- match.arg(objective)

  model <- fastgbm::fastgbm(
    formula, data = data, objective = objective,
    ntrees = ntrees, learning_rate = learning_rate, max_depth = max_depth,
    pexp_bins = pexp_bins, verbose = FALSE, ...
  )

  time_status <- all.vars(formula[[2]])
  structure(list(
    model = model,
    learner = "fastgbm",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2],
    objective = objective
  ), class = "mlsurv_model", engine = "fastgbm")
}

#' Predict Survival Probabilities from a fastgbm Model
#'
#' Generates predicted survival probabilities at the requested \code{times} for a
#' model fitted with \code{fit_fastgbm}, for any of its objectives.
#'
#' @param object An \code{"mlsurv_model"} from \code{fit_fastgbm}.
#' @param newdata A \code{data.frame} of predictors.
#' @param times Numeric vector of evaluation time points.
#'
#' @return A \code{data.frame} of survival probabilities with one row per
#'   observation in \code{newdata} and one column per time point (columns named
#'   \code{"t=<time>"}).
#'
#' @seealso \code{fit_fastgbm}
#'
#' @examplesIf requireNamespace("fastgbm", quietly = TRUE)
#' mod <- fit_fastgbm(Surv(time, status) ~ age + karno, veteran, ntrees = 50)
#' predict_fastgbm(mod, veteran[1:5, ], times = c(100, 200, 300))
#'
#' @keywords internal
#' @noRd
predict_fastgbm <- function(object, newdata, times) {
  if (!is.null(object$learner) && object$learner != "fastgbm") {
    warning("Object passed to predict_fastgbm() may not come from fit_fastgbm().")
  }
  stopifnot(requireNamespace("fastgbm", quietly = TRUE))
  newdata <- as.data.frame(newdata)
  parsed <- .parse_surv_formula(object$formula, object$data)
  # The engine rebuilds the model frame from the formula, so the outcome columns
  # must exist; add dummies when newdata carries predictors only.
  for (col in c(parsed$time_col, parsed$status_col)) {
    if (!col %in% names(newdata)) newdata[[col]] <- 1
  }
  sv <- function(t) stats::predict(object$model, newdata = newdata, type = "survival", times = t)
  # fastgbm < 0.7.5 scrambles multi-time "aft" predictions; one time at a time is correct
  pred <- if (identical(object$objective, "aft")) vapply(times, sv, numeric(nrow(newdata))) else sv(times)
  pred <- matrix(pred, nrow = nrow(newdata))
  .finalize_survmat(pred, times = times)
}

#' Tune fastgbm survival hyperparameters by cross-validation
#'
#' Scores every row of `param_grid` by `folds`-fold cross-validation and returns
#' the grid with one column per metric, best first. The default grid crosses the
#' number of boosting rounds, the learning rate, the tree depth and the minimum
#' leaf size: `ntrees` 100, 300; `learning_rate` 0.03, 0.1; `max_depth` 3, 5;
#' `min_node_size` 10, 30 (16 configurations). Any other `fastgbm::fastgbm()`
#' argument can be a grid column (for example `max_leaves`, `subsample`,
#' `lambda`, `pexp_bins`), and fixed arguments such as `objective` go in `...`.
#'
#' @param formula,data,times Survival formula, data and evaluation times.
#' @param param_grid Data frame of configurations (one row each).
#' @param metrics Metrics to report; the first orders the result.
#' @param folds,seed Number of folds and fold seed.
#' @param ncores Workers for the fold loop.
#' @param ... Fixed arguments passed to `fit_fastgbm()`.
#' @return A data frame: the grid columns and one column per metric.
#' @keywords internal
#' @noRd
tune_fastgbm <- function(formula, data, times,
                         param_grid = expand.grid(ntrees = c(100, 300), learning_rate = c(0.03, 0.1),
                                                  max_depth = c(3, 5), min_node_size = c(10, 30)),
                         metrics = c("cindex", "ibs"), folds = 5, seed = 123, ncores = 1, ...) {
  # spell out a `.` so cv_survlearner() does not warn about it
  formula <- stats::formula(stats::terms(formula, data = data))
  fixed <- list(...)
  res <- .pmap_rbind_dt(param_grid, function(...) {
    p <- list(...)
    cvr <- do.call(cv_survlearner, c(list(formula = formula, data = data, fit_fun = fit_fastgbm,
                                          pred_fun = predict_fastgbm, times = times, metrics = metrics,
                                          folds = folds, seed = seed, ncores = ncores), p, fixed))
    .wide_metric_row(p, cv_summary(cvr))
  }, learner = "fastgbm")
  .arrange_by_metric_dt(res, metrics[1], .metric_maximize(metrics[1]))
}
