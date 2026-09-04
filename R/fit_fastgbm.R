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
#'   \code{"cox"}, or \code{"aft"}. Only \code{"pexp"} supports
#'   \code{predict(type = "survival")} at arbitrary times and is required for the
#'   survalis prediction contract.
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
#' @seealso \code{\link{predict_fastgbm}}, \code{\link{tune_fastgbm}},
#'   \code{\link[fastgbm]{fastgbm}}
#'
#' @examplesIf requireNamespace("fastgbm", quietly = TRUE)
#' mod <- fit_fastgbm(Surv(time, status) ~ age + karno + celltype, veteran,
#'                    ntrees = 50, pexp_bins = 8)
#' summary(mod)
#'
#' @export
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
#' model fitted with \code{\link{fit_fastgbm}} (piecewise-exponential objective).
#'
#' @param object An \code{"mlsurv_model"} from \code{\link{fit_fastgbm}}.
#' @param newdata A \code{data.frame} of predictors.
#' @param times Numeric vector of evaluation time points.
#'
#' @return A \code{data.frame} of survival probabilities with one row per
#'   observation in \code{newdata} and one column per time point (columns named
#'   \code{"t=<time>"}).
#'
#' @seealso \code{\link{fit_fastgbm}}
#'
#' @examplesIf requireNamespace("fastgbm", quietly = TRUE)
#' mod <- fit_fastgbm(Surv(time, status) ~ age + karno, veteran, ntrees = 50)
#' predict_fastgbm(mod, veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export
predict_fastgbm <- function(object, newdata, times) {
  if (!is.null(object$learner) && object$learner != "fastgbm") {
    warning("Object passed to predict_fastgbm() may not come from fit_fastgbm().")
  }
  stopifnot(requireNamespace("fastgbm", quietly = TRUE))
  if (!is.null(object$objective) && object$objective != "pexp") {
    stop("predict_fastgbm() needs objective = 'pexp'; refit with fit_fastgbm(..., objective = 'pexp').")
  }

  pred <- stats::predict(object$model, newdata = newdata, type = "survival", times = times)
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)
  .finalize_survmat(pred, times = times)
}
