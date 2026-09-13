#' Fit a Conditional Naive Bayes Survival Model
#'
#' Fits a conditional naive Bayes survival model with the \pkg{nbsurv}
#' package and returns an object compatible with the `mlsurv_model`
#' interface. Frames survival at each horizon as a binary classification
#' problem ("survive beyond t" vs. "fail before t") and applies Bayes' rule
#' under a conditional-independence assumption over predictors, with IPCW
#' weighting for censoring.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing the formula variables.
#' @param scale Logical; standardize continuous predictors using training
#'   mean/sd (default \code{TRUE}).
#' @param laplace Laplace smoothing constant for categorical predictors
#'   (default \code{1}).
#' @param min_sd Floor for estimated SDs of Gaussian class-conditional
#'   distributions (default \code{0.05}).
#' @param ... Additional arguments passed to \code{nbsurv::nbsurv()} (for
#'   example \code{time_grid}, \code{cov_structure}, \code{shrinkage},
#'   \code{time_smooth}, \code{bandwidth}).
#'
#' @return An object of class \code{"mlsurv_model"} with elements \code{model},
#'   \code{learner} (\code{"nbsurv"}), \code{formula}, \code{data}, \code{time},
#'   and \code{status}.
#'
#' @seealso \code{\link{predict_nbsurv}}, \code{\link[nbsurv]{nbsurv}},
#'   \code{\link{fit}} for the verb interface
#'
#' @examplesIf requireNamespace("nbsurv", quietly = TRUE)
#' mod <- fit_nbsurv(Surv(time, status) ~ age + karno + celltype, veteran)
#' summary(mod)
#'
#' @keywords internal
#' @export
fit_nbsurv <- function(formula, data, scale = TRUE, laplace = 1, min_sd = 0.05, ...) {
  stopifnot(requireNamespace("nbsurv", quietly = TRUE))

  model <- nbsurv::nbsurv(
    formula, data = data, scale = scale, laplace = laplace, min_sd = min_sd, ...
  )

  time_status <- all.vars(formula[[2]])
  structure(list(
    model = model,
    learner = "nbsurv",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model", engine = "nbsurv")
}

#' Predict Survival Probabilities from a nbsurv Model
#'
#' Generates predicted survival probabilities at the requested \code{times}
#' for a model fitted with \code{\link{fit_nbsurv}}.
#'
#' @param object An \code{"mlsurv_model"} from \code{\link{fit_nbsurv}}.
#' @param newdata A \code{data.frame} of predictors.
#' @param times Numeric vector of evaluation time points.
#'
#' @return A \code{data.frame} of survival probabilities with one row per
#'   observation in \code{newdata} and one column per time point (columns named
#'   \code{"t=<time>"}).
#'
#' @seealso \code{\link{fit_nbsurv}}
#'
#' @examplesIf requireNamespace("nbsurv", quietly = TRUE)
#' mod <- fit_nbsurv(Surv(time, status) ~ age + karno, veteran)
#' predict_nbsurv(mod, veteran[1:5, ], times = c(100, 200, 300))
#'
#' @keywords internal
#' @export
predict_nbsurv <- function(object, newdata, times) {
  if (!is.null(object$learner) && object$learner != "nbsurv") {
    warning("Object passed to predict_nbsurv() may not come from fit_nbsurv().")
  }
  stopifnot(requireNamespace("nbsurv", quietly = TRUE))

  pred <- stats::predict(object$model, newdata = newdata, times = times, type = "survival")
  .finalize_survmat(pred, times = times)
}
