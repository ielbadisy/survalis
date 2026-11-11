#' Fit an Accelerated Failure Time Model Using Generalized Estimating Equations
#'
#' Fits an accelerated failure time (AFT) model via \pkg{aftgee}, which implements
#' GEE methodology for censored survival data. Returns an \code{mlsurv_model}-compatible
#' object used throughout the package.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A data.frame containing the variables in the model.
#' @param corstr Working correlation structure; one of \code{"independence"},
#'   \code{"exchangeable"}, or \code{"ar1"}. Default: \code{"independence"}.
#'
#' @details
#' The AFT model assumes \eqn{\log(T) = X \beta + \epsilon}, where \eqn{T} is survival time,
#' \eqn{X} is the covariate matrix, and \eqn{\epsilon} is an error term whose distribution
#' determines the baseline survival. \pkg{aftgee} estimates \eqn{\beta} using GEE.
#' In this wrapper we set \code{id = NULL}, assuming one observation per subject.
#'
#' @return An object of class \code{mlsurv_model} with components:
#' \itemize{
#'   \item \code{model}: the fitted \code{aftgee} model.
#'   \item \code{learner}: \code{"aftgee"}.
#'   \item \code{formula}: the survival formula.
#'   \item \code{data}: the training data.
#' }
#'
#' @references
#' Jin, Z., Lin, D. Y., Wei, L. J., & Ying, Z. (2003). Rank-based inference for the
#' accelerated failure time model. \emph{Biometrika}, 90(2), 341-353.
#'
#' @examplesIf requireNamespace("aftgee", quietly = TRUE)
#' \donttest{
#' mod <- fit_aftgee(
#'   Surv(time, status) ~ trt + karno + age,
#'   data = veteran
#' )
#' head(predict_aftgee(mod, newdata = veteran[1:5, ], times = c(20, 60, 120)))
#' }

#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @keywords internal
#' @noRd

fit_aftgee <- function(formula, data, corstr = "independence") {
  stopifnot(requireNamespace("aftgee", quietly = TRUE))

  # expand `.` and record the training encoding (terms, factor levels,
  # contrasts, columns) so predict_aftgee() rebuilds the same design matrix
  design <- .survalis_design_fit(formula, data, intercept = TRUE)
  formula <- design$formula

  model <- aftgee::aftgee(
    formula = formula,
    data = data,
    id = NULL, ## explicitly no clustering supported <=> one observation per subject
    corstr = corstr
  )

  out <- structure(list(
    model = model,
    learner = "aftgee",
    formula = formula,
    data = data,
    design = design[c("terms", "xlevels", "x_contrasts", "x_cols")]
  ), class = "mlsurv_model", engine = "aftgee")
  out$baseline <- .aftgee_baseline(out)
  out
}

# Kaplan-Meier estimate of the error distribution from the censored
# residuals e = log(T) - x'beta on the training data. This is the
# semiparametric baseline of the AFT model (no distribution assumed).
.aftgee_baseline <- function(object) {
  # complete cases only, so the outcome and the linear predictor align
  d <- object$data[stats::complete.cases(object$data[, all.vars(object$formula), drop = FALSE]), , drop = FALSE]
  y <- stats::model.response(stats::model.frame(object$formula, d))
  lp <- .aftgee_lp(object, d)
  ok <- is.finite(lp) & y[, 1] > 0
  km <- survival::survfit(survival::Surv(log(y[ok, 1]) - lp[ok], y[ok, 2]) ~ 1)
  list(resid = km$time, surv = km$surv)
}

.aftgee_lp <- function(object, newdata) {
  beta <- object$model$coef.res
  if (!is.null(object$design)) {
    X <- .survalis_design_new(object$design, newdata)
    if (!is.null(names(beta)) && all(names(beta) %in% colnames(X))) {
      X <- X[, names(beta), drop = FALSE]
    }
  } else {
    X <- model.matrix(delete.response(terms(object$formula)), newdata)
  }
  as.numeric(X %*% beta)
}


#' Predict Survival Probabilities from an \code{aftgee} Model
#'
#' Computes survival probabilities \eqn{S(t \mid x)} at specified time points from a fitted
#' semiparametric AFT model.
#'
#' @param object An \code{mlsurv_model} produced by \code{fit_aftgee}.
#' @param newdata A data.frame of predictor values for prediction.
#' @param times Numeric vector of time points at which to estimate survival probabilities.
#'   If \code{NULL}, a sequence of 20 values from \code{1} to the maximum observed time in
#'   \code{object$data} is used.
#'
#' @details
#' With \eqn{\log T = x^\top \beta + \epsilon} and \eqn{\hat S_0} the Kaplan-Meier
#' estimate of the survival function of the censored training residuals
#' \eqn{\log T_i - x_i^\top \hat\beta},
#' \deqn{S(t \mid x) = \hat S_0(\log t - x^\top \hat{\beta}),}
#' read as a step function (1 below the smallest residual, constant beyond the largest).
#'
#' @return A data.frame of survival probabilities with one row per \code{newdata} observation
#'   and one column per requested time (\code{"t=<time>"}).
#'
#' @examplesIf requireNamespace("aftgee", quietly = TRUE)
#'   mod <- fit_aftgee(Surv(time, status) ~ trt + karno + age, data = veteran)
#'   predict_aftgee(mod, newdata = veteran[1:5, ], times = c(20, 60, 120))
#'
#' @importFrom stats delete.response model.matrix terms
#' @keywords internal
#' @noRd

predict_aftgee <- function(object, newdata, times = NULL) {


  if (!is.null(object$learner) && object$learner != "aftgee") {
    warning("Object passed to predict_aftgee() may not come from fit_aftgee().")
  }

  stopifnot(requireNamespace("aftgee", quietly = TRUE))

  log_time_pred <- .aftgee_lp(object, newdata)

  if (is.null(times)) {
    t_max <- max(object$data[[all.vars(object$formula[[2]])[1]]], na.rm = TRUE)
    times <- seq(1, t_max, length.out = 20)
  }

  # models fitted before the KM baseline was stored: rebuild it
  base <- object$baseline
  if (is.null(base)) base <- .aftgee_baseline(object)

  survmat <- vapply(times, function(t) {
    k <- findInterval(log(t) - log_time_pred, base$resid)
    c(1, base$surv)[k + 1L]
  }, numeric(length(log_time_pred)))
  if (!is.matrix(survmat)) survmat <- matrix(survmat, ncol = length(times))

  .finalize_survmat(survmat, times = times)
}
