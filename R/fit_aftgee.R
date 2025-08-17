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
#' @examples
#' \donttest{
#' mod <- fit_aftgee(
#'   Surv(time, status) ~ trt + karno + age,
#'   data = veteran
#' )
#' head(predict_aftgee(mod, newdata = veteran[1:5, ], times = c(20, 60, 120)))
#' }

#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @export

fit_aftgee <- function(formula, data, corstr = "independence") {
  stopifnot(requireNamespace("aftgee", quietly = TRUE))


  model <- aftgee::aftgee(
    formula = formula,
    data = data,
    id = NULL, ## explicitly no clustering supported <=> one observation per subject
    corstr = corstr
  )

  structure(list(
    model = model,
    learner = "aftgee",
    formula = formula,
    data = data
  ), class = "mlsurv_model", engine = "aftgee")
}


#' Predict Survival Probabilities from an \code{aftgee} Model
#'
#' Computes survival probabilities \eqn{S(t \mid x)} at specified time points from a fitted
#' AFT model using a log-normal approximation for the error distribution.
#'
#' @param object An \code{mlsurv_model} produced by \code{\link{fit_aftgee}}.
#' @param newdata A data.frame of predictor values for prediction.
#' @param times Numeric vector of time points at which to estimate survival probabilities.
#'   If \code{NULL}, a sequence of 20 values from \code{1} to the maximum observed time in
#'   \code{object$data} is used.
#'
#' @details
#' We use the approximation
#' \deqn{S(t \mid x) \approx 1 - \Phi\!\left(\frac{\log t - x^\top \hat{\beta}}{\sigma}\right),}
#' where \eqn{\Phi} is the standard normal CDF. Here we set \eqn{\sigma = 1} as a simple,
#' distribution-agnostic proxy; this yields a monotone-in-time score useful for benchmarking.
#' For production use, prefer a parametric AFT fit where \eqn{\sigma} is estimated.
#'
#' @return A data.frame of survival probabilities with one row per \code{newdata} observation
#'   and one column per requested time (\code{"t=<time>"}).
#'
#' @examples
#'   mod <- fit_aftgee(Surv(time, status) ~ trt + karno + age, data = veteran)
#'   predict_aftgee(mod, newdata = veteran[1:5, ], times = c(20, 60, 120))
#'
#' @importFrom stats delete.response model.matrix terms pnorm
#' @export

predict_aftgee <- function(object, newdata, times = NULL) {


  if (!is.null(object$learner) && object$learner != "aftgee") {
    warning("Object passed to predict_aftgee() may not come from fit_aftgee().")
  }

  stopifnot(requireNamespace("aftgee", quietly = TRUE))

  beta <- object$model$coef.res
  X <- model.matrix(delete.response(terms(object$formula)), newdata)
  log_time_pred <- as.numeric(X %*% beta)

  if (is.null(times)) {
    t_max <- max(newdata[[all.vars(object$formula[[2]])[1]]], na.rm = TRUE)
    times <- seq(1, t_max, length.out = 20)
  }

  # approximate survival using lognormal-like assumption
  survmat <- sapply(times, function(t) {
    1 - pnorm((log(t) - log_time_pred))  # approximate S(t) from log-T model
  })

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}

