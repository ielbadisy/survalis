#' Fit an Additive Hazards (Aalen) Model
#'
#' Fits an additive hazards regression model using \pkg{timereg}'s
#' \code{\link[timereg]{aalen}} and returns an \code{mlsurv_model} compatible
#' with the \pkg{survalis} workflow.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param max.time Optional maximum follow-up time used by the fitting routine.
#' @param n.sim Integer; number of simulations for variance estimation (default \code{0}).
#' @param resample.iid Integer; indicator for iid resampling (passed to \code{aalen}).
#'
#' @details
#' The Aalen model assumes an additive hazard:
#' \deqn{\lambda(t \mid X) = \beta_0(t) + X^\top \beta(t),}
#' with nonparametric cumulative coefficient functions.
#'
#' @return A list of class \code{"mlsurv_model"} with elements:
#' \code{model}, \code{learner="aalen"}, \code{engine="timereg"}, \code{formula},
#' \code{data}, \code{time}, \code{status}.
#'
#' @examples
#' \donttest{
#'   mod_aalen <- fit_aalen(
#'     Surv(time, status) ~ trt + karno + age,
#'     data = veteran
#'   )
#'   head(predict_aalen(mod_aalen, newdata = veteran[1:5, ], times = c(50, 100, 150)))
#' }
#' @export
fit_aalen <- function(formula, data, max.time = NULL, n.sim = 0, resample.iid = 1) {
  stopifnot(requireNamespace("timereg", quietly = TRUE))

  model <- timereg::aalen(
    formula = formula,
    data = data,
    max.time = max.time,
    n.sim = n.sim,
    resample.iid = resample.iid
  )

  time_status <- all.vars(formula[[2]]) # c(time, status)

  structure(list(
    model   = model,
    learner = "aalen",
    engine  = "timereg",
    formula = formula,
    data    = data,
    time    = time_status[1],
    status  = time_status[2]
  ), class = "mlsurv_model")
}

#' Predict Survival from an Aalen Additive Hazards Model
#'
#' Computes survival probabilities at specified time points from a model fitted
#' with \code{\link{fit_aalen}}.
#'
#' @param object An \code{"mlsurv_model"} returned by \code{\link{fit_aalen}}.
#' @param newdata A data frame of new observations.
#' @param times Numeric vector of time points at which to evaluate survival probabilities.
#'
#' @return A data frame (rows = observations, columns = \code{"t=<time>"}).
#'
#' @examples
#' \donttest{
#'   mod <- fit_aalen(Surv(time, status) ~ trt + karno + age, data = veteran, max.time = 600)
#'   head(predict_aalen(mod, newdata = veteran[1:5, ], times = 0:10))
#' }
#' @export
predict_aalen <- function(object, newdata, times) {
  stopifnot(inherits(object, "mlsurv_model"))
  stopifnot(identical(object$learner, "aalen"))

  pred <- timereg::predict.aalen(
    object$model,
    newdata = as.data.frame(newdata),
    times   = times
  )

  if (is.null(pred$S0))
    stop("timereg::predict.aalen() did not return S0 (survival probabilities).")

  surv_probs <- pmin(pmax(pred$S0, 0), 1)
  colnames(surv_probs) <- paste0("t=", pred$time)
  as.data.frame(surv_probs)
}


