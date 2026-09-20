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
#' @keywords internal
#' @export
fit_aalen <- function(formula, data, max.time = NULL, n.sim = 0, resample.iid = 1) {
  stopifnot(requireNamespace("timereg", quietly = TRUE))

  time_status <- all.vars(formula[[2]]) # c(time, status)

  # timereg::aalen() silently returns an all-zero fit when the design matrix is badly scaled
  # (for example a covariate near 1000 next to a 0/1 covariate): the predicted survival is then
  # 1 for everyone. Numeric covariates are therefore standardized inside the learner, and the
  # same scaling is applied in predict_aalen().
  scaling <- .aalen_scaling(formula, data, time_status)
  model <- timereg::aalen(
    formula = formula,
    data = .aalen_scale(data, scaling),
    max.time = max.time,
    n.sim = n.sim,
    resample.iid = resample.iid
  )

  structure(list(
    model   = model,
    learner = "aalen",
    engine  = "timereg",
    formula = formula,
    data    = data,
    time    = time_status[1],
    status  = time_status[2],
    scaling = scaling
  ), class = "mlsurv_model")
}

.aalen_scaling <- function(formula, data, time_status) {
  vars <- if ("." %in% all.vars(formula)) names(data) else all.vars(formula)
  vars <- setdiff(intersect(vars, names(data)), time_status)
  vars <- vars[vapply(data[vars], is.numeric, logical(1))]
  center <- vapply(data[vars], function(v) mean(v, na.rm = TRUE), numeric(1))
  scale <- vapply(data[vars], function(v) stats::sd(v, na.rm = TRUE), numeric(1))
  keep <- is.finite(scale) & scale > 0
  list(center = center[keep], scale = scale[keep])
}

.aalen_scale <- function(df, scaling) {
  df <- as.data.frame(df)
  for (v in intersect(names(scaling$center), names(df))) {
    df[[v]] <- (df[[v]] - scaling$center[[v]]) / scaling$scale[[v]]
  }
  df
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
#' @keywords internal
#' @export
predict_aalen <- function(object, newdata, times) {
  stopifnot(inherits(object, "mlsurv_model"))
  stopifnot(identical(object$learner, "aalen"))

  pred <- timereg::predict.aalen(
    object$model,
    newdata = .aalen_scale(newdata, object$scaling),
    times   = times
  )

  if (is.null(pred$S0))
    stop("timereg::predict.aalen() did not return S0 (survival probabilities).")

  .finalize_survmat(pred$S0, times = pred$time)
}

