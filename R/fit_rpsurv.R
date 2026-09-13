#' Fit a Royston-Parmar Flexible Parametric Survival Model
#'
#' Fits a Royston-Parmar flexible parametric survival model with the
#' \pkg{rpsurv} package and returns an object compatible with the
#' `mlsurv_model` interface. Restricted cubic splines model the baseline on
#' the chosen \code{scale}, optionally with time-varying covariate effects
#' (\code{tve}), so predicted survival curves are available at arbitrary
#' time points without a proportional-hazards assumption.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing the formula variables.
#' @param df Degrees of freedom for the baseline spline (default \code{4}).
#' @param knots Optional explicit knot vector (log time scale); default
#'   \code{NULL} derives knots from \code{df}.
#' @param tve Optional character vector of covariate names given
#'   time-varying effects via their own restricted cubic spline.
#' @param tve.df Degrees of freedom for each time-varying-effect spline
#'   (default \code{3}).
#' @param scale One of \code{"hazard"} (proportional hazards, default),
#'   \code{"odds"}, or \code{"normal"}.
#' @param ... Additional arguments passed to \code{rpsurv::rpsurv()} (for
#'   example \code{control}).
#'
#' @return An object of class \code{"mlsurv_model"} with elements \code{model},
#'   \code{learner} (\code{"rpsurv"}), \code{formula}, \code{data}, \code{time},
#'   \code{status}, plus \code{rhs_terms}/\code{cov_terms}/\code{xlev} used to
#'   rebuild the design matrix for arbitrary \code{newdata} at predict time.
#'
#' @seealso \code{\link{predict_rpsurv}}, \code{\link[rpsurv]{rpsurv}},
#'   \code{\link{fit}} for the verb interface
#'
#' @examplesIf requireNamespace("rpsurv", quietly = TRUE)
#' mod <- fit_rpsurv(Surv(time, status) ~ age + karno + celltype, veteran, df = 3)
#' summary(mod)
#'
#' @keywords internal
#' @export
fit_rpsurv <- function(formula, data, df = 4, knots = NULL, tve = NULL,
                       tve.df = 3, scale = c("hazard", "odds", "normal"), ...) {
  stopifnot(requireNamespace("rpsurv", quietly = TRUE))
  scale <- match.arg(scale)

  model <- rpsurv::rpsurv(
    formula, data = data, df = df, knots = knots, tve = tve, tve.df = tve.df,
    scale = scale, ...
  )

  time_status <- all.vars(formula[[2]])
  rhs_terms <- stats::delete.response(stats::terms(formula))
  cov_terms <- attr(rhs_terms, "term.labels")

  # rpsurv::predict.rpsurv() needs `newdata` already expanded to the same
  # model.matrix columns rpsurv() built internally (e.g. factor dummies), not
  # the raw covariate columns. Record each raw factor's training levels so
  # predict_rpsurv() can rebuild that expansion identically for any newdata.
  rhs_mf <- stats::model.frame(rhs_terms, data)
  xlev <- lapply(rhs_mf, function(col) if (is.factor(col)) levels(col) else NULL)
  xlev <- xlev[!vapply(xlev, is.null, logical(1))]

  structure(list(
    model = model,
    learner = "rpsurv",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2],
    rhs_terms = rhs_terms,
    cov_terms = cov_terms,
    xlev = xlev
  ), class = "mlsurv_model", engine = "rpsurv")
}

#' Predict Survival Probabilities from a rpsurv Model
#'
#' Generates predicted survival probabilities at the requested \code{times}
#' for a model fitted with \code{\link{fit_rpsurv}}.
#'
#' @param object An \code{"mlsurv_model"} from \code{\link{fit_rpsurv}}.
#' @param newdata A \code{data.frame} of predictors (raw covariate columns,
#'   not pre-expanded dummies).
#' @param times Numeric vector of evaluation time points.
#'
#' @return A \code{data.frame} of survival probabilities with one row per
#'   observation in \code{newdata} and one column per time point (columns named
#'   \code{"t=<time>"}).
#'
#' @seealso \code{\link{fit_rpsurv}}
#'
#' @examplesIf requireNamespace("rpsurv", quietly = TRUE)
#' mod <- fit_rpsurv(Surv(time, status) ~ age + karno, veteran, df = 3)
#' predict_rpsurv(mod, veteran[1:5, ], times = c(100, 200, 300))
#'
#' @keywords internal
#' @export
predict_rpsurv <- function(object, newdata, times) {
  if (!is.null(object$learner) && object$learner != "rpsurv") {
    warning("Object passed to predict_rpsurv() may not come from fit_rpsurv().")
  }
  stopifnot(requireNamespace("rpsurv", quietly = TRUE))

  mf_new <- stats::model.frame(object$rhs_terms, newdata, xlev = object$xlev)
  cov_expanded <- as.data.frame(
    stats::model.matrix(stats::reformulate(object$cov_terms, intercept = FALSE), mf_new)
  )

  long <- stats::predict(object$model, newdata = cov_expanded, times = times, type = "survival")
  pred <- matrix(long$est, nrow = nrow(newdata), ncol = length(times), byrow = TRUE)
  .finalize_survmat(pred, times = times)
}
