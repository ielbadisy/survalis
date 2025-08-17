#' Fit a Flexible Parametric Survival Model (rstpm2, mlsurv_model-compatible)
#'
#' Fits a parametric/penalised generalised survival model using
#' \pkg{rstpm2}'s \code{stpm2()}, returning an `mlsurv_model` object that
#' integrates with the `survalis` evaluation and cross‑validation pipeline.
#' The baseline hazard is modeled with restricted cubic splines controlled by
#' the degrees of freedom.
#'
#' @section Engine:
#' Uses \strong{rstpm2::stpm2}. Spline complexity is governed by \code{df};
#' larger values allow more flexibility in the baseline hazard. Additional
#' engine arguments can be passed via \code{...}.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left‑hand side must be a `Surv()` object from the \pkg{survival} package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param df Integer degrees of freedom for the restricted cubic spline baseline.
#'   Default is `3`.
#' @param ... Additional arguments forwarded to \code{rstpm2::stpm2()}.
#'
#' @return An object of class `mlsurv_model`, a named list with elements:
#' \describe{
#'   \item{model}{The fitted \code{rstpm2::stpm2} object.}
#'   \item{learner}{Character scalar identifying the learner (`"stpm2"`).}
#'   \item{engine}{Character scalar naming the engine (`"rstpm2"`).}
#'   \item{formula}{The original survival formula.}
#'   \item{data}{The training dataset (or a minimal subset needed for prediction).}
#'   \item{time}{Name of the survival time variable.}
#'   \item{status}{Name of the event indicator (1 = event, 0 = censored).}
#' }
#'
#' @details
#' This wrapper standardizes the model object to the `mlsurv_model` contract so
#' downstream prediction and evaluation behave consistently across learners.
#'
#' @references
#' Royston P, Parmar MKB (2002). Flexible parametric proportional‑hazards and
#' proportional‑odds models for censored survival data. \emph{Statistics in Medicine}.  
#' Lambert PC, Royston P, Crowther MJ (2009). Restricted cubic splines for
#' non‑proportional hazards. \emph{Statistics in Medicine}.
#'
#' @examples
#'
#' mod_stpm2 <- fit_stpm2(Surv(time, status) ~ age + karno + celltype,
#'                        data = veteran, df = 4)
#' predict_stpm2(mod_stpm2, newdata = veteran[1:5, ], times = c(100, 200, 300))
#' summary(mod_stpm2)
#'
#' @export

fit_stpm2 <- function(formula, data, df = 4, ...) {
  stopifnot(requireNamespace("rstpm2", quietly = TRUE))

  # Inject nsx into formula environment so it's always found
  fenv <- new.env(parent = environment(formula))
  fenv$nsx <- getFromNamespace("nsx", "rstpm2")
  environment(formula) <- fenv

  model <- rstpm2::stpm2(formula = formula, data = data, df = df, ...)
  
  time_status <- all.vars(formula[[2]])
  structure(list(
    model = model,
    learner = "stpm2",
    engine = "rstpm2",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model")
}

#' Predict Survival Probabilities with an rstpm2 Model
#'
#' Generates survival probabilities at specified time points using a fitted
#' \pkg{rstpm2} model wrapped as an `mlsurv_model`.
#'
#' @param object A fitted `mlsurv_model` returned by [fit_stpm2()].
#' @param newdata A `data.frame` of new observations for prediction. Must include
#'   the same predictor variables used at training time.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   training survival time). Must be non‑negative and finite.
#' @param ... Additional arguments forwarded to \code{predict()}.
#'
#' @return A base `data.frame` with one row per observation in `newdata` and
#' columns named `t={time}` (character), containing numeric survival
#' probabilities in \[0, 1].
#'
#' @details
#' Internally expands \code{newdata} over the requested \code{times} and calls
#' \code{predict(object$model, type = "surv", newtime = times)}; results are
#' reshaped to a wide format with one column per time point.
#'
#' @seealso [fit_stpm2()]
#'
#' @examples
#'
#' mod_stpm2 <- fit_stpm2(Surv(time, status) ~ age + karno + celltype,
#'                        data = veteran, df = 4)
#' predict_stpm2(mod_stpm2, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export

predict_stpm2 <- function(object, newdata, times, ...) {
  if (!inherits(object, "mlsurv_model") || is.null(object$learner) || object$learner != "stpm2") {
    stop("`object` must be an mlsurv_model from fit_stpm2().")
  }
  stopifnot(requireNamespace("rstpm2", quietly = TRUE))

  newdata <- as.data.frame(newdata)
  Tvar <- object$time  # usually "time"

  # Build a long newdata with a `time` column
  nd_long <- newdata[rep(seq_len(nrow(newdata)), each = length(times)), , drop = FALSE]
  nd_long[[Tvar]] <- rep(times, times = nrow(newdata))

  # Predict survival on that grid
  # predict.stpm2 will dispatch here because we call rstpm2::predict()
  surv_vec <- rstpm2::predict(object$model, newdata = nd_long, type = "surv", se.fit = FALSE, ...)

  # Reshape back to nrow(newdata) x length(times)
  surv_mat <- matrix(as.numeric(surv_vec), nrow = nrow(newdata), byrow = TRUE)
  colnames(surv_mat) <- paste0("t=", times)
  as.data.frame(surv_mat)
}
