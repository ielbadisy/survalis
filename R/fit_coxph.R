

#' Fit a Cox Proportional Hazards Model
#'
#' Fits a Cox proportional hazards regression model using the
#' `survival::coxph()` function, and returns an object compatible
#' with the `mlsurv_model` interface.
#'
#' @param formula A survival formula of the form
#'   \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing the variables in the model.
#' @param ... Additional arguments passed to \code{survival::coxph()}.
#'
#' @details
#' The fitted object is stored along with metadata such as the
#' learner name (\code{"coxph"}), original formula, data, and names
#' of the time and status variables. The function requires the
#' \pkg{survival} and \pkg{pec} packages.
#'
#' @return An object of class \code{"mlsurv_model"} containing:
#' \itemize{
#'   \item \code{model} -- the fitted \code{coxph} object
#'   \item \code{learner} -- the string \code{"coxph"}
#'   \item \code{formula} -- the survival formula used
#'   \item \code{data} -- the training dataset
#'   \item \code{time}, \code{status} -- names of the survival outcome variables
#' }
#'
#' @examples
#' mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' summary(mod_cox)
#'
#' @export

fit_coxph <- function(formula, data, ...) {
  stopifnot(requireNamespace("survival", quietly = TRUE))
  stopifnot(requireNamespace("pec", quietly = TRUE))

  model <- survival::coxph(formula, data = data, x = TRUE, y = TRUE, ...)

  time_status <- all.vars(formula[[2]])
  structure(list(
    model = model,
    learner = "coxph",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model", engine = "survival")
}

#' Predict Survival Probabilities from a Cox PH Model
#'
#' Generates predicted survival probabilities at specified time points
#' for new data, using a fitted Cox proportional hazards model.
#'
#' @param object An \code{"mlsurv_model"} object returned by
#'   \code{fit_coxph()}.
#' @param newdata A data frame containing the predictor variables for
#'   which to compute predictions.
#' @param times A numeric vector of time points at which to evaluate
#'   survival probabilities.
#'
#' @details
#' Predictions are computed using \code{pec::predictSurvProb()}.
#' The output is formatted as a data frame with one row per observation
#' in \code{newdata} and one column per time point.
#'
#' @return A data frame of survival probabilities with columns named
#' \code{"t=<time>"} for each requested time.
#'
#' @examples
#' mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' predict_coxph(mod_cox, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export

predict_coxph <- function(object, newdata, times) {

  if (!is.null(object$learner) && object$learner != "coxph") {warning("Object passed to predict_coxph() may not come from fit_coxph().")}

  stopifnot(requireNamespace("pec", quietly = TRUE))

  survmat <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}

