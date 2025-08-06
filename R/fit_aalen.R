fit_aalen <- function(formula, data, max.time = NULL, n.sim = 0, resample.iid = 1, ...) {

  stopifnot(requireNamespace("timereg", quietly = TRUE))

  model <- timereg::aalen(
    formula = formula,
    data = data,
    max.time = max.time,
    n.sim = n.sim,
    resample.iid = resample.iid,
    ...
  )

  time_status <- all.vars(formula[[2]])

  structure(list(
    model = model,
    learner = "aalen",
    engine = "timereg",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model")
}

predict_aalen <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "aalen")

  pout <- predict(object$model, newdata = newdata, times = times, ...)

  if (is.null(pout$S0)) {
    stop("Survival probabilities (S0) not returned by predict().")
  }

  survmat <- pout$S0
  colnames(survmat) <- paste0("t=", pout$time)
  as.data.frame(survmat)
}


predict_aalen <- function(object, newdata, times, ...) {
  stopifnot(inherits(object, "mlsurv_model"))
  stopifnot(object$learner == "aalen")
  stopifnot(requireNamespace("pec", quietly = TRUE))

  survmat <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  colnames(survmat) <- paste0("t=", times)
  rownames(survmat) <- NULL

  # Add time 0 = 1 as initial survival
  survmat <- cbind("t=0" = 1, survmat)
  survmat
}


#' Predict survival probabilities from Aalen Additive Hazards Model
#'
#' @param object Fitted mlsurv_model from fit_aalen()
#' @param newdata New data to predict on
#' @param times Vector of time points
#' @param ... Additional arguments
#' @return A data frame of survival probabilities
#' @export
predict_aalen <- function(object, newdata, times, ...) {
  stopifnot(inherits(object, "mlsurv_model"))
  stopifnot(object$learner == "aalen")

  pred <- timereg::predict.aalen(
    object$model,
    newdata = newdata,
    times = times,
    ...
  )

  if (is.null(pred$S0)) {
    stop("S0 (survival probability) not found in timereg::predict.aalen output.")
  }

  # Clip probabilities to [0, 1]
  surv_probs <- pmin(pmax(pred$S0, 0), 1)
  colnames(surv_probs) <- paste0("t=", pred$time)
  as.data.frame(surv_probs)
}


library(survival)
library(timereg)

veteran <- survival::veteran
times <- 0:100

mod <- fit_aalen(Surv(time, status) ~ trt + karno + age,
                 data = veteran, max.time = 600)

pred <- predict_aalen(mod, newdata = veteran[1:100, ], times = times)
pred


cv_survlearner(
  formula = Surv(time, status) ~ trt + karno + age,
  data = veteran,
  fit_fun = fit_aalen,
  pred_fun = predict_aalen,
  times = times,
  metrics = c("cindex", "ibs"),
  folds = 5
)


