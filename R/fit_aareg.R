#' Fit an Aalen Additive Survival Model (mlsurv_model-compatible)
#'
#' @param formula A survival formula
#' @param data Data frame
#' @param max.time Max follow-up time (optional)
#' @param n.sim Number of simulations
#' @param resample.iid Type of resampling (default = 1)
#' @param ... Passed to timereg::aalen
#' @return A mlsurv_model object
#' @export
fit_aareg <- function(formula, data, max.time = NULL, n.sim = 0, resample.iid = 1, ...) {
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
    learner = "aareg",
    engine = "timereg",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model")
}


#' Predict Survival Probabilities for Aalen Model
#'
#' @param object A mlsurv_model object from fit_aareg()
#' @param newdata New data to predict on
#' @param times Vector of times
#' @param ... Ignored
#' @return A data.frame of survival probabilities
#' @export
predict_aareg <- function(object, newdata, times, ...) {
  stopifnot(inherits(object, "mlsurv_model"))
  stopifnot(object$learner == "aareg")

  pout <- predict(object$model, newdata = newdata, times = times, ...)

  if (is.null(pout$S0)) {
    stop("Survival probabilities (S0) not returned by predict().")
  }

  surv_probs <- pout$S0
  colnames(surv_probs) <- paste0("t=", pout$time)
  #rownames(survmat) <- paste0("id", seq_len(nrow(newdata)))
  as.data.frame(surv_probs)
}



library(survival)
library(timereg)

data(veteran)
times <- c(0, 100, 300, 500)

mod <- fit_aareg(Surv(time, status) ~ trt + karno + age,
                 data = veteran, max.time = 600)

pred <- predict_aareg(mod, newdata = veteran[1:5, ], times = times)
print(round(pred, 3))


cv_survlearner(
  formula = Surv(time, status) ~ trt + karno + age,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = times,
  metrics = c("cindex", "ibs"),
  folds = 5
)


