#' Fit a Survival Forest using aorsf
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data.frame containing the variables
#' @param ... Additional arguments passed to `aorsf::orsf()`
#'
#' @return A list of class 'mlsurv_model' with attributes for unified prediction
#' @export
fit_orsf <- function(formula, data, ...) {
  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  model <- aorsf::orsf(
    formula = formula,
    data = data,
    na_action = "omit",
    ...
  )

  structure(list(
    model = model,
    learner = "aorsf",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "aorsf")
}


#' Predict survival probabilities using ORSF model
#'
#' @param object An object from `fit_orsf()`
#' @param newdata New data for prediction
#' @param times Time points at which to predict survival probabilities
#' @param ... Additional arguments passed to `aorsf::predict()`
#'
#' @return A data.frame: rows = observations, columns = time points
#' @export
predict_orsf <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "aorsf")
  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  pred <- predict(
    object$model,
    new_data = newdata,
    pred_horizon = times,
    pred_type = "surv",
    na_action = "pass",
    boundary_checks = FALSE,
    ...
  )

  if (is.vector(pred)) {
    pred <- matrix(pred, nrow = 1)
  }

  colnames(pred) <- paste0("t=", times)
  rownames(pred) <- paste0("ID_", seq_len(nrow(newdata)))

  as.data.frame(pred)
}

library(aorsf)
library(survival)

mod_orsf <- fit_orsf(Surv(time, status) ~ age + ph.ecog, data = lung)
pred_orsf <- predict_orsf(mod_orsf, newdata = lung[1:5, ], times = c(100, 200, 300))
print(round(pred_orsf, 3))
