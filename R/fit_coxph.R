
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
  ), class = "mlsurv_model", engine = "coxph")
}


#' Predict Survival Probabilities from Cox Model
#'
#' @param object An object from fit_cox().
#' @param newdata A data.frame with predictors.
#' @param times A numeric vector of times at which to estimate survival.
#'
#' @return A data.frame of survival probabilities (rows: subjects, cols: times).
#' @export
predict_coxph <- function(object, newdata, times) {
  stopifnot(object$learner == "coxph")
  stopifnot(requireNamespace("pec", quietly = TRUE))

  survmat <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  colnames(survmat) <- paste0("t=", times)
  rownames(survmat) <- paste0("ID_", seq_len(nrow(newdata)))
  as.data.frame(survmat)
}


library(survival)
data(veteran)

mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred_cox <- predict_coxph(mod_cox, newdata = veteran[1:5, ], times = c(100, 200, 300))
print(round(pred_cox, 3))


