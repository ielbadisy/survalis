
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



predict_coxph <- function(object, newdata, times) {

  if (!is.null(object$learner) && object$learner != "coxph") {warning("Object passed to predict_coxph() may not come from fit_coxph().")}

  stopifnot(requireNamespace("pec", quietly = TRUE))

  survmat <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}


library(survival)
data(veteran)

mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred_cox <- predict_coxph(mod_cox, newdata = veteran[1:5, ], times = 0:100)
pred_cox


summary(mod_cox)
