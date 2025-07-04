
fit_cox <- function(formula, data, ...) {
  stopifnot(requireNamespace("survival", quietly = TRUE))
  stopifnot(requireNamespace("pec", quietly = TRUE))

  model <- survival::coxph(formula, data = data, x = TRUE, y = TRUE, ...)

  time_status <- all.vars(formula[[2]])
  structure(list(
    model = model,
    learner = "cox",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model", engine = "cox")
}




predict_cox <- function(object, newdata, times) {
  stopifnot(object$learner == "cox")
  stopifnot(requireNamespace("pec", quietly = TRUE))

  survmat <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  colnames(survmat) <- paste0("t=", times)
  rownames(survmat) <- paste0("ID_", seq_len(nrow(newdata)))
  as.data.frame(survmat)
}


library(survival)
data(veteran)

mod_cox <- fit_cox(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred_cox <- predict_cox(mod_cox, newdata = veteran[1:5, ], times = c(100, 200, 300))
print(round(pred_cox, 3))


