
fit_coxaalen <- function(formula, data, max.time = NULL, n.sim = 100, ...) {
  stopifnot(requireNamespace("timereg", quietly = TRUE))

  terms_rhs <- attr(terms(formula), "term.labels")
  new_rhs <- paste0("prop(", terms_rhs, ")", collapse = " + ")
  new_formula <- as.formula(paste(deparse(formula[[2]]), "~", new_rhs))

  model <- timereg::cox.aalen(
    formula = new_formula,
    data = data,
    max.time = max.time,
    n.sim = n.sim,
    robust = 1,
    ...
  )

  structure(list(
    model = model,
    learner = "cox.aalen",
    formula = new_formula,
    data = data
  ), class = "mlsurv_model", engine = "cox.aalen")
}




predict_coxaalen <- function(object, newdata, times = NULL, ...) {
  stopifnot(object$learner == "cox.aalen")
  stopifnot(requireNamespace("timereg", quietly = TRUE))

  pred <- predict(object$model, newdata = newdata, times = times, se = FALSE, uniform = FALSE, ...)

  surv_probs <- pred$S0
  if (is.null(surv_probs)) stop("No survival probabilities returned.")

  if (!is.null(times)) {
    keep_cols <- which(round(pred$time, 4) %in% round(times, 4))
    if (length(keep_cols) == 0) stop("None of the requested times matched prediction grid.")
    surv_probs <- surv_probs[, keep_cols, drop = FALSE]
    colnames(surv_probs) <- paste0("t=", round(pred$time[keep_cols], 3))
  } else {
    colnames(surv_probs) <- paste0("t=", round(pred$time, 3))
  }

  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))
  return(as.data.frame(surv_probs))
}


library(timereg)
data(sTRACE)

mod_caalen <- fit_coxaalen(
  Surv(time, status == 9) ~ age + sex + vf + chf + diabetes,
  data = sTRACE
)

pred_surv <- predict_coxaalen(mod_caalen, newdata = sTRACE[1:5, ], times = 4:7)

print(round(pred_surv, 3))
