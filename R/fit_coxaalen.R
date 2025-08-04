
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
    learner = "coxaalen",
    formula = new_formula,
    data = data
  ), class = "mlsurv_model", engine = "timreg")
}




predict_coxaalen <- function(object, newdata, times = NULL, ...) {

  if (!is.null(object$learner) && object$learner != "coxaalen") {
    warning("Object passed to predict_coxaalen() may not come from fit_coxaalen().")
  }

  stopifnot(requireNamespace("timereg", quietly = TRUE))

  pred <- predict(object$model, newdata = newdata, times = times, se = FALSE, uniform = FALSE, ...)

  survmat <- pred$S0
  if (is.null(survmat)) stop("No survival probabilities returned.")

  if (!is.null(times)) {
    keep_cols <- which(pred$time %in% times)
    if (length(keep_cols) == 0) stop("None of the requested times matched prediction grid.")
    survmat <- survmat[, keep_cols, drop = FALSE]
    colnames(survmat) <- paste0("t=", pred$time[keep_cols])
  } else {
    colnames(survmat) <- paste0("t=", pred$time)
  }

  as.data.frame(survmat)
}


library(timereg)
data(sTRACE)

mod_caalen <- fit_coxaalen(
  Surv(time, status == 9) ~ age + sex + vf + chf + diabetes,
  data = sTRACE
)

pred_surv <- predict_coxaalen(mod_caalen, newdata = sTRACE[1:5, ], times = 0:10)

pred_surv
