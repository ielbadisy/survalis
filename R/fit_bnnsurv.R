
fit_bnnsurv <- function(formula, data, ...) {
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  model <- bnnSurvival::bnnSurvival(
    formula = formula,
    data = data,
    ...
  )

  structure(list(
    model = model,
    learner = "bnnsurv",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model")
}


predict_bnnsurv <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "bnnsurv")
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  result <- predict(object$model, newdata)

  all_times <- bnnSurvival::timepoints(result)
  pred_mat <- bnnSurvival::predictions(result)

  # interpolate to requested times
  interp_mat <- t(apply(pred_mat, 1, function(row) {
    approx(x = all_times, y = row, xout = times, method = "linear", rule = 2)$y
  }))

  colnames(interp_mat) <- paste0("t=", times)
  rownames(interp_mat) <- paste0("ID_", seq_len(nrow(newdata)))

  as.data.frame(interp_mat)
}


library(survival)
library(bnnSurvival)

set.seed(123)
mod_bnnsurv <- fit_bnnsurv(Surv(time, status) ~ age + trt + karno + diagtime + prior, data = veteran)

pred <- predict_bnnsurv(mod_bnnsurv, newdata = veteran[1:5, ], times = c(100, 200, 300))
print(round(pred, 3))

