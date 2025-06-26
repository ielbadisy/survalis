fit_ranger <- function(formula, data, ...) {
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  model <- ranger::ranger(
    formula = formula,
    data = data,
    ...,
    respect.unordered.factors = "order",
    case.weights = NULL
  )

  structure(list(
    model = model,
    learner = "ranger",
    formula = formula,
    data = data
  ), class = "mlsurv_model", engine = "ranger")
}


predict_ranger <- function(model, newdata, times) {
  requireNamespace("ranger")

  pred_obj <- predict(model$model, data = newdata)
  pred <- pred_obj$survival
  model_times <- model$model$unique.death.times

  if (is.null(pred)) stop("Prediction returned NULL survival matrix.")
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)  # ensure matrix for single row

  # preallocate output
  surv_mat <- matrix(NA_real_, nrow = nrow(pred), ncol = length(times))
  for (i in seq_len(nrow(pred))) {
    y <- pred[i, ]
    surv_mat[i, ] <- stats::approx(
      x = model_times, y = y, xout = times,
      method = "linear", rule = 2, ties = "ordered"
    )$y
  }

  colnames(surv_mat) <- paste0("t=", times)
  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(surv_mat)))
  return(as.data.frame(surv_mat))
}


# example test
library(survival)
data(veteran, package = "survival")

mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred <- predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))

pred
