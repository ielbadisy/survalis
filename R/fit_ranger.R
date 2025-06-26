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

  # get full survival curve estimates
  pred <- predict(model$model, data = newdata)$survival
  model_times <- model$model$unique.death.times

  # handle missing survival matrix (can happen with 1-row input)
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)

  # interpolate survival probabilities at desired time points
  surv_mat <- t(apply(pred, 1, function(row_surv) {
    stats::approx(x = model_times, y = row_surv, xout = times,
                  method = "linear", rule = 2)$y
  }))

  colnames(surv_mat) <- paste0("t=", times)
  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(surv_mat)))
  return(surv_mat)
}

# example test
library(survival)
data(veteran, package = "survival")

mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred <- predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))

pred
