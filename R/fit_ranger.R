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


predict_ranger <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "ranger")
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  pred <- predict(object$model, data = newdata, type = "response")

  surv_mat <- pred$survival
  model_times <- pred$unique.death.times

  if (is.null(times)) {
    times <- model_times
  }

  # Interpolate survival probabilities for requested times
  surv_probs <- apply(surv_mat, 1, function(row_surv) {
    stats::approx(
      x = model_times,
      y = row_surv,
      xout = times,
      method = "linear",
      rule = 2
    )$y
  })

  surv_probs <- t(surv_probs)  # rows = observations
  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))

  as.data.frame(surv_probs)
}

# example test
library(survival)
data(veteran, package = "survival")

mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred <- predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))

pred
