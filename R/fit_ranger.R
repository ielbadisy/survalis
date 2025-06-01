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
  # check that the model is from the expected 'ranger' learner
  stopifnot(object$learner == "ranger")
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  # predict survival curve using ranger
  pred <- predict(object$model, data = newdata, type = "response")

  surv_mat <- pred$survival            # matrix of survival probabilities over time
  model_times <- pred$unique.death.times  # corresponding time points

  if (is.null(times)) {
    times <- model_times               # if not specified, return all default prediction times
  }

  # interpolate survival probabilities at requested times using linear interpolation
  surv_probs <- apply(surv_mat, 1, function(row_surv) {
    stats::approx(
      x = model_times,
      y = row_surv,
      xout = times,
      method = "linear",
      rule = 2                     # rule = 2 ensures extrapolation is allowed if time > max
    )$y
  })

  # --- FIX: Ensure consistent output shape for 1 time point ---
  # when 'times' has length 1, 'apply' returns a vector instead of a matrix.
  # wrap it into a matrix to maintain consistent dimensions (rows = observations)
  if (length(times) == 1) {
    surv_probs <- matrix(surv_probs, ncol = 1)
  } else {
    surv_probs <- t(surv_probs)   # otherwise, transpose the result to get one row per obs
  }

  # add readable column and row names
  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))

  # return as a tidy data.frame
  as.data.frame(surv_probs)
}



# example test
library(survival)
data(veteran, package = "survival")

mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred <- predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))

pred
