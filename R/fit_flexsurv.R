fit_flexsurv <- function(formula, data, dist = "weibull", ...) {
  stopifnot(requireNamespace("flexsurv", quietly = TRUE))

  model <- flexsurv::flexsurvreg(formula = formula, data = data, dist = dist, ...)

  structure(list(
    model = model,
    learner = "flexsurv",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "flexsurv")
}

predict_flexsurv <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "flexsurv")
  stopifnot(requireNamespace("flexsurv", quietly = TRUE))

  pred <- predict(
    object$model,
    newdata = newdata,
    type = "survival",
    times = times
  )

  # Unpack .pred (a list-column of tibbles)
  .pred <- pred$.pred
  n <- nrow(newdata)
  k <- length(times)
  mat <- matrix(NA_real_, nrow = n, ncol = k)
  rownames(mat) <- paste0("ID_", seq_len(n))
  colnames(mat) <- paste0("t=", times)

  for (i in seq_len(n)) {
    mat[i, ] <- .pred[[i]]$.pred_survival
  }

  return(as.data.frame(mat))
}

library(survival)
library(flexsurv)

data(veteran, package = "survival")

# Fit the model
mod_flex <- fit_flexsurv(Surv(time, status) ~ age + celltype + karno,
                         data = veteran,
                         dist = "weibull")

# Predict survival probabilities at given time points
times <- c(100, 200, 300)
predicted_surv <- predict_flexsurv(mod_flex, newdata = veteran[1:5, ], times = times)

# View result
print(round(predicted_surv, 3))


