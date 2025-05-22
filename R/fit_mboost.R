
fit_mboost <- function(formula, data, weights = NULL, mstop = 100, nu = 0.1,
                       minsplit = 10, minbucket = 4, maxdepth = 2, ...) {
  stopifnot(requireNamespace("mboost", quietly = TRUE))
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  data <- as.data.frame(data)
  family <- mboost::CoxPH()

  control <- mboost::boost_control(mstop = mstop, nu = nu, ...)
  tree_control <- partykit::ctree_control(minsplit = minsplit, minbucket = minbucket, maxdepth = maxdepth)

  # conditionally include weights
  model <- if (!is.null(weights)) {
    mboost::blackboost(
      formula = formula,
      data = data,
      weights = weights,
      family = family,
      control = control,
      tree_controls = tree_control
    )
  } else {
    mboost::blackboost(
      formula = formula,
      data = data,
      family = family,
      control = control,
      tree_controls = tree_control
    )
  }

  structure(
    list(
      model = model,
      learner = "mboost",
      formula = formula,
      data = data,
      time = all.vars(formula)[[2]],
      status = all.vars(formula)[[3]]
    ),
    class = "mlsurv_model",
    engine = "mboost"
  )
}



survival_curve_to_prob <- function(eval_time, event_times, survival_prob) {

    if (event_times[1] != -Inf) {
    event_times <- c(-Inf, event_times)
    survival_prob <- rbind(1, survival_prob)
  }
  if (event_times[length(event_times)] != Inf) {
    event_times <- c(event_times, Inf)
    survival_prob <- rbind(survival_prob, 0)
  }

  index <- findInterval(eval_time, event_times)
  survival_prob[index, , drop = FALSE]
}


predict_mboost <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "mboost")
  stopifnot(requireNamespace("mboost", quietly = TRUE))

  surv_fit <- mboost::survFit(object$model, newdata = newdata)

  # ensure `times` is within the range of surv_fit$time
  max_time <- max(surv_fit$time)
  if (any(times > max_time)) {
    warning("Some prediction times exceed the range of survFit$time. Truncating.")
    times <- times[times <= max_time]
  }

  if (length(times) == 0) stop("No valid `times` for prediction after filtering.")

  # transpose so individuals are rows
  surv_probs <- t(survival_curve_to_prob(
    eval_time = times,
    event_times = surv_fit$time,
    survival_prob = surv_fit$surv
  ))

  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- rownames(newdata)

  return(as.data.frame(surv_probs))
}




library(survival)
data("veteran", package = "survival")

mod_mboost <- fit_mboost(Surv(time, status) ~ age + karno + celltype, data = veteran)

predict_mboost(mod_mboost, newdata = veteran[1:5, ], times = c(100, 200, 300))




# Load data and required packages
library(survival)
data("veteran", package = "survival")

# Fit model with CV
cv_results_mboost <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  fit_fun = fit_mboost,
  pred_fun = predict_mboost,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42
)

# Summary and visualization
print(cv_results_mboost)
cv_summary_
