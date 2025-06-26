
fit_bnnsurv <- function(formula, data,
                        k = 2,
                        num_base_learners = 50,
                        num_features_per_base_learner = NULL,
                        metric = "mahalanobis",
                        weighting_function = function(x) x * 0 + 1,
                        replace = TRUE,
                        sample_fraction = NULL,
                        ...) {
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  model <- bnnSurvival::bnnSurvival(
    formula = formula,
    data = data,
    k = k,
    num_base_learners = num_base_learners,
    num_features_per_base_learner = num_features_per_base_learner,
    metric = metric,
    weighting_function = weighting_function,
    replace = replace,
    sample_fraction = sample_fraction,
    ...
  )

  structure(list(
    model = model,
    learner = "bnnsurv",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "bnnsurv")
}


predict_bnnsurv <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "bnnsurv")
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  # predict full survival curves
  pred <- predict(object$model, newdata)

  # extract prediction grid
  pred_times <- bnnSurvival::timepoints(pred)
  pred_matrix <- bnnSurvival::predictions(pred)

  surv_probs <- t(apply(pred_matrix, 1, function(row) {
    approx(x = pred_times, y = row, xout = times, method = "linear", rule = 2)$y
  }))

  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))

  return(as.data.frame(surv_probs))
}


library(survival)
library(bnnSurvival)

mod_bnnsurv <- fit_bnnsurv(Surv(time, status) ~ age + karno + diagtime + prior, data = veteran)

pred_bnnsurv <- predict_bnnsurv(mod_bnnsurv, newdata = veteran[1:5, ], times = c(100, 200, 300))
print(round(pred_bnnsurv, 3))
