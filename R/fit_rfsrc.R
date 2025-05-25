fit_rfsrc <- function(formula, data,
                      ntree = 500,
                      mtry = NULL,
                      nodesize = 15,
                      ...) {
  stopifnot(requireNamespace("randomForestSRC", quietly = TRUE))

  model <- randomForestSRC::rfsrc(
    formula = formula,
    data = data,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    ...
  )

  structure(list(
    model = model,
    learner = "rfsrc",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "rfsrc")
}


predict_rfsrc <- function(object, newdata, times = NULL, ...) {
  stopifnot(object$learner == "rfsrc")
  stopifnot(requireNamespace("randomForestSRC", quietly = TRUE))

  pred <- randomForestSRC::predict.rfsrc(object$model, newdata = newdata)
  surv_mat <- pred$survival
  time_points <- pred$time.interest

  # interpolation to align predictions to `times`
  if (!is.null(times)) {
    idx <- sapply(times, function(t) which.min(abs(time_points - t)))
    surv_mat <- surv_mat[, idx, drop = FALSE]
    colnames(surv_mat) <- paste0("t=", times)
  } else {
    colnames(surv_mat) <- paste0("t=", time_points)
  }

  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(newdata)))
  return(as.data.frame(surv_mat))
}


library(survival)
library(randomForestSRC)

data(veteran, package = "survival")

mod_rsf <- fit_rfsrc(Surv(time, status) ~ age + celltype + karno, data = veteran, ntree = 200)

times <- c(100, 200, 300)
pred_probs <- predict_rfsrc(mod_rsf, newdata = veteran[1:5, ], times = times)

print(round(pred_probs, 3))
