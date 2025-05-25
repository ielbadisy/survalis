fit_survivalsvm <- function(formula, data,
                            type = "regression",
                            gamma.mu = 0.1,
                            opt.meth = "quadprog",
                            kernel = "add_kernel",
                            diff.meth = NULL,
                            ...) {
  stopifnot(requireNamespace("survivalsvm", quietly = TRUE))

  model <- survivalsvm::survivalsvm(
    formula = formula,
    data = data,
    type = type,
    gamma.mu = gamma.mu,
    opt.meth = opt.meth,
    kernel = kernel,
    diff.meth = diff.meth,
    ...
  )

  structure(list(
    model = model,
    learner = "survivalsvm",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "survivalsvm")
}





predict_survivalsvm <- function(object, newdata, times, dist = "exp", shape = 1) {
  stopifnot(object$learner == "survivalsvm")
  stopifnot(requireNamespace("survivalsvm", quietly = TRUE))

  newdata <- as.data.frame(newdata)
  pred <- predict(object$model, newdata = newdata)
  predicted_times <- as.numeric(pred$predicted)

  dist <- match.arg(dist, choices = c("exp", "weibull"))

  probs <- switch(dist,
                  "exp" = sapply(times, function(t) exp(-t / predicted_times)),
                  "weibull" = sapply(times, function(t) exp(- (t / predicted_times)^shape))
  )

  colnames(probs) <- paste0("t=", times)
  rownames(probs) <- paste0("ID_", seq_along(predicted_times))
  return(as.data.frame(probs))
}

library(survival)
library(survivalsvm)

mod_svm <- fit_survivalsvm(Surv(time, status) ~ age + celltype + karno,
                           data = veteran,
                           type = "regression",
                           gamma.mu = 0.1,
                           kernel = "lin_kernel")

times <- c(100, 300, 500)
predict_survivalsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "exp")
predict_survivalsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "weibull", shape = 1.5)

