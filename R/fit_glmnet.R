
fit_glmnet <- function(formula, data, alpha = 1, ...) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE))

  # Extract survival outcome and predictor matrix
  x <- model.matrix(formula, data = data)[, -1, drop = FALSE]
  y <- survival::Surv(data[[all.vars(formula)[1]]], data[[all.vars(formula)[2]]])

  model <- glmnet::cv.glmnet(x, y, family = "cox", alpha = alpha, ...)

  structure(list(
    model = model,
    learner = "glmnet",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2]
  ), class = "mlsurv_model", engine = "glmnet")
}


# predict survival probabilities using a glmnet cox model

predict_glmnet <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "glmnet")
  stopifnot(requireNamespace("glmnet", quietly = TRUE))

  x_new   <- model.matrix(object$formula, newdata)[, -1, drop = FALSE]
  x_train <- model.matrix(object$formula, object$data)[, -1, drop = FALSE]
  y_train <- survival::Surv(object$data[[object$time]], object$data[[object$status]])

  # predict linear predictors
  lp_train <- predict(object$model, newx = x_train, s = "lambda.min", type = "link")[, 1]
  lp_new   <- predict(object$model, newx = x_new,   s = "lambda.min", type = "link")[, 1]

  # estimate baseline hazard
  base_haz <- survival::basehaz(survival::coxph(y_train ~ lp_train), centered = FALSE)

  # interpolate cumulative hazard at requested times and compute survival
  surv_probs <- sapply(times, function(t) {
    H0_t <- approx(base_haz$time, base_haz$hazard, xout = t, rule = 2)$y
    exp(-H0_t * exp(lp_new))
  })

  if (is.vector(surv_probs)) {
    surv_probs <- matrix(surv_probs, ncol = length(times))
  }

  colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))

  return(as.data.frame(surv_probs))
}


library(survival)
data(veteran)

mod_glmnet <- fit_glmnet(Surv(time, status) ~ age + karno + celltype, data = veteran)

pred_probs <- predict_glmnet(mod_glmnet, newdata = veteran[1:5, ], times = c(100, 200, 300))

print(round(pred_probs, 3))

