
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


#--------- add tuner

#' Tune glmnet Cox model over alpha
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data.frame
#' @param times Vector of time points for evaluation
#' @param param_grid A data.frame with values of alpha (default = 0 to 1)
#' @param metrics Vector of metrics (e.g., "cindex", "ibs")
#' @param folds Number of CV folds
#' @param seed Random seed for reproducibility
#' @param ... Additional arguments passed to fit_glmnet
#' @return A tibble with average metric scores per alpha
#' @export
tune_glmnet <- function(formula, data, times,
                        param_grid = expand.grid(alpha = seq(0, 1, by = 0.25)),
                        metrics = c("cindex", "ibs"),
                        folds = 5,
                        seed = 123, ...) {

  purrr::pmap_dfr(param_grid, function(alpha) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_glmnet,
      pred_fun = predict_glmnet,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      alpha = alpha,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(alpha = alpha) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  }) |> arrange(dplyr::across(any_of(metrics[1])))
}





res_glmnet <- tune_glmnet(
  formula = Surv(time, status) ~ age + karno + celltype,
  param_grid = expand.grid(alpha = seq(0, 5, by = 0.25)),
  data = veteran,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 5
)
