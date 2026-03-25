
#' Fit a Penalized Cox Proportional Hazards Model (glmnet)
#'
#' Fits a Cox proportional hazards model with elastic net regularization
#' using \pkg{glmnet}. The function wraps \code{\link[glmnet]{cv.glmnet}} to
#' select the optimal penalty parameter \eqn{\lambda} by cross-validation.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing all variables referenced in \code{formula}.
#' @param alpha Numeric value in \code{[0, 1]} specifying the elastic net
#'   mixing parameter:
#'   \itemize{
#'     \item \code{1} for LASSO penalty
#'     \item \code{0} for ridge penalty
#'     \item between 0 and 1 for elastic net
#'   }
#' @param ... Additional arguments passed to \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return An object of class \code{"mlsurv_model"} containing:
#' \itemize{
#'   \item \code{model} -- the fitted \code{cv.glmnet} object
#'   \item \code{learner} -- the string \code{"glmnet"}
#'   \item \code{formula}, \code{data}, \code{time}, and \code{status} metadata
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link{predict_glmnet}}
#'
#' @examples
#'
#' mod_glmnet <- fit_glmnet(
#'   Surv(time, status) ~ age + karno + celltype,
#'   data = veteran
#' )
#'
#' summary(mod_glmnet$model)
#'
#' @export

fit_glmnet <- function(formula, data, alpha = 1, ...) {
  stopifnot(requireNamespace("glmnet", quietly = TRUE))

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



#' Predict Survival Probabilities from a Penalized Cox Model (glmnet)
#'
#' Predicts survival probabilities at specified time points from a fitted
#' penalized Cox proportional hazards model produced by \code{\link{fit_glmnet}}.
#'
#' @param object A fitted \code{"mlsurv_model"} object from \code{\link{fit_glmnet}}.
#' @param newdata A \code{data.frame} of new observations for prediction.
#' @param times Numeric vector of time points at which to estimate survival
#'   probabilities. Must be in the same scale as the training time variable.
#' @param ... Not used.
#'
#' @return A \code{data.frame} with one row per observation in \code{newdata}
#'   and one column per requested time point (\code{"t=<time>"}).
#'
#' @details
#' Predictions are computed by:
#' \enumerate{
#'   \item Computing the linear predictors for training and new data at
#'     \code{s = "lambda.min"}.
#'   \item Fitting a Cox PH model on the training linear predictors to
#'     estimate the baseline cumulative hazard.
#'   \item Interpolating the baseline hazard at each \code{times} and
#'     transforming via \eqn{S(t \mid x) = \exp(-H_0(t) \exp(\eta))}.
#' }
#'
#' @seealso \code{\link{fit_glmnet}}
#'
#' @examples
#'
#' mod_glmnet <- fit_glmnet(
#'   Surv(time, status) ~ age + karno + celltype,
#'   data = veteran
#' )
#'
#' predict_glmnet(
#'   mod_glmnet,
#'   newdata = veteran[1:5, ],
#'   times = c(100, 200, 300)
#' )
#'
#' @export

predict_glmnet <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "glmnet") {
    warning("Object passed to predict_glmnet() may not come from fit_glmnet().")
  }

  stopifnot(requireNamespace("glmnet", quietly = TRUE))

  x_new   <- model.matrix(object$formula, newdata)[, -1, drop = FALSE]
  x_train <- model.matrix(object$formula, object$data)[, -1, drop = FALSE]
  y_train <- survival::Surv(object$data[[object$time]], object$data[[object$status]])

  # predict linear predictors
  lp_train <- predict(object$model, newx = x_train, s = "lambda.min", type = "link")[, 1]
  lp_new   <- predict(object$model, newx = x_new,   s = "lambda.min", type = "link")[, 1]
  .predict_cox_from_lp(lp_train = lp_train, lp_new = lp_new, y_train = y_train, times = times)
}

#' Tune Penalized Cox Proportional Hazards Model via Cross-Validation
#'
#' Performs hyperparameter tuning for penalized Cox models (\pkg{glmnet}) over
#' a grid of \code{alpha} values using cross-validation on one or more metrics.
#' Optionally refits the best model on the full dataset.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing all variables referenced in \code{formula}.
#' @param times Numeric vector of time points at which to evaluate performance.
#' @param param_grid A \code{data.frame} or list specifying candidate \code{alpha} values.
#'   Typically created with \code{\link[base]{expand.grid}}.
#' @param metrics Character vector of performance metrics to compute. The first
#'   entry is used as the primary selection metric.
#' @param folds Integer; number of cross-validation folds. Default is \code{5}.
#' @param seed Integer random seed for reproducibility. Default is \code{123}.
#' @param refit_best Logical; if \code{TRUE}, fits the best configuration on the
#'   full dataset and returns the fitted model. If \code{FALSE} (default), returns
#'   a results table.
#' @param ... Additional arguments passed to \code{\link{fit_glmnet}}.
#'
#' @return
#' If \code{refit_best = FALSE}, returns a \code{data.frame} (class \code{"tuned_surv"})
#' with hyperparameters and metric values for each grid combination, sorted by
#' the primary metric.  
#' If \code{refit_best = TRUE}, returns a fitted \code{"mlsurv_model"} from
#' \code{\link{fit_glmnet}}.
#'
#' @seealso \code{\link{fit_glmnet}}, \code{\link{predict_glmnet}}
#'
#' @examples
#'
#' param_grid <- expand.grid(alpha = seq(0, 1, by = 0.25))
#'
#' # Evaluate grid without refitting
#' res_glmnet <- tune_glmnet(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   param_grid = param_grid,
#'   metrics = c("cindex", "ibs"),
#'   folds = 3
#' )
#' print(res_glmnet)
#'
#' # Refit best model directly
#' mod_glmnet_best <- tune_glmnet(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   param_grid = param_grid,
#'   refit_best = TRUE
#' )
#'
#' summary(mod_glmnet_best)
#'
#' @export

tune_glmnet <- function(formula, data, times,
                        param_grid = c(alpha = seq(0, 1, by = 0.25)),
                        metrics = c("cindex", "ibs"),
                        folds = 5,
                        seed = 123,
                        ncores = 1,
                        refit_best = FALSE,
                        ...) {

  grid_df <- tidyr::crossing(!!!param_grid)

  results <- purrr::pmap_dfr(grid_df, function(alpha) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_glmnet,
      pred_fun = predict_glmnet,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ncores = ncores,
      alpha = alpha,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(alpha = alpha) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(
          summary[, c("metric", "mean")],
          names_from = metric,
          values_from = mean
        )
      )
  })

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))
  } else {
    results <- results |> dplyr::arrange(!!rlang::sym(metrics[1]))
  }

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_glmnet(
      formula = formula,
      data = data,
      alpha = best_row$alpha,
      ...
    )
    return(best_model)
  }
}
