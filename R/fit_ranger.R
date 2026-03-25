
#' Fit a Survival Random Forest Model Using ranger
#'
#' This function fits a survival random forest model using the
#' \pkg{ranger} package and returns an object compatible with the
#' `mlsurv_model` class.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param ... Additional arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @details
#' This function wraps \code{\link[ranger]{ranger}} for survival analysis and
#' stores the result in a standardized `mlsurv_model` object.
#'
#' @return An object of class \code{"mlsurv_model"} containing:
#' \itemize{
#'   \item \code{model} - the fitted \code{ranger} model
#'   \item \code{learner} - character string "ranger"
#'   \item \code{formula} - the model formula
#'   \item \code{data} - training data used to fit the model
#' }
#'
#' @seealso \code{\link{predict_ranger}}, \code{\link{tune_ranger}}, \code{\link[ranger]{ranger}}
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' summary(mod)
#'
#' @export

fit_ranger <- function(formula, data, ...) {
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  model <- ranger::ranger(
    formula = formula,
    data = data,
    ...,
    respect.unordered.factors = "order",
    case.weights = NULL
  )

  structure(
    list(
      model = model,
      learner = "ranger",
      formula = formula,
      data = data
    ),
    class = "mlsurv_model",
    engine = "ranger"
  )
}


#' Predict Survival Probabilities from a ranger Model
#'
#' Generates predicted survival probabilities for given time points
#' from a model fitted with \code{\link{fit_ranger}}.
#'
#' @param object An \code{"mlsurv_model"} object returned by \code{\link{fit_ranger}}.
#' @param newdata A \code{data.frame} containing the same predictors as the training data.
#' @param times Numeric vector of time points at which to estimate survival probabilities.
#'
#' @details
#' Predictions are obtained from \code{\link[ranger]{predict.ranger}} and
#' survival curves are interpolated to match the requested \code{times}.
#'
#' @return A \code{data.frame} with one row per observation in \code{newdata}
#' and one column per time point (columns named \code{"t=<time>"}).
#'
#' @seealso \code{\link{fit_ranger}}, \code{\link{tune_ranger}}
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export

predict_ranger <- function(object, newdata, times) {



  if (!is.null(object$learner) && object$learner != "ranger") {
    warning("Object passed to predict_ranger() may not come from fit_ranger().")
  }

  stopifnot(requireNamespace("ranger", quietly = TRUE))

  pred_obj <- predict(object$model, data = newdata)
  pred <- pred_obj$survival
  model_times <- object$model$unique.death.times

  if (is.null(pred)) stop("Prediction returned NULL survival matrix.")
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)

  survmat <- matrix(NA_real_, nrow = nrow(pred), ncol = length(times))
  for (i in seq_len(nrow(pred))) {
    y <- pred[i, ]
    survmat[i, ] <- stats::approx(
      x = model_times, y = y, xout = times,
      method = "linear", rule = 2, ties = "ordered"
    )$y
  }

  .finalize_survmat(survmat, times = times)
}



#' Hyperparameter Tuning for ranger Survival Models
#'
#' Performs grid search tuning of a survival random forest model using
#' \pkg{ranger} over a set of hyperparameter combinations.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing the variables in the model.
#' @param times Numeric vector of time points at which to evaluate performance.
#' @param param_grid A \code{data.frame} of hyperparameter combinations. Must contain
#'   columns \code{num.trees}, \code{mtry}, and \code{min.node.size}.
#' @param metrics Character vector of evaluation metrics (e.g., \code{"cindex"}, \code{"ibs"}).
#' @param folds Number of cross-validation folds.
#' @param seed Random seed for reproducibility.
#' @param refit_best Logical; if \code{TRUE}, refits the model using the best parameters.
#' @param ... Additional arguments passed to \code{\link{fit_ranger}}.
#'
#' @details
#' Uses \code{\link{cv_survlearner}} to perform cross-validation for each
#' parameter combination. If \code{refit_best = TRUE}, the function returns the
#' best-fitting model; otherwise, it returns a tuning results table.
#'
#' @return
#' If \code{refit_best = FALSE}, returns a \code{data.frame} of tuning results
#' sorted by the primary metric. If \code{refit_best = TRUE}, returns an
#' \code{"mlsurv_model"} object with the best parameters and an attribute
#' \code{"tuning_results"}.
#'
#' @seealso \code{\link{fit_ranger}}, \code{\link{predict_ranger}}
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' mod_ranger_best <- tune_ranger(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   param_grid = expand.grid(
#'     num.trees = c(100, 300),
#'     mtry = c(1, 2),
#'     min.node.size = c(3, 5)
#'   ),
#'   metrics = c("cindex", "ibs"),
#'   folds = 3,
#'   refit_best = TRUE
#' )
#' summary(mod_ranger_best)
#'
#' @export

tune_ranger <- function(formula, data, times,
                        param_grid = expand.grid(
                          num.trees = c(100, 300),
                          mtry = c(1, 2, 3),
                          min.node.size = c(3, 5)
                        ),
                        metrics = c("cindex", "ibs"),
                        folds = 5, seed = 123,
                        ncores = 1,
                        refit_best = FALSE,
                        ...) {

  # cv loop
  res <- purrr::pmap_dfr(param_grid, function(num.trees, mtry, min.node.size) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_ranger,
      pred_fun = predict_ranger,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ncores = ncores,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(num.trees, mtry, min.node.size) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    res <- dplyr::arrange(res, dplyr::desc(.data[[metrics[1]]]))
  } else {
    res <- dplyr::arrange(res, .data[[metrics[1]]])
  }

  if (refit_best) {
    best <- res[1, ]
    mod <- fit_ranger(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      num.trees = best$num.trees,
      mtry = best$mtry,
      min.node.size = best$min.node.size,
      ...
    )
    attr(mod, "tuning_results") <- res
    class(mod) <- c("tune_surv", class(mod))  # this line enables summary()
    return(mod)
  }

  return(res)
}
