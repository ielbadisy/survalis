
#' Fit an Oblique Random Survival Forest (ORSF) Model
#'
#' Fits an Oblique Random Survival Forest using the \pkg{aorsf} package.
#' This method builds an ensemble of oblique decision trees for survival analysis,
#' where splits are based on linear combinations of features, allowing for
#' improved performance in high-dimensional or correlated feature settings.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing the variables specified in \code{formula}.
#' @param ... Additional arguments passed to \code{aorsf::orsf()}.
#'
#' @details
#' ORSF models extend traditional Random Survival Forests by allowing oblique splits,
#' which can improve prediction accuracy when predictors are correlated.
#' Missing data are omitted by default (\code{na_action = "omit"}).
#'
#' @return An object of class \code{"mlsurv_model"} containing:
#' \itemize{
#'   \item \code{model}: The fitted \pkg{aorsf} ORSF model object.
#'   \item \code{learner}: Character string, always \code{"orsf"}.
#'   \item \code{formula}: The survival formula used.
#'   \item \code{data}: The training dataset.
#'   \item \code{time}: Name of the survival time variable.
#'   \item \code{status}: Name of the event indicator variable.
#' }
#'
#' @references
#' Jaeger BC, Long DL, Long DM, Sims M, Szychowski JM, Min YI, Bandyopadhyay D.
#' Oblique random survival forests. \emph{Annals of Applied Statistics}. 2019;13(3):1847-1883.
#' doi:10.1214/19-AOAS1261
#'
#' @examplesIf requireNamespace("aorsf", quietly = TRUE)
#' mod <- fit_orsf(Surv(time, status) ~ age + karno, data = veteran)
#' summary(mod)
#'
#' @export

fit_orsf <- function(formula, data, ...) {
  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  model <- aorsf::orsf(
    formula = formula,
    data = data,
    na_action = "omit",
    ...
  )

  structure(list(
    model = model,
    learner = "orsf", # Oblique Random Survival Forests

    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "aorsf")
}

#' Predict Survival Probabilities from an ORSF Model
#'
#' Generates survival probability predictions at specified time points
#' from a fitted Oblique Random Survival Forest model.
#'
#' @param object A fitted ORSF model object from \code{\link{fit_orsf}}.
#' @param newdata A data frame of new observations for prediction.
#' @param times A numeric vector of time points at which to estimate survival probabilities.
#' @param ... Additional arguments passed to \code{predict.aorsf()}.
#'
#' @details
#' Predictions are computed using the \code{pred_type = "surv"} option from \pkg{aorsf},
#' returning estimated survival probabilities for each observation at the specified time points.
#'
#' @return A data frame where each row corresponds to an observation in \code{newdata}
#' and each column corresponds to a requested prediction time (\code{"t=<time>"}).
#'
#' @examplesIf requireNamespace("aorsf", quietly = TRUE)
#' mod <- fit_orsf(Surv(time, status) ~ age + karno, data = veteran)
#' pred <- predict_orsf(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#' head(pred)
#'
#' @export

predict_orsf <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "orsf") {
    warning("Object passed to predict_orsf() may not come from fit_orsf().")
  }

  stopifnot(requireNamespace("aorsf", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  survmat <- predict(
    object$model,
    new_data = newdata,
    pred_horizon = times,
    pred_type = "surv",
    na_action = "pass",
    boundary_checks = FALSE,
    ...
  )

  if (is.vector(survmat)) {
    survmat <- matrix(survmat, nrow = 1)
  }

  .finalize_survmat(survmat, times = times)
}


#' Tune Oblique Random Survival Forests (ORSF) via Cross-Validation
#'
#' Performs grid-search hyperparameter tuning for \pkg{aorsf} models using
#' cross-validation and selects the best configuration by the primary metric.
#' Optionally refits the best model on the full dataset.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing all variables referenced in \code{formula}.
#' @param times Numeric vector of evaluation time points (same scale as the training
#'   survival time). Must be positive and finite.
#' @param param_grid A \code{data.frame} (e.g., from \code{expand.grid()}) of candidate
#'   hyperparameters. Typical columns include:
#'   \itemize{
#'     \item \code{n_tree}: number of trees (e.g., \code{c(100, 300)})
#'     \item \code{mtry}: number of variables tried at each split
#'     \item \code{min_events}: minimum number of events required to attempt a split
#'   }
#' @param metrics Character vector of metrics to evaluate/optimize
#'   (e.g., \code{"cindex"}, \code{"ibs"}). The first entry is used as the primary
#'   selection metric.
#' @param folds Integer; number of cross-validation folds. Default is \code{5}.
#' @param seed Integer random seed for reproducibility. Default is \code{123}.
#' @param refit_best Logical; if \code{TRUE}, refits the best configuration on the
#'   full data and returns it. Default is \code{FALSE}.
#' @param ... Additional arguments forwarded to the underlying engine via
#'   \code{\link{fit_orsf}}.
#'
#' @return If \code{refit_best = FALSE}, a \code{data.frame} (class \code{"tuned_surv"})
#'   with one row per grid combination and columns for hyperparameters and metrics,
#'   ordered by the first metric. If \code{refit_best = TRUE}, a fitted
#'   \code{mlsurv_model} from \code{\link{fit_orsf}} using the best settings.
#'
#' @details
#' Internally calls \code{cv_survlearner()} with \code{\link{fit_orsf}} /
#' \code{\link{predict_orsf}} so tuning uses the exact same code paths as production.
#' The \code{min_events} column of \code{param_grid} is passed to the engine as
#' \code{n_split} (minimum events for a candidate split).
#'
#' @seealso \code{\link{fit_orsf}}, \code{\link{predict_orsf}}, \pkg{aorsf}
#'
#' @examplesIf requireNamespace("aorsf", quietly = TRUE)
#' \donttest{
#'   # Cross-validated tuning (small grid to keep runtime light)
#'   res_orsf <- tune_orsf(
#'     formula    = Surv(time, status) ~ age + karno,
#'     data       = veteran,
#'     times      = c(100, 200, 300),
#'     param_grid = expand.grid(
#'       n_tree     = c(100, 200),
#'       mtry       = c(1, 2),
#'       min_events = c(5, 10)
#'     ),
#'     metrics    = c("cindex", "ibs"),
#'     folds      = 2
#'   )
#'   print(res_orsf)
#'
#'   # Refit the best model on the full data
#'   mod_orsf_best <- tune_orsf(
#'     formula    = Surv(time, status) ~ age + karno,
#'     data       = veteran,
#'     times      = c(100, 200, 300),
#'     param_grid = expand.grid(
#'       n_tree     = c(100, 200),
#'       mtry       = c(1, 2),
#'       min_events = c(5, 10)
#'     ),
#'     metrics    = c("cindex", "ibs"),
#'     folds      = 2,
#'     refit_best = TRUE
#'   )
#'
#'   summary(mod_orsf_best)
#' }
#' @export

tune_orsf <- function(formula, data, times,
                      param_grid = expand.grid(
                        n_tree = c(100, 300),
                        mtry = c(2, 3),
                        min_events = c(5, 10)
                      ),
                      metrics = c("cindex", "ibs"),
                      folds = 5,
                      seed = 123,
                      ncores = 1,
                      refit_best = FALSE,
                      ...) {

  results <- purrr::pmap_dfr(param_grid, function(n_tree, mtry, min_events) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_orsf,
      pred_fun = predict_orsf,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ncores = ncores,
      n_tree = n_tree,
      mtry = mtry,
      n_split = min_events,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      n_tree = n_tree,
      mtry = mtry,
      min_events = min_events
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
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
    best_model <- fit_orsf(
      formula = formula,
      data = data,
      n_tree = best_row$n_tree,
      mtry = best_row$mtry,
      n_split = best_row$min_events,
      ...
    )
    return(best_model)
  }
}
