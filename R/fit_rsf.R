#' Fit a Random Survival Forest (RSF) Model
#'
#' Fits a Random Survival Forest for right-censored time-to-event data using
#' \pkg{randomForestSRC} and returns a standardized `mlsurv_model` compatible
#' with the `survalis` framework.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#' @param data A `data.frame` containing the variables referenced in `formula`.
#' @param ntree Integer; number of trees to grow (default: 500).
#' @param mtry Integer or `NULL`; number of variables randomly selected at each split.
#'   If `NULL`, the engine's default is used.
#' @param nodesize Integer; minimum terminal node size (default: 15).
#' @param ... Additional arguments forwarded to \code{\link[randomForestSRC]{rfsrc}}.
#'
#' @return An object of class `mlsurv_model`, a named list with elements:
#' \describe{
#'   \item{model}{The fitted \code{randomForestSRC::rfsrc} object.}
#'   \item{learner}{Character scalar identifying the learner (`"rsf"`).}
#'   \item{engine}{Character scalar naming the engine (`"randomForestSRC"`).}
#'   \item{formula}{The original survival formula.}
#'   \item{data}{The training dataset (or a minimal subset needed for prediction).}
#'   \item{time}{Name of the survival time variable.}
#'   \item{status}{Name of the event indicator (1 = event, 0 = censored).}
#' }
#'
#' @details
#' RSF extends random forests to survival data by growing an ensemble of survival
#' trees on bootstrap samples and aggregating survival functions across trees.
#'
#' @references
#' Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS (2008).
#' Random survival forests. \emph{Annals of Applied Statistics}, 2(3):841-860.
#'
#' @examples
#' mod_rsf <- fit_rsf(Surv(time, status) ~ age + celltype + karno,
#'                    data = veteran, ntree = 200)
#'
#' times <- c(100, 200, 300)
#' pred_probs <- predict_rsf(mod_rsf, newdata = veteran[1:5, ], times = times)
#' print(round(pred_probs, 3))
#'
#' @export

fit_rsf <- function(formula, data,
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
    learner = "rsf",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "randomForestSRC")
}

#' Predict Survival Probabilities from an RSF Model
#'
#' Generates survival probability predictions for new data using a model
#' fitted by [fit_rsf()].
#'
#' @param object A fitted `mlsurv_model` returned by [fit_rsf()].
#' @param newdata A `data.frame` of new observations for prediction. Must include
#'   the same predictor variables used at training time.
#' @param times Optional numeric vector of time points at which to return survival
#'   probabilities. If `NULL`, predictions are returned at the model's internal
#'   time grid.
#' @param ... Additional arguments forwarded to
#'   \code{\link[randomForestSRC]{predict.rfsrc}}.
#'
#' @return A base `data.frame` of survival probabilities with one row per
#' observation in `newdata` and columns named `t={time}` (character), containing
#' numeric values in \[0, 1].
#'
#' @details
#' If `times` is provided, the function aligns predictions by selecting, for each
#' requested time, the closest available time from the RSF prediction object's
#' `time.interest` grid (nearest-neighbor matching).
#'
#' @seealso [fit_rsf()], \code{\link[randomForestSRC]{predict.rfsrc}}
#'
#' @examples
#' mod_rsf <- fit_rsf(Surv(time, status) ~ age + celltype + karno,
#'                    data = veteran, ntree = 200)
#'
#' times <- c(100, 200, 300)
#' pred_probs <- predict_rsf(mod_rsf, newdata = veteran[1:5, ], times = times)
#' print(round(pred_probs, 3))
#'
#' @export

predict_rsf <- function(object, newdata, times = NULL, ...) {


  if (!is.null(object$learner) && object$learner != "rsf") {
    warning("Object passed to predict_rsf() may not come from fit_rsf().")
  }

  stopifnot(requireNamespace("randomForestSRC", quietly = TRUE))

  pred <- randomForestSRC::predict.rfsrc(object$model, newdata = newdata)
  surv_mat <- pred$survival
  time_points <- pred$time.interest

  # interpolation to align predictions to `times`
  if (!is.null(times)) {
    idx <- sapply(times, function(t) which.min(abs(time_points - t)))
    survmat <- surv_mat[, idx, drop = FALSE]
    colnames(survmat) <- paste0("t=", times)
    } else {
    colnames(survmat) <- paste0("t=", time_points)
    }

  as.data.frame(survmat)
}




#' Tune Random Survival Forest Hyperparameters (Cross-Validation)
#'
#' Cross-validates RSF models over a specified parameter grid and selects the
#' best configuration according to the primary metric. Optionally refits the
#' best model on the full dataset.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   survival time used at training). Must be non-negative and finite.
#' @param param_grid A data frame (e.g., from `expand.grid()`) containing columns
#'   for RSF hyperparameters such as `ntree`, `mtry`, and `nodesize`.
#' @param metrics Character vector of metrics to evaluate/optimize
#'   (e.g., `"cindex"`, `"ibs"`). The first entry is the primary selection metric.
#' @param folds Integer; number of cross-validation folds.
#' @param seed Integer random seed for reproducibility.
#' @param refit_best Logical; if `TRUE`, refits the best configuration on the
#'   full data and returns it, with tuning results attached.
#' @param ... Additional arguments forwarded to the underlying engine where applicable.
#'
#' @return If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) of grid
#' results with metric columns and hyperparameters, ordered by the first metric.
#' If `refit_best = TRUE`, a fitted `mlsurv_model` from [fit_rsf()] with
#' attribute `"tuning_results"` containing the full grid results.
#'
#' @details
#' Internally calls `cv_survlearner()` with `fit_rsf()`/`predict_rsf()` so tuning
#' uses the same code paths as production. Typical grids vary `ntree`, `mtry`,
#' and `nodesize`.
#'
#' @seealso [fit_rsf()], [predict_rsf()]
#'
#' @examples
#' \donttest{
#' mod_rsf_best <- tune_rsf(
#'   formula = Surv(time, status) ~ age + celltype + karno,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   param_grid = expand.grid(
#'     ntree = c(200, 500),
#'     mtry = c(1, 2, 3),
#'     nodesize = c(5, 15)
#'   ),
#'   metrics = c("cindex", "ibs"),
#'   folds = 3,
#'   refit_best = TRUE
#' )
#'
#' summary(mod_rsf_best)
#' }
#' @export

tune_rsf <- function(formula, data, times,
                       param_grid = expand.grid(
                         ntree = c(200, 500),
                         mtry = c(1, 2, 3),
                         nodesize = c(5, 15)
                       ),
                       metrics = c("cindex", "ibs"),
                       folds = 5,
                       seed = 123,
                       refit_best = FALSE,
                       ...) {

  results <- purrr::pmap_dfr(param_grid, function(ntree, mtry, nodesize) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_rsf,
      pred_fun = predict_rsf,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))  # like other learners
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_rsf(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      ntree = best_row$ntree,
      mtry = best_row$mtry,
      nodesize = best_row$nodesize,
      ...
    )
    attr(best_model, "tuning_results") <- results  # optional metadata
    return(best_model)
  }
}
