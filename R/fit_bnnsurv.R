#' Fit a kNN–Ensemble Survival Model (bnnSurvival)
#'
#' Fits an ensemble of k–nearest neighbour survival learners via
#' \pkg{bnnSurvival}. The fitted object is standardized to the
#' \code{mlsurv_model} contract used in \pkg{survalis}.
#'
#' @section Engine:
#' Uses \strong{bnnSurvival::bnnSurvival}. This wrapper calls the engine via
#' \code{requireNamespace("bnnSurvival", quietly = TRUE)} and stores the native
#' model in \code{$model}.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing all variables referenced in \code{formula}.
#' @param k Integer, number of neighbours for each base learner. Default \code{20}.
#' @param num_base_learners Integer, number of base learners in the ensemble.
#'   Default \code{10}.
#' @param num_features_per_base_learner Integer or \code{NULL}. If not \code{NULL},
#'   each base learner samples this many features. Passed to the engine.
#' @param metric Character distance metric for neighbour search (for example,
#'   \code{"mahalanobis"}). Passed to the engine.
#' @param weighting_function Function used to weight neighbours. Defaults to a
#'   constant weighting \code{function(x) x * 0 + 1}. Passed to the engine.
#' @param replace Logical; sample with replacement when drawing observations for
#'   a base learner. Passed to the engine.
#' @param sample_fraction Optional numeric in \code{(0, 1]}; fraction of rows used
#'   by each base learner. Passed to the engine.
#'
#' @return An object of class \code{"mlsurv_model"} with elements:
#' \describe{
#'   \item{model}{The underlying \pkg{bnnSurvival} fit.}
#'   \item{learner}{Scalar \code{"bnnsurv"}.}
#'   \item{engine}{Scalar \code{"bnnSurvival"}.}
#'   \item{formula, data}{Inputs preserved for downstream use.}
#'   \item{time, status}{Character names of the survival outcome fields.}
#' }
#'
#' @details
#' The native engine returns full survival curves on an internal time grid.
#' See \code{\link{predict_bnnsurv}} for how these are post–processed and
#' (optionally) interpolated to user–requested times.
#'
#' @seealso [predict_bnnsurv()], [tune_bnnsurv()]
#'
#' @examplesIf requireNamespace("bnnSurvival", quietly = TRUE)
#'
#' mod <- fit_bnnsurv(
#'   Surv(time, status) ~ age + karno + diagtime + prior,
#'   data = veteran
#' )
#' head(predict_bnnsurv(mod, newdata = veteran[1:5, ], times = c(50, 100)))
#'
#' @export

fit_bnnsurv <- function(formula, data,
                        k = 5,
                        num_base_learners = 10,
                        num_features_per_base_learner = NULL,
                        metric = "mahalanobis",
                        weighting_function = function(x) x * 0 + 1,
                        replace = TRUE,
                        sample_fraction = NULL) {
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
    sample_fraction = sample_fraction
  )

  structure(list(
    model   = model,
    learner = "bnnsurv",
    engine  = "bnnSurvival",
    formula = formula,
    data    = data,
    time    = all.vars(formula)[[2]],
    status  = all.vars(formula)[[3]]
  ), class = "mlsurv_model")
}


#' Predict Survival with a bnnSurvival Model
#'
#' Generates survival probabilities from an \code{mlsurv_model} fitted by
#' \code{\link{fit_bnnsurv}}. If \code{times} is supplied, survival curves are
#' linearly interpolated from the engine's internal time grid to those times.
#'
#' @param object A fitted \code{mlsurv_model} returned by \code{fit_bnnsurv()}.
#' @param newdata A data frame of new observations for prediction.
#' @param times Optional numeric vector of evaluation time points. If \code{NULL}
#'   (default), the function returns survival on the engine's native grid.
#'
#' @return A base \code{data.frame} with one row per observation and columns
#'   named \code{"t=<time>"} containing survival probabilities in \code{[0, 1]}.
#'   Probabilities are clipped to \code{[0, 1]} and made non–increasing over time.
#'
#' @details
#' Internally, predictions are obtained via \code{bnnSurvival::predict()}.
#' The returned survival matrix is post–processed to enforce monotonicity
#' (cumulative minimum over time). When \code{times} is provided, values are
#' obtained by \code{stats::approx()} (linear interpolation, rule = 2).
#'
#' @seealso [fit_bnnsurv()], [tune_bnnsurv()]
#'
#' @examplesIf requireNamespace("bnnSurvival", quietly = TRUE)
#'
#' mod <- fit_bnnsurv(Surv(time, status) ~ age + karno + diagtime + prior, data = veteran)
#' pred <- predict_bnnsurv(mod, newdata = veteran[1:3, ], times = c(50, 100, 200))
#' pred
#'
#' @export

predict_bnnsurv <- function(object, newdata, times = NULL) {
  stopifnot(object$learner == "bnnsurv")
  stopifnot(requireNamespace("bnnSurvival", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  # Get full survival curves from bnnSurvival
  pr   <- bnnSurvival::predict(object$model, newdata)
  tg   <- bnnSurvival::timepoints(pr)        # numeric time grid
  Smat <- bnnSurvival::predictions(pr)       # nrow(newdata) x length(tg)

  # ensure probabilities are within [0,1] and non-increasing over time grid
  Smat <- pmin(pmax(Smat, 0), 1)
  Smat <- t(apply(Smat, 1L, cummin))

  if (is.null(times)) {
    return(.finalize_survmat(Smat, times = tg))
  }

  stopifnot(is.numeric(times), length(times) > 0L, all(is.finite(times)))

  # interpolate from model grid (tg) to requested times
  survmat <- t(apply(Smat, 1L, function(row)
    stats::approx(x = tg, y = row, xout = times, method = "linear", rule = 2)$y
  ))
  .finalize_survmat(survmat, times = times)
}


#' Tune bnnSurvival Hyperparameters (Cross-Validation)
#'
#' Cross-validates \code{\link{fit_bnnsurv}} over a user-supplied grid and
#' aggregates metrics (for example, \code{"cindex"}, \code{"ibs"}). Optionally
#' refits and returns the best configuration.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data Training data frame.
#' @param times Numeric vector of evaluation time points used during CV.
#' @param param_grid A data frame (for example from \code{expand.grid()}) with
#'   columns \code{k}, \code{num_base_learners}, and \code{sample_fraction}.
#'   Defaults cover a small illustrative search.
#' @param metrics Character vector of metric names to compute and summarize.
#'   The first entry is treated as the primary ranking metric.
#' @param folds Integer number of cross-validation folds. Default \code{5}.
#' @param seed Integer random seed for reproducibility.
#' @param refit_best Logical; if \code{TRUE}, refits the best configuration on
#'   the full data and returns that model. Otherwise returns the results table.
#'
#' @return If \code{refit_best = FALSE}, a \code{tibble} with one row per
#'   configuration containing \code{metrics} and the tuning parameters, with a
#'   \code{failed} flag for combinations that errored during CV. If
#'   \code{refit_best = TRUE}, an \code{mlsurv_model} refit at the selected
#'   hyperparameters.
#'
#' @details
#' Evaluation is performed by \code{cv_survlearner()} using
#' \code{fit_bnnsurv()} and \code{predict_bnnsurv()} to match production code
#' paths. Results are ordered by the first entry in \code{metrics}.
#'
#' @seealso [fit_bnnsurv()], [predict_bnnsurv()], [cv_survlearner()]
#'
#' @examplesIf requireNamespace("bnnSurvival", quietly = TRUE)
#' \donttest{
#' grid <- expand.grid(
#'   k = c(2, 3),
#'   num_base_learners = c(5, 10),
#'   sample_fraction = c(0.5, 1),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Performance table
#' res <- tune_bnnsurv(
#'   formula = Surv(time, status) ~ age + karno + diagtime + prior,
#'   data = veteran,
#'   times = c(100, 200),
#'   param_grid = grid,
#'   refit_best = FALSE
#' )
#' res
#'
#' # Best refitted model
#' best <- tune_bnnsurv(
#'   formula = Surv(time, status) ~ age + karno + diagtime + prior,
#'   data = veteran,
#'   times = c(100, 200),
#'   param_grid = grid,
#'   refit_best = TRUE
#' )
#' head(predict_bnnsurv(best, newdata = veteran[1:5, ], times = c(50, 100)))
#' }
#' 
#' @export

tune_bnnsurv <- function(formula, data, times,
                         param_grid = expand.grid(
                           k = c(2, 3),
                           num_base_learners = c(30, 50),
                           sample_fraction = c(0.5, 1),
                           stringsAsFactors = FALSE
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5,
                         seed = 123,
                         ncores = 1,
                         refit_best = FALSE) {
  results <- purrr::pmap_dfr(param_grid, function(k, num_base_learners, sample_fraction) {
    cv_result <- tryCatch({
      cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_bnnsurv,
        pred_fun = predict_bnnsurv,
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        ncores = ncores,
        k = k,
        num_base_learners = num_base_learners,
        sample_fraction = sample_fraction
      )
    }, error = function(e) NULL)

    if (is.null(cv_result)) {
      return(tibble::tibble(
        k = k,
        num_base_learners = num_base_learners,
        sample_fraction = sample_fraction,
        failed = TRUE
      ))
    }

    summary <- cv_summary(cv_result)

    tibble::tibble(
      k = k,
      num_base_learners = num_base_learners,
      sample_fraction = sample_fraction,
      failed = FALSE
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  results <- results[!results[["failed"]], , drop = FALSE]

  if (nrow(results) == 0) {
    warning("All tuning combinations failed.")
    return(tibble::tibble())
  }

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    results <- dplyr::arrange(results, dplyr::desc(.data[[metrics[1]]]))
  } else {
    results <- dplyr::arrange(results, .data[[metrics[1]]])
  }

  if (refit_best) {
    best <- results[1, ]
    model <- fit_bnnsurv(
      formula = formula,
      data = data,
      k = best$k,
      num_base_learners = best$num_base_learners,
      sample_fraction = best$sample_fraction
    )
    return(model)
  }

  return(results)
}
