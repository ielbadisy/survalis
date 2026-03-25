#' Fit a Survival SVM Model (mlsurv_model-compatible)
#'
#' Fits a survival support vector machine using \pkg{survivalsvm}. The default
#' setup uses the regression-type loss with quadratic programming optimization
#' and an additive or linear kernel, but all \pkg{survivalsvm} options are exposed.
#' Returns an `mlsurv_model` object that integrates with the `survalis`
#' evaluation and cross-validation pipeline.
#'
#' @section Engine:
#' Uses \strong{survivalsvm::survivalsvm}. Some optimization methods (e.g.,
#' `"quadprog"`) require additional system dependencies and the \pkg{quadprog}
#' package to be installed.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left-hand side must be a `Surv()` object from the `survival` package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param type SVM loss type used by \code{survivalsvm} (e.g., \code{"regression"}).
#' @param gamma.mu Regularization parameter for the margin/hinge component.
#' @param opt.meth Optimization method (e.g., \code{"quadprog"}).
#' @param kernel Kernel type (e.g., \code{"lin_kernel"}, \code{"add_kernel"}).
#' @param diff.meth Optional differentiation method passed to the engine.
#'
#' @return An object of class `mlsurv_model`, a named list with elements:
#' \describe{
#'   \item{model}{The underlying fitted \pkg{survivalsvm} model.}
#'   \item{learner}{Character scalar identifying the learner (`"survsvm"`).}
#'   \item{engine}{Character scalar naming the engine (`"survivalsvm"`).}
#'   \item{formula}{The original survival formula.}
#'   \item{data}{The training dataset (or a minimal subset needed for prediction).}
#'   \item{time}{Name of the survival time variable.}
#'   \item{status}{Name of the event indicator (1 = event, 0 = censored).}
#' }
#'
#' @details
#' \strong{Design contract.} All `fit_*()` functions in `survalis`:
#' (i) return a named list with `model`, `learner`, `engine`, `formula`,
#' `data`, `time`, and `status`; (ii) preserve `terms(formula)` or equivalent
#' for consistent prediction design; and (iii) keep engine arguments required
#' downstream by `predict_*()`.
#'
#' @seealso [predict_survsvm()], [tune_survsvm()]
#'
#' @references
#' Binder H, et al. (2009). Survival Support Vector Machines (and related work).
#' \pkg{survivalsvm} package documentation.
#'
#' @examplesIf requireNamespace("survivalsvm", quietly = TRUE)
#' mod_svm <- fit_survsvm(Surv(time, status) ~ age + celltype + karno,
#'                            data = veteran,
#'                            type = "regression",
#'                            gamma.mu = 0.1,
#'                            kernel = "lin_kernel")
#'
#' times <- c(100, 300, 500)
#' predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "exp")
#' predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "weibull", shape = 1.5)
#'
#' #-------------------
#'
#' cv_results_svm <- cv_survlearner(
#'   formula = Surv(time, status) ~ age + celltype + karno,
#'   data = veteran,
#'   fit_fun = fit_survsvm,
#'   pred_fun = predict_survsvm,
#'   times = c(100, 300, 500),
#'   metrics = c("cindex", "ibs"),
#'   folds = 5,
#'   seed = 42,
#'   gamma.mu = 0.1,
#'   kernel = "lin_kernel")
#'
#' print(cv_results_svm)
#' cv_summary(cv_results_svm)
#' cv_plot(cv_results_svm)
#'
#' @export

fit_survsvm <- function(formula, data,
                            type = "regression",
                            gamma.mu = 0.1,
                            opt.meth = "quadprog",
                            kernel = "add_kernel",
                            diff.meth = NULL) {



  stopifnot(requireNamespace("survivalsvm", quietly = TRUE)) ## need also to install quadprog

  model <- survivalsvm::survivalsvm(
    formula = formula,
    data = data,
    type = type,
    gamma.mu = gamma.mu,
    opt.meth = opt.meth,
    kernel = kernel,
    diff.meth = diff.meth
  )

  structure(list(
    model = model,
    learner = "survsvm",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "survivalsvm")
}



#' Predict Survival Probabilities with Survival SVM
#'
#' Generates survival probabilities at specified times from a fitted survival
#' SVM (`mlsurv_model`). Since many SVM variants output a single predicted time
#' (or rank), this function maps predicted times to survival curves using either
#' an Exponential or Weibull parametric assumption.
#'
#' @param object A fitted `mlsurv_model` returned by [fit_survsvm()].
#' @param newdata A `data.frame` of new observations for prediction. Must include
#'   the same predictor variables used at training time.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   survival time used at training). Must be non-negative and finite.
#' @param dist Parametric family used to map predicted times to survival curves:
#'   \code{"exp"} (exponential) or \code{"weibull"}.
#' @param shape Weibull shape parameter (used only if \code{dist = "weibull"}).
#'
#' @return A base `data.frame` with one row per observation in `newdata` and
#' columns named `t={time}` (character), containing numeric values in \[0, 1].
#'
#' @details
#' \strong{Parametric mapping.} Let \eqn{\hat{T}} be the predicted time from
#' the SVM model (per row). For a requested evaluation time \eqn{t}:
#' \itemize{
#'   \item Exponential: \eqn{S(t) = \exp(-t / \hat{T})}
#'   \item Weibull:     \eqn{S(t) = \exp\{-(t / \hat{T})^{\kappa}\}}, where \eqn{\kappa} is \code{shape}
#' }
#' This mapping provides calibrated-looking survival probabilities but is a
#' modeling assumption external to the SVM fit; verify adequacy in practice.
#'
#' @seealso [fit_survsvm()], [tune_survsvm()]
#'
#' @examplesIf requireNamespace("survivalsvm", quietly = TRUE)
#' mod_svm <- fit_survsvm(Surv(time, status) ~ age + celltype + karno,
#'                            data = veteran,
#'                            type = "regression",
#'                            gamma.mu = 0.1,
#'                            kernel = "lin_kernel")
#'
#' times <- c(100, 300, 500)
#' predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "exp")
#' predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "weibull", shape = 1.5)
#'
#' @export

predict_survsvm <- function(object, newdata, times, dist = "exp", shape = 1) {

  if (!is.null(object$learner) && object$learner != "survsvm") {
    warning("Object passed to predict_survsvm() may not come from fit_survsvm().")
  }

  stopifnot(requireNamespace("survivalsvm", quietly = TRUE))

  newdata <- as.data.frame(newdata)
  pred <- predict(object$model, newdata = newdata)
  predicted_times <- as.numeric(pred$predicted)

  dist <- match.arg(dist, choices = c("exp", "weibull"))

  survmat <- switch(dist,
                  "exp" = sapply(times, function(t) exp(-t / predicted_times)),
                  "weibull" = sapply(times, function(t) exp(- (t / predicted_times)^shape))
  )

  .finalize_survmat(survmat, times = times)
}

#' Tune Survival SVM Hyperparameters (Cross-Validation)
#'
#' Cross-validates Survival SVM models over a user-specified grid and selects
#' the best configuration based on the chosen metric. Optionally refits the
#' best model on the full dataset.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left-hand side must be a `Surv()` object from the `survival` package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   survival time used at training). Must be non-negative and finite.
#' @param metrics Character vector of metrics to evaluate/optimize
#'   (e.g., `"cindex"`, `"ibs"`). The first entry is the primary selection metric.
#' @param param_grid A named list of candidate hyperparameters; typical entries
#'   include `gamma.mu` and `kernel`. Use \code{list(gamma.mu = c(...), kernel = c(...))}.
#' @param folds Number of cross-validation folds.
#' @param seed Integer random seed for reproducibility.
#' @param refit_best Logical; if `TRUE`, refits the best configuration on the
#'   full data and returns it as the result (with tuning attributes).
#' @param dist Parametric mapping for predictions during CV: `"exp"` or `"weibull"`.
#' @param shape Weibull shape parameter if \code{dist = "weibull"}.
#'
#' @return If `refit_best = FALSE`, a `data.frame` (class `"tune_surv"`) of grid
#' results with metric columns and tuning parameters. If `refit_best = TRUE`, a
#' fitted `mlsurv_model` (class augmented with `"tune_surv"`) with the full
#' results attached in `attr(, "tuning_results")`.
#'
#' @details
#' \strong{Evaluation.} Internally calls `cv_survlearner()` with
#' `fit_survsvm()`/`predict_survsvm()` to keep code paths consistent with
#' production usage. The prediction step applies the specified \code{dist}/\code{shape}
#' mapping to convert predicted times into survival probabilities at `times`.
#'
#' @seealso [fit_survsvm()], [predict_survsvm()]
#'
#' @examplesIf requireNamespace("survivalsvm", quietly = TRUE)
#' grid <- list(
#'   gamma.mu = c(0.01, 0.1),
#'   kernel = c("lin_kernel", "add_kernel")
#' )
#'
#' res_svm <- tune_survsvm(
#'   formula = Surv(time, status) ~ age + celltype + karno,
#'   data = veteran,
#'   times = c(100, 300, 500),
#'   metrics = c("cindex", "ibs"),
#'   param_grid = grid,
#'   folds = 3,
#'   refit_best = TRUE
#' )
#'
#' summary(res_svm)
#'
#' res_svm <- tune_survsvm(
#'   formula = Surv(time, status) ~ age + celltype + karno,
#'   data = veteran,
#'   times = c(100, 300, 500),
#'   metrics = c("cindex", "ibs"),
#'   param_grid = grid,
#'   folds = 3,
#'   refit_best = FALSE
#' )
#'
#' res_svm
#'
#' @export

tune_survsvm <- function(formula, data, times,
                             metrics = "cindex",
                             param_grid,
                             folds = 5, seed = 42, ncores = 1,
                             refit_best = FALSE,
                             dist = "exp", shape = 1) {
  stopifnot(is.list(param_grid))
  param_df <- tidyr::crossing(!!!param_grid)

  # Fixed safe defaults
  fixed_type <- "regression"
  fixed_opt_meth <- "quadprog"
  fixed_diff_meth <- NULL

  results <- purrr::pmap_dfr(param_df, function(gamma.mu, kernel) {
    set.seed(seed)
    tryCatch({
      res_cv <- cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_survsvm,
        pred_fun = function(...) predict_survsvm(..., dist = dist, shape = shape),
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        ncores = ncores,
        gamma.mu = gamma.mu,
        kernel = kernel,
        type = fixed_type,
        opt.meth = fixed_opt_meth,
        diff.meth = fixed_diff_meth
      )

      summary <- cv_summary(res_cv)

      tibble::tibble(
        gamma.mu = gamma.mu,
        kernel = kernel
      ) |>
        dplyr::bind_cols(
          tidyr::pivot_wider(summary[, c("metric", "mean")],
                             names_from = metric, values_from = mean)
        )
    }, error = function(e) {
      invisible(cat(glue("Skipping (gamma.mu={gamma.mu}, kernel={kernel}): {e$message}\n"),
                    file = "errors.log", append = TRUE)) ## keep a trace without showing it on the console
      NULL
    })
  })

  if (nrow(results) == 0) {
    warning("All tuning combinations failed.")
    return(NULL)
  }

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    results <- dplyr::arrange(results, dplyr::desc(.data[[metrics[1]]]))
  } else {
    results <- dplyr::arrange(results, .data[[metrics[1]]])
  }

  if (refit_best) {
    best <- results[1, ]
    model <- fit_survsvm(
      formula = formula,
      data = data,
      gamma.mu = best$gamma.mu,
      kernel = best$kernel,
      type = fixed_type,
      opt.meth = fixed_opt_meth,
      diff.meth = fixed_diff_meth
    )
    attr(model, "tuning_results") <- results
    class(model) <- c("tune_surv", class(model))
    return(model)
  }

  class(results) <- c("tune_surv", class(results))
  return(results)
}
