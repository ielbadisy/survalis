#' Fit a Deep Neural Network Survival Model (mlsurv_model-compatible)
#'
#' Fits a deep neural network survival model using \pkg{survdnn} (backed by
#' \pkg{torch}). Supports common survival losses (e.g., `"cox"`, `"aft"`,
#' `"coxtime"` depending on your \pkg{survdnn} build). Returns an
#' `mlsurv_model` object that integrates with the `survalis` evaluation and
#' cross-validation pipeline.
#'
#' @section Engine:
#' Uses \strong{survdnn::survdnn} with \pkg{torch}. Hidden layer widths,
#' activation, learning rate, epochs, and verbosity are passed through.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left-hand side must be a `Surv()` object from the `survival` package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param loss Loss function name understood by \pkg{survdnn} (e.g., `"cox"`,
#'   `"aft"`, `"coxtime"`).
#' @param hidden Integer vector of hidden layer sizes (e.g., `c(32L, 32L, 16L)`).
#' @param activation Activation function name (e.g., `"relu"`, `"gelu"`, `"tanh"`),
#'   as supported by \pkg{survdnn}.
#' @param lr Learning rate.
#' @param epochs Number of training epochs.
#' @param verbose Logical; print training progress
#'
#' @param ... Additional arguments forwarded to the underlying engine where applicable.
#'
#' @return An object of class `mlsurv_model`, a named list with elements:
#' \describe{
#'   \item{model}{The underlying fitted \pkg{survdnn} model.}
#'   \item{learner}{Character scalar identifying the learner (`"survdnn"`).}
#'   \item{engine}{Character scalar naming the engine (`"survdnn"`).}
#'   \item{formula}{The original survival formula.}
#'   \item{data}{The training dataset (or a minimal subset needed for prediction).}
#'   \item{time}{Name of the survival time variable.}
#'   \item{status}{Name of the event indicator (1 = event, 0 = censored).}
#' }
#'
#' @details
#' \strong{Design contract.} All `fit_*()` functions in `survalis`:
#' (i) return a named list with `model`, `learner`, `engine`, `formula`,
#' `data`, `time`, and `status`; and (ii) retain information required by
#' `predict_*()` to build a consistent design for new data.
#'
#' @seealso [predict_survdnn()], [tune_survdnn()]
#'
#' @references
#' \pkg{survdnn} documentation; \pkg{torch} for deep learning in R.
#'
#' @examples
#'
#' mod <- fit_survdnn(Surv(time, status) ~ age + karno + celltype,
#'                    data = veteran, loss = "cox", epochs = 50, verbose = FALSE)
#'
#' pred <- predict_survdnn(mod, newdata = veteran[1:5, ], times = c(30, 90, 180))
#' print(pred)
#'
#' @export

fit_survdnn <- function(formula, data,
                        loss = "cox",
                        hidden = c(32L, 32L, 16L),
                        activation = "relu",
                        lr = 1e-4,
                        epochs = 300L,
                        verbose = FALSE,
                        ...) {

  stopifnot(requireNamespace("survdnn", quietly = TRUE))
  stopifnot(requireNamespace("torch", quietly = TRUE))

  model <- survdnn::survdnn(
    formula = formula,
    data = data,
    loss = loss,
    hidden = hidden,
    activation = activation,
    lr = lr,
    epochs = epochs,
    verbose = verbose
  )

  time_status <- all.vars(formula[[2]])

  structure(list(
    model = model,
    learner = "survdnn",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model", engine = "survdnn")
  }




#' Predict Survival Probabilities with a DNN Survival Model
#'
#' Generates survival probabilities at specified time points using a fitted
#' deep neural network survival model (`mlsurv_model`).
#'
#' @param object A fitted `mlsurv_model` returned by [fit_survdnn()].
#' @param newdata A `data.frame` of new observations for prediction. Must include
#'   the same predictor variables used at training time.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   survival time used at training). Must be non-negative and finite.
#' @param ... Additional arguments forwarded to the underlying engine where applicable.
#'
#' @return A base `data.frame` with one row per observation in `newdata` and
#' columns named `t={time}` (character), containing numeric values in \[0, 1].
#'
#' @details
#' Delegates to `predict(object$model, type = "survival", times = times)`.
#'
#' @seealso [fit_survdnn()], [tune_survdnn()]
#'
#' @examples
#'
#' mod <- fit_survdnn(Surv(time, status) ~ age + karno + celltype,
#'                    data = veteran, loss = "cox", epochs = 50, verbose = FALSE)
#'
#' pred <- predict_survdnn(mod, newdata = veteran[1:5, ], times = c(30, 90, 180))
#' print(pred)
#'
#' @export

predict_survdnn <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "survdnn") {
    warning("Object passed to predict_survdnn() may not come from fit_survdnn().")
  }

  predict(object$model, newdata = newdata, type = "survival", times = times)
}


#' Tune Deep Neural Network Survival Models (Cross-Validation)
#'
#' Cross-validates \pkg{survdnn}-based models over a user-specified grid and
#' selects the best configuration by the primary metric. Optionally refits the
#' best model on the full dataset.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left-hand side must be a `Surv()` object from the `survival` package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   survival time used at training). Must be non-negative and finite.
#' @param param_grid A named list of candidate hyperparameters. Typical entries:
#'   `hidden` (list of integer vectors), `lr` (numeric), `activation` (character),
#'   `epochs` (integer), and `loss` (character).
#' @param metrics Character vector of metrics to evaluate/optimize
#'   (e.g., `"cindex"`, `"ibs"`). The first entry is used as the primary
#'   selection metric.
#' @param folds Number of cross-validation folds.
#' @param seed Integer random seed for reproducibility.
#' @param refit_best Logical; if `TRUE`, refits the best configuration on the
#'   full data and returns it as the result (with tuning attributes).
#' @param ... Additional arguments forwarded to the underlying engine where applicable.
#'
#' @return If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) with
#' hyperparameter settings and metric columns, ordered by the first metric.
#' If `refit_best = TRUE`, a fitted `mlsurv_model` from [fit_survdnn()] with
#' attribute `"tuning_results"` containing the full grid results.
#'
#' @details
#' \strong{Evaluation.} Internally calls `cv_survlearner()` with
#' `fit_survdnn()`/`predict_survdnn()` so tuning uses the same code paths as
#' production. Hyperparameters are combined via `tidyr::crossing()` and
#' iterated with `purrr::pmap_dfr()`.
#'
#' @seealso [fit_survdnn()], [predict_survdnn()]
#'
#' @examples
#' \donttest{
#'
#' grid <- list(
#'   hidden = list(c(16), c(32, 16)),
#'   lr = c(1e-4, 5e-4),
#'   activation = c("relu", "tanh"),
#'   epochs = c(300),
#'   loss = c("cox", "coxtime")
#' )
#'
#' mod <- tune_survdnn(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(90),
#'   metrics = c("cindex", "ibs"),
#'   param_grid = grid,
#'   refit_best = TRUE
#' )
#'
#' summary(mod)
#' }
#'
#' @export

tune_survdnn <- function(formula, data, times,
                         param_grid = list(
                           hidden = list(c(32, 16), c(64, 32, 16)),
                           lr = c(1e-3, 5e-4),
                           activation = c("relu", "gelu"),
                           epochs = c(100, 200),
                           loss = c("cox", "aft")
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 3,
                         seed = 42,
                         refit_best = FALSE,
                         ...) {

  stopifnot(!missing(formula), !missing(data), !missing(times), !missing(param_grid))
  stopifnot(is.list(param_grid))

  param_df <- tidyr::crossing(!!!param_grid)

  results <- purrr::pmap_dfr(param_df, function(hidden, lr, activation, epochs, loss) {
    set.seed(seed)
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_survdnn,
      pred_fun = predict_survdnn,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      hidden = hidden,
      lr = lr,
      activation = activation,
      epochs = epochs,
      loss = loss,
      verbose = FALSE,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      hidden = list(hidden),
      lr = lr,
      activation = activation,
      epochs = epochs,
      loss = loss
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric,
                           values_from = mean)
      )
  })

  results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_survdnn(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      hidden = best_row$hidden[[1]],
      lr = best_row$lr,
      activation = best_row$activation,
      epochs = best_row$epochs,
      loss = best_row$loss,
      verbose = FALSE,
      ...
    )
    attr(best_model, "tuning_results") <- results
    return(best_model)
  }
}

