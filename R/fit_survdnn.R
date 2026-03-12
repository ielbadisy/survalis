#' Fit a Deep Neural Network Survival Model (mlsurv_model-compatible)
#'
#' Fits a deep neural network survival model using \pkg{survdnn} (backed by
#' \pkg{torch}). Supports the losses currently exposed by \pkg{survdnn}
#' (`"cox"`, `"cox_l2"`, `"aft"`, and `"coxtime"`). Returns an
#' `mlsurv_model` object that integrates with the `survalis` evaluation and
#' cross-validation pipeline.
#'
#' @section Engine:
#' Uses \strong{survdnn::survdnn} with \pkg{torch}. The wrapper exposes the
#' core training arguments used by the engine, including optimizer choice,
#' dropout, batch normalization, callbacks, device selection, and missing-value
#' handling.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left-hand side must be a `Surv()` object from the `survival` package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param loss Loss function name understood by \pkg{survdnn}. One of
#'   `"cox"`, `"cox_l2"`, `"aft"`, or `"coxtime"`.
#' @param hidden Integer vector of hidden layer sizes (e.g., `c(32L, 32L, 16L)`).
#' @param activation Activation function name. Supported options depend on the
#'   installed \pkg{survdnn} version and include `"relu"`, `"leaky_relu"`,
#'   `"tanh"`, `"sigmoid"`, `"gelu"`, `"elu"`, and `"softplus"`.
#' @param lr Learning rate.
#' @param epochs Number of training epochs.
#' @param optimizer Optimizer name. One of `"adam"`, `"adamw"`, `"sgd"`,
#'   `"rmsprop"`, or `"adagrad"`.
#' @param optim_args Optional named list of extra optimizer arguments passed to
#'   \pkg{torch} via \pkg{survdnn}.
#' @param verbose Logical; print training progress.
#' @param dropout Numeric dropout rate in `[0, 1]`.
#' @param batch_norm Logical; whether to use batch normalization in hidden layers.
#' @param callbacks Optional list of callback functions used by \pkg{survdnn}.
#' @param .seed Optional integer random seed passed through to \pkg{survdnn}.
#' @param .device Computation device. One of `"auto"`, `"cpu"`, or `"cuda"`.
#' @param na_action How to handle missing values. One of `"omit"` or `"fail"`.
#' @param ... Additional arguments forwarded to the underlying engine.
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
                        optimizer = "adam",
                        optim_args = list(),
                        verbose = FALSE,
                        dropout = 0.3,
                        batch_norm = TRUE,
                        callbacks = NULL,
                        .seed = NULL,
                        .device = "auto",
                        na_action = "omit",
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
    optimizer = optimizer,
    optim_args = optim_args,
    verbose = verbose,
    dropout = dropout,
    batch_norm = batch_norm,
    callbacks = callbacks,
    .seed = .seed,
    .device = .device,
    na_action = na_action,
    ...
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
#'   survival time used at training). Required for `type = "survival"` and
#'   `type = "risk"`. Ignored for `type = "lp"`.
#' @param type Prediction type. `"survival"` returns survival probabilities,
#'   `"lp"` returns the engine's linear predictor / latent score, and `"risk"`
#'   returns `1 - S(t)` at a single requested time.
#' @param ... Additional arguments forwarded to `predict.survdnn()`.
#'
#' @return
#' - If `type = "survival"`, a base `data.frame` with one row per observation
#'   in `newdata` and columns named `t={time}` containing values in `[0, 1]`.
#' - If `type = "lp"` or `type = "risk"`, a numeric vector.
#'
#' @details
#' Delegates to the installed \pkg{survdnn} prediction method. The default
#' remains `type = "survival"` so the result stays compatible with
#' `cv_survlearner()` and the rest of the `survalis` evaluation pipeline.
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

predict_survdnn <- function(object, newdata, times = NULL,
                            type = c("survival", "lp", "risk"), ...) {


  if (!is.null(object$learner) && object$learner != "survdnn") {
    warning("Object passed to predict_survdnn() may not come from fit_survdnn().")
  }

  type <- match.arg(type)

  predict(object$model, newdata = newdata, type = type, times = times, ...)
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
#'   `hidden` (list of integer vectors), `lr` (numeric), `activation`
#'   (character), `epochs` (integer), `loss` (character), `optimizer`
#'   (character), `dropout` (numeric), and `batch_norm` (logical).
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
#' production. Hyperparameters are combined via `tidyr::crossing()` and each
#' row is passed through to `fit_survdnn()`, so the grid can include any
#' supported engine argument exposed by this wrapper.
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
#'   loss = c("cox", "coxtime"),
#'   optimizer = "adam",
#'   dropout = c(0.1, 0.3)
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
                           loss = c("cox", "aft"),
                           optimizer = "adam",
                           dropout = c(0.1, 0.3),
                           batch_norm = c(TRUE)
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 3,
                         seed = 42,
                         refit_best = FALSE,
                         ...) {

  stopifnot(!missing(formula), !missing(data), !missing(times), !missing(param_grid))
  stopifnot(is.list(param_grid))

  param_df <- tidyr::crossing(!!!param_grid)
  extra_args <- list(...)
  higher_is_better <- metrics[1] %in% c("cindex", "auc", "accuracy")

  as_tibble_row <- function(params) {
    tibble::as_tibble(purrr::imap(params, function(x, nm) {
      if (is.list(param_df[[nm]]) || length(x) != 1L || is.list(x)) {
        list(x)
      } else {
        x
      }
    }))
  }

  results <- purrr::map_dfr(seq_len(nrow(param_df)), function(i) {
    params <- purrr::map(param_df[i, , drop = FALSE], ~ .x[[1]])
    set.seed(seed)
    cv_args <- c(
      list(
        formula = formula,
        data = data,
        fit_fun = fit_survdnn,
        pred_fun = predict_survdnn,
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        verbose = FALSE
      ),
      params,
      extra_args
    )
    cv_results <- do.call(cv_survlearner, cv_args)

    summary <- cv_summary(cv_results)

    as_tibble_row(params) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric,
                           values_from = mean)
      )
  })

  if (higher_is_better) {
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
    best_params <- purrr::map(best_row[names(param_df)], ~ .x[[1]])
    fit_args <- c(
      list(
        formula = formula,
        data = tidyr::drop_na(data, all.vars(formula)),
        verbose = FALSE
      ),
      best_params,
      extra_args
    )
    best_model <- do.call(fit_survdnn, fit_args)
    attr(best_model, "tuning_results") <- results
    return(best_model)
  }
}
