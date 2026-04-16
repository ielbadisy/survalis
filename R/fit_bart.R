#' Fit a Bayesian Additive Regression Trees (BART) Survival Model
#'
#' Fits a right‐censored survival model using \pkg{BART}'s
#' \code{\link[BART]{surv.bart}} and returns an `mlsurv_model` compatible with
#' the \pkg{survalis} pipeline.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#'   Only right‐censored responses are supported.
#' @param data A \code{data.frame} containing variables referenced in \code{formula}.
#' @param K Integer; number of internal time grid intervals used by the BART
#'   survival engine (default \code{3}). Larger values yield a finer internal
#'   evaluation grid.
#' @param ... Further args to BART::surv.bart

#'
#' @details
#' The response must be of class \code{Surv(..., type = "right")}. The fitted
#' object stores the engine's internal evaluation times in \code{$eval_times}
#' (from \code{bart_fit$times}) for downstream prediction alignment.
#'
#' @return An object of class \code{mlsurv_model} with elements:
#' \describe{
#'   \item{model}{Fitted \code{BART::surv.bart} object.}
#'   \item{learner}{\code{"bart"}.}
#'   \item{formula, data}{Original inputs.}
#'   \item{eval_times}{Engine's internal time grid used for prediction.}
#' }
#'
#' @seealso \code{\link{predict_bart}}, \code{\link[BART]{surv.bart}}
#'
#' @examplesIf requireNamespace("BART", quietly = TRUE)
#' ex_data <- veteran[1:40, c("time", "status", "age", "karno", "celltype")]
#' mod_bart <- fit_bart(
#'   Surv(time, status) ~ age + karno + celltype,
#'   data = ex_data,
#'   K = 1,
#'   ntree = 5,
#'   ndpost = 20,
#'   nskip = 5,
#'   mc.cores = 1,
#'   seed = 42
#' )
#' @export


fit_bart <- function(formula, data, K = 3, ...) {

  stopifnot(requireNamespace("BART", quietly = TRUE))

  mf <- model.frame(formula, data)
  y <- model.response(mf)

  if (!inherits(y, "Surv")) stop("Response must be a 'Surv' object.")
  if (attr(y, "type") != "right") stop("Only right-censored data supported.")

  time_vec   <- y[, "time"]
  status_vec <- y[, "status"]
  terms_obj <- stats::terms(formula, data = data)
  xlevels <- stats::.getXlevels(terms_obj, mf)
  x <- stats::model.matrix(terms_obj, data = mf)[, -1, drop = FALSE]
  x_cols <- colnames(x)
  x_contrasts <- attr(x, "contrasts")

  # silence verbose output
  invisible(
    capture.output({
      bart_fit <- BART::surv.bart(
        x.train = x,
        times = time_vec,
        delta = status_vec,
        K = K,
        printevery = 1
      )
    })
  )

  structure(
    list(
      model = bart_fit,
      learner = "bart",
      formula = formula,
      data = data,
      eval_times = bart_fit$times,
      terms = stats::delete.response(terms_obj),
      xlevels = xlevels,
      x_cols = x_cols,
      x_contrasts = x_contrasts
    ),
    class = "mlsurv_model",
    engine = "BART"
  )
}



#' Predict Survival Probabilities from a BART Survival Model
#'
#' Generates survival probability predictions at requested times for new data
#' using a model fitted by \code{\link{fit_bart}}.
#'
#' @param object An \code{mlsurv_model} returned by \code{\link{fit_bart}}.
#' @param newdata A \code{data.frame} of new observations for prediction.
#' @param times Numeric vector of time points at which to estimate survival
#'   probabilities (same scale as the training time).
#'
#' @details
#' The BART engine predicts survival on its internal grid \code{object$eval_times}.
#' Requested \code{times} are aligned to that grid by nearest‐neighbor matching,
#' returning one survival estimate per requested time.
#'
#' @return A base \code{data.frame} with one row per observation in \code{newdata}
#' and columns named \code{"t=<time>"} (character), containing survival
#' probabilities in \[0, 1].
#'
#' @seealso \code{\link{fit_bart}}, \code{\link[BART]{surv.bart}}
#'
#' @examplesIf requireNamespace("BART", quietly = TRUE)
#' ex_data <- veteran[1:40, c("time", "status", "age", "karno", "celltype")]
#' mod_bart <- fit_bart(
#'   Surv(time, status) ~ age + karno + celltype,
#'   data = ex_data,
#'   K = 1,
#'   ntree = 5,
#'   ndpost = 20,
#'   nskip = 5,
#'   mc.cores = 1,
#'   seed = 42
#' )
#' predict_bart(mod_bart, newdata = ex_data[1:5, ], times = c(10, 30, 60))
#' @export


predict_bart <- function(object, newdata, times) {


  if (!is.null(object$learner) && object$learner != "bart") {
    warning("Object passed to predict_bart() may not come from fit_bart().")
  }


  stopifnot(requireNamespace("BART", quietly = TRUE))

  if (!is.null(object$terms) && !is.null(object$x_cols)) {
    mf_new <- stats::model.frame(object$terms, data = newdata,
                                 xlev = object$xlevels, na.action = stats::na.pass)
    newx <- stats::model.matrix(object$terms, data = mf_new,
                                contrasts.arg = object$x_contrasts)
    if ("(Intercept)" %in% colnames(newx)) {
      newx <- newx[, colnames(newx) != "(Intercept)", drop = FALSE]
    }
    missing_cols <- setdiff(object$x_cols, colnames(newx))
    if (length(missing_cols) > 0) {
      for (nm in missing_cols) {
        newx <- cbind(newx, setNames(data.frame(rep(0, nrow(newx))), nm))
      }
    }
    newx <- as.matrix(newx[, object$x_cols, drop = FALSE])
  } else {
    newx <- model.matrix(object$formula, data = newdata)[, -1, drop = FALSE]
  }

  n_obs <- nrow(newx)
  K <- length(object$eval_times)
  if (K < 1) stop("Fitted BART model does not contain valid evaluation times.")

  if (K == 1) {
    newx_with_t <- cbind(t = rep(object$eval_times[1], n_obs), newx)
  } else {
    newx_with_t <- cbind(
      t = object$eval_times,
      newx[rep(seq_len(n_obs), each = K), , drop = FALSE]
    )
  }

  newx_pred <- newx_with_t
  if (!is.null(object$model$rm.const)) {
    newx_pred <- newx_pred[, object$model$rm.const, drop = FALSE]
  }

  expected_p <- length(object$model$treedraws$cutpoints)
  if (ncol(newx_pred) != expected_p) {
    stop("Prediction design matrix has unexpected number of columns for the fitted BART model.")
  }

  invisible(capture.output(
    pred_values <- predict(object$model, newdata = as.matrix(newx_pred))$surv.test.mean
  ))

  if (length(pred_values) == n_obs) {
    prob_matrix <- matrix(pred_values, nrow = n_obs, ncol = 1, byrow = TRUE)
    eval_times <- object$eval_times[1]
  } else if (length(pred_values) == n_obs * K) {
    prob_matrix <- matrix(pred_values, nrow = n_obs, ncol = K, byrow = TRUE)
    eval_times <- object$eval_times
  } else {
    stop("Unexpected prediction shape returned by BART.")
  }

  # align internal BART times to requested times via nearest neighbor match
  mapped_times <- sapply(times, function(t) which.min(abs(eval_times - t)))
  survmat <- prob_matrix[, mapped_times, drop = FALSE]
  .finalize_survmat(survmat, times = times)

}

#' Tune BART Survival Hyperparameters (Cross-Validation)
#'
#' Cross-validates BART survival models over a user‐supplied grid and selects the
#' best configuration by the primary metric. Optionally refits the best model on
#' the full dataset.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} with variables referenced in \code{formula}.
#' @param times Numeric vector of evaluation time points.
#' @param param_grid A \code{data.frame} (e.g., from \code{expand.grid()}) of candidate
#'   hyperparameters (e.g., \code{K}, \code{ntree}, \code{power}, \code{base}).
#' @param metrics Character vector of evaluation metrics (default \code{c("cindex","ibs")}).
#'   The first entry is used as the primary selection metric.
#' @param folds Integer; number of cross‐validation folds (default \code{5}).
#' @param seed Integer random seed for reproducibility (default \code{123}).
#' @param ncores Integer number of CPU cores passed to \code{\link{cv_survlearner}}
#'   for fold evaluation (default \code{1}).
#' @param refit_best Logical; if \code{TRUE}, refits the best configuration on the
#'   full data and returns it. If \code{FALSE}, returns a results table.
#'
#' @return If \code{refit_best = FALSE}, a \code{data.frame} (class \code{"tuned_surv"})
#' sorted by the primary metric with one row per grid combination.
#' If \code{refit_best = TRUE}, a fitted \code{mlsurv_model} returned by \code{\link{fit_bart}}.
#'
#' @details
#' Internally calls \code{cv_survlearner()} with \code{fit_bart()}/\code{predict_bart()}
#' so tuning mirrors the production prediction path.
#'
#' @seealso \code{\link{fit_bart}}, \code{\link{predict_bart}}, \code{\link[BART]{surv.bart}}
#'
#' @examplesIf requireNamespace("BART", quietly = TRUE)
#' \donttest{
#' # Evaluate a grid without refitting
#' grid <- expand.grid(K = c(3), ntree = c(50), power = c(2), base = c(0.75, 0.95))
#' res <- tune_bart(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   param_grid = grid,
#'   times = c(10, 60),
#'   refit_best = FALSE
#' )
#' print(res)
#'
#' # Refit best model
#' mod_bart <- tune_bart(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   param_grid = grid,
#'   times = c(10, 60),
#'   refit_best = TRUE
#' )
#' summary(mod_bart)
#' }
#' @export

tune_bart <- function(formula, data, times,
                      param_grid = expand.grid(
                        K = c(3, 5),
                        ntree = c(50, 100),
                        power = c(2, 2.5),
                        base = c(0.75, 0.95)
                      ),
                      metrics = c("cindex", "ibs"),
                      folds = 5, seed = 123,
                      ncores = 1,
                      refit_best = FALSE) {

  results <- purrr::pmap_dfr(param_grid, function(K, ntree, power, base) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_bart,
      pred_fun = predict_bart,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ncores = ncores,
      K = K, ntree = ntree, power = power, base = base
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(K = K, ntree = ntree, power = power, base = base) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric,
                           values_from = mean)
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
    best_model <- fit_bart(
      formula = formula,
      data = data,
      K = best_row$K,
      ntree = best_row$ntree,
      power = best_row$power,
      base = best_row$base
    )
    return(best_model)
  }
}
