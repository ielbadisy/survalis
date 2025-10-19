#' Fit an XGBoost Survival Model (mlsurv_model-compatible)
#'
#' Fits a gradient-boosted tree model for time-to-event outcomes using
#' \pkg{xgboost}. Supports `"survival:aft"` (default) and `"survival:cox"`.
#' Returns an `mlsurv_model` object compatible with the `survalis` pipeline.
#'
#' @section Engine:
#' Uses \strong{xgboost}. For `"survival:aft"`, interval labels are set via
#' \code{label_lower_bound}/\code{label_upper_bound} with right-censoring
#' represented by \code{Inf} on the upper bound. For `"survival:cox"`, labels
#' follow the Cox objective convention.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing variables in \code{formula}.
#' @param booster XGBoost booster type (default \code{"gbtree"}).
#' @param objective One of \code{"survival:aft"} (default) or \code{"survival:cox"}.
#' @param aft_loss_distribution AFT error distribution: \code{"extreme"} (default),
#'   \code{"logistic"}, or \code{"normal"}.
#' @param aft_loss_distribution_scale Positive numeric scale for the AFT loss (default \code{1}).
#' @param nrounds Integer number of boosting iterations (default \code{100}).
#'
#' @return An object of class \code{mlsurv_model}:
#' \describe{
#'   \item{model}{Fitted \code{xgboost} model.}
#'   \item{learner}{\code{"xgboost"}.}
#'   \item{formula, data}{Original inputs.}
#'   \item{time, status}{Names of the survival time and event variables.}
#'   \item{objective}{Objective used.}
#'   \item{dist, scale}{AFT distribution and scale (populated even if Cox is used).}
#'   \item{terms}{Preserved \code{terms} for consistent prediction design matrices.}
#' }
#'
#' @details
#' Design contract: returns a named list with engine metadata and preserves
#' \code{terms(formula)} so prediction uses the exact same encoding as training.
#'
#' @seealso [predict_xgboost()], [tune_xgboost()]
#'
#' @examples
#' \donttest{
#'   mod_xgb <- fit_xgboost(
#'     Surv(time, status) ~ age + karno + celltype,
#'     data = veteran
#'   )
#'   }
#'
#' @export

fit_xgboost <- function(formula, data,
                        booster = "gbtree",
                        objective = "survival:aft",
                        aft_loss_distribution = "extreme",
                        aft_loss_distribution_scale = 1,
                        nrounds = 100) {
  stopifnot(requireNamespace("xgboost", quietly = TRUE))

  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(formula, data = mf)[, -1, drop = FALSE]  # remove intercept

  y_time <- y[, "time"]
  y_event <- y[, "status"] == 1

  dmat <- xgboost::xgb.DMatrix(data = x)
  if (objective == "survival:aft") {
    y_lower <- y_time
    y_upper <- y_time
    y_upper[!y_event] <- Inf
    xgboost::setinfo(dmat, "label_lower_bound", y_lower)
    xgboost::setinfo(dmat, "label_upper_bound", y_upper)
  } else if (objective == "survival:cox") {
    xgboost::setinfo(dmat, "label", y_time)
  } else {
    stop("Unsupported objective.")
  }

  params <- list(
    booster = booster,
    objective = objective,
    aft_loss_distribution = aft_loss_distribution,
    aft_loss_distribution_scale = aft_loss_distribution_scale,
    eval_metric = "cox-nloglik"
  )

  model <- xgboost::xgboost(
    params = params,
    data = dmat,
    nrounds = nrounds,
    verbose = 0
  )

  structure(list(
    model = model,
    learner = "xgboost",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2],
    objective = objective,
    dist = aft_loss_distribution,
    scale = aft_loss_distribution_scale,
    terms = terms(formula)
  ), class = "mlsurv_model", engine = "xgboost")
}

#' Predict Survival with XGBoost
#'
#' Generates predictions from an XGBoost survival \code{mlsurv_model}.
#' If \code{times} is \code{NULL}, returns the raw linear predictor. If \code{times}
#' is provided, returns survival probabilities computed via the AFT mapping
#' using \code{object$dist} and \code{object$scale}.
#'
#' @param object A fitted \code{mlsurv_model} returned by [fit_xgboost()].
#' @param newdata A \code{data.frame} of new observations for prediction.
#' @param times Optional numeric vector of evaluation time points (same scale as training time).
#'   If \code{NULL}, returns the linear predictor.
#'
#' @return
#' \itemize{
#'   \item If \code{times} is \code{NULL}: a numeric vector of linear predictors (one per row of \code{newdata}).
#'   \item If \code{times} is provided: a base \code{data.frame} of survival probabilities with
#'         columns named \code{"t=<time>"} and row names \code{"ID_<i>"}.
#' }
#'
#' @details
#' AFT mapping with \eqn{q = (\log t - \eta)/\sigma}:
#' \itemize{
#'   \item Normal: \eqn{S(t) = 1 - \Phi(q)}
#'   \item Extreme (Gumbel): \eqn{S(t) = \exp\{-\exp(q)\}}
#'   \item Logistic: \eqn{S(t) = 1/(1 + \exp(q))}
#' }
#' Note: If the model was trained with \code{objective = "survival:cox"}, survival
#' probabilities are still computed using the supplied AFT distribution/scale
#' stored in the object; interpret with caution.
#'
#' @seealso [fit_xgboost()], [tune_xgboost()]
#'
#' @examples
#' \donttest{
#'   mod_xgb <- fit_xgboost(Surv(time, status) ~ age + karno + celltype,
#'                          data = veteran, nrounds = 20)
#'   predict_xgboost(mod_xgb, newdata = veteran[1:5, ], times = c(100, 200, 300))
#' }
#'
#' @export

predict_xgboost <- function(object, newdata, times = NULL) {
  stopifnot(object$learner == "xgboost")
  stopifnot(requireNamespace("xgboost", quietly = TRUE))

  mf <- model.frame(object$terms, newdata)
  x_new <- model.matrix(object$terms, data = mf)[, -1, drop = FALSE]

  preds <- predict(object$model, newdata = x_new, outputmargin = TRUE)

  if (is.null(times)) {
    return(preds)  # return linear predictor
  }

  log_times <- outer(rep(1, length(preds)), log(times))
  quants <- (log_times - preds) / object$scale

  surv_probs <- switch(object$dist,
                       "normal" = 1 - pnorm(quants),
                       "extreme" = exp(-exp(quants)),
                       "logistic" = 1 / (1 + exp(quants)),
                       stop("Unsupported distribution.")
  )

  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(x_new)))
  colnames(surv_probs) <- paste0("t=", times)

  as.data.frame(surv_probs)
}


#' Tune XGBoost Survival Hyperparameters (Cross-Validation)
#'
#' Cross-validates XGBoost survival models over a user-specified grid and
#' returns a results table with metric summaries per configuration. Each row
#' that errors during CV is marked \code{failed = TRUE}.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing variables referenced in \code{formula}.
#' @param times Numeric vector of evaluation time points.
#' @param param_grid A \code{data.frame} (e.g., from \code{expand.grid()}) with columns:
#'   \code{nrounds}, \code{max_depth}, \code{eta}, \code{aft_loss_distribution},
#'   \code{aft_loss_distribution_scale}, and \code{objective}. Values are passed
#'   through to [fit_xgboost()] and \code{xgboost::xgboost()}.
#' @param metrics Character vector of metrics to evaluate (e.g., \code{"cindex"}, \code{"ibs"}).
#' @param folds Integer number of CV folds.
#' @param seed Integer random seed for reproducibility.
#'
#' @return A \code{tibble} with one row per grid configuration, containing:
#' \describe{
#'   \item{nrounds, max_depth, eta, aft_loss_distribution, aft_loss_distribution_scale, objective}{Grid values.}
#'   \item{failed}{Logical; \code{TRUE} if the configuration errored.}
#'   \item{<metric columns>}{One column per metric from \code{metrics}, when available.}
#' }
#' The table is arranged by the first metric in \code{metrics} (ascending as implemented).
#'
#' @details
#' Internally calls \code{cv_survlearner()} with \code{fit_xgboost()}/\code{predict_xgboost()}.
#' Any configuration that errors (e.g., due to invalid parameters or data issues) is
#' recorded with \code{failed = TRUE} and omitted from metric summarization.
#'
#' @seealso [fit_xgboost()], [predict_xgboost()]
#'
#' @examples
#' \donttest{
#'   grid <- expand.grid(
#'     nrounds = c(20, 40),
#'     max_depth = c(2, 3),
#'     eta = c(0.1, 0.2),
#'     aft_loss_distribution = c("extreme", "logistic"),
#'     aft_loss_distribution_scale = c(0.5, 1),
#'     objective = "survival:aft",
#'     stringsAsFactors = FALSE
#'   )
#'
#'   res_xgb <- tune_xgboost(
#'     formula = Surv(time, status) ~ age + karno + celltype,
#'     data    = veteran,
#'     times   = c(100, 200),
#'     metrics = c("cindex", "ibs"),
#'     param_grid = grid,
#'     folds   = 2,
#'     seed    = 123
#'   )
#'   head(res_xgb)
#' }
#'
#' @export
tune_xgboost <- function(formula, data, times,
                         param_grid = expand.grid(
                           nrounds = c(50, 100),
                           max_depth = c(3, 6),
                           eta = c(0.01, 0.1),
                           aft_loss_distribution = c("extreme", "logistic"),
                           aft_loss_distribution_scale = c(0.5, 1),
                           objective = "survival:aft",
                           stringsAsFactors = FALSE
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5,
                         seed = 123
                         ) {

  stopifnot(is.numeric(times), all(is.finite(times)), all(times >= 0))

  # allow a named-list grid too
  if (is.list(param_grid) && !is.data.frame(param_grid)) {
    param_grid <- do.call(expand.grid, c(param_grid, stringsAsFactors = FALSE))
  }

  .higher_is_better <- function(metric) metric %in% c("cindex", "auc", "accuracy")

  results <- purrr::pmap_dfr(
    param_grid,
    function(nrounds, max_depth, eta,
             aft_loss_distribution,
             aft_loss_distribution_scale,
             objective) {

      cv_results <- tryCatch({
        cv_survlearner(
          formula  = formula,
          data     = data,
          fit_fun  = fit_xgboost,
          pred_fun = predict_xgboost,
          times    = times,
          metrics  = metrics,
          folds    = folds,
          seed     = seed,
          booster  = "gbtree",
          objective = objective,
          aft_loss_distribution        = aft_loss_distribution,
          aft_loss_distribution_scale  = aft_loss_distribution_scale,
          nrounds  = nrounds,
          max_depth = max_depth,
          eta       = eta
        )
      }, error = function(e) NULL)

      if (is.null(cv_results)) {
        return(tibble::tibble(
          nrounds = nrounds,
          max_depth = max_depth,
          eta = eta,
          aft_loss_distribution = aft_loss_distribution,
          aft_loss_distribution_scale = aft_loss_distribution_scale,
          objective = objective,
          failed = TRUE
        ))
      }

      smry <- cv_summary(cv_results)[, c("metric", "mean")]
      tibble::tibble(
        nrounds = nrounds,
        max_depth = max_depth,
        eta = eta,
        aft_loss_distribution = aft_loss_distribution,
        aft_loss_distribution_scale = aft_loss_distribution_scale,
        objective = objective,
        failed = FALSE
      ) |>
        dplyr::bind_cols(
          tidyr::pivot_wider(
            smry,
            names_from  = "metric",
            values_from = "mean"
          )
        )
    }
  )

  primary <- metrics[1]
  if (primary %in% names(results)) {
    if (.higher_is_better(primary)) {
      results <- dplyr::arrange(results, dplyr::desc(rlang::.data[[primary]]))
    } else {
      results <- dplyr::arrange(results, rlang::.data[[primary]])
    }
  }

  results
}
