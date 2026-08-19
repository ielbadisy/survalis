#' Fit a CoxBoost Survival Model (mlsurv_model-compatible)
#'
#' Fits a likelihood-based boosting model for time-to-event outcomes using
#' \pkg{CoxBoost} (Binder et al. 2009). Returns an \code{mlsurv_model} object
#' compatible with the \code{survalis} pipeline.
#'
#' @section Engine:
#' Uses \strong{CoxBoost}, a componentwise boosting algorithm for the Cox
#' partial likelihood that updates one covariate per boosting step
#' (penalized score statistic selection), unlike gradient-boosted trees.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing variables in \code{formula}.
#' @param stepno Integer number of boosting steps (default \code{100}).
#' @param penalty Numeric penalty applied at each boosting step (default
#'   \code{9 * number of events}, following \code{CoxBoost}'s own default).
#' @param criterion Selection criterion for each boosting step: one of
#'   \code{"pscore"} (default), \code{"score"}, \code{"hpscore"}, \code{"hscore"}.
#' @param stepsize.factor Numeric step-size factor (default \code{1}).
#'
#' @return An object of class \code{mlsurv_model}:
#' \describe{
#'   \item{model}{Fitted \code{CoxBoost} model.}
#'   \item{learner}{\code{"coxboost"}.}
#'   \item{formula, data}{Original inputs.}
#'   \item{time, status}{Names of the survival time and event variables.}
#'   \item{terms}{Preserved \code{terms} for consistent prediction design matrices.}
#' }
#'
#' @details
#' Design contract: returns a named list with engine metadata and preserves
#' \code{terms(formula)} so prediction uses the exact same encoding as training.
#'
#' @seealso [predict_coxboost()], [tune_coxboost()]
#'
#' @examples
#' \donttest{
#'   mod_cb <- fit_coxboost(
#'     Surv(time, status) ~ age + karno + celltype,
#'     data = veteran,
#'     stepno = 50
#'   )
#'   }
#'
#' @export

fit_coxboost <- function(formula, data,
                         stepno = 100,
                         penalty = NULL,
                         criterion = "pscore",
                         stepsize.factor = 1) {
  stopifnot(requireNamespace("CoxBoost", quietly = TRUE))

  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(formula, data = mf)[, -1, drop = FALSE]  # remove intercept

  y_time <- y[, "time"]
  y_status <- as.integer(y[, "status"] == 1)

  if (is.null(penalty)) {
    penalty <- 9 * sum(y_status == 1)
  }

  model <- CoxBoost::CoxBoost(
    time = y_time,
    status = y_status,
    x = x,
    stepno = stepno,
    penalty = penalty,
    criterion = criterion,
    stepsize.factor = stepsize.factor
  )

  structure(list(
    model = model,
    learner = "coxboost",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2],
    terms = terms(formula)
  ), class = "mlsurv_model", engine = "coxboost")
}

#' Predict Survival with CoxBoost
#'
#' Generates survival-probability predictions from a CoxBoost
#' \code{mlsurv_model}.
#'
#' @param object A fitted \code{mlsurv_model} returned by [fit_coxboost()].
#' @param newdata A \code{data.frame} of new observations for prediction.
#' @param times Optional numeric vector of evaluation time points. If
#'   \code{NULL}, returns the linear predictor.
#'
#' @return
#' \itemize{
#'   \item If \code{times} is \code{NULL}: a numeric vector of linear predictors
#'         (one per row of \code{newdata}).
#'   \item If \code{times} is provided: a base \code{data.frame} of survival
#'         probabilities with columns named \code{"t=<time>"}.
#' }
#'
#' @seealso [fit_coxboost()], [tune_coxboost()]
#'
#' @examples
#' \donttest{
#'   mod_cb <- fit_coxboost(Surv(time, status) ~ age + karno + celltype,
#'                          data = veteran, stepno = 30)
#'   predict_coxboost(mod_cb, newdata = veteran[1:5, ], times = c(100, 200, 300))
#' }
#'
#' @export

predict_coxboost <- function(object, newdata, times = NULL) {
  stopifnot(object$learner == "coxboost")
  stopifnot(requireNamespace("CoxBoost", quietly = TRUE))

  mf <- model.frame(object$terms, newdata)
  x_new <- model.matrix(object$terms, data = mf)[, -1, drop = FALSE]

  if (is.null(times)) {
    lp <- predict(object$model, newdata = x_new, type = "lp")
    return(as.numeric(lp))
  }

  surv_probs <- predict(object$model, newdata = x_new, times = times, type = "risk")

  .finalize_survmat(surv_probs, times = times)
}

#' Tune CoxBoost Survival Hyperparameters (Cross-Validation)
#'
#' Cross-validates CoxBoost survival models over a user-specified grid and
#' returns a results table with metric summaries per configuration. Any row
#' that errors during CV is marked \code{failed = TRUE}.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing variables referenced in \code{formula}.
#' @param times Numeric vector of evaluation time points.
#' @param param_grid A \code{data.frame} (e.g., from \code{expand.grid()}) with
#'   columns: \code{stepno}, \code{penalty}, \code{criterion}, \code{stepsize.factor}.
#'   Values are passed through to \code{\link{fit_coxboost}}.
#' @param metrics Character vector of metrics to evaluate (e.g., \code{"cindex"}, \code{"ibs"}).
#'   The first entry is treated as the primary selection metric for ordering.
#' @param folds Integer number of CV folds.
#' @param seed Integer random seed for reproducibility.
#' @param ncores Integer number of CPU cores passed to \code{\link{cv_survlearner}}
#'   for fold evaluation (default \code{1}).
#'
#' @return A \code{data.table} with one row per grid configuration, containing
#'   the grid values, a \code{failed} column, and one column per metric.
#'
#' @details
#' Internally calls \code{\link{cv_survlearner}} with \code{\link{fit_coxboost}} /
#' \code{\link{predict_coxboost}}. Any configuration that errors is recorded
#' with \code{failed = TRUE} and omitted from metric summarization.
#'
#' @seealso \code{\link{fit_coxboost}}, \code{\link{predict_coxboost}}
#'
#' @examples
#' \donttest{
#'   grid <- expand.grid(
#'     stepno = c(50, 100),
#'     penalty = c(100, 200),
#'     criterion = "pscore",
#'     stepsize.factor = 1,
#'     stringsAsFactors = FALSE
#'   )
#'
#'   res_cb <- tune_coxboost(
#'     formula = survival::Surv(time, status) ~ age + karno + celltype,
#'     data    = survival::veteran,
#'     times   = c(100, 200),
#'     metrics = c("cindex", "ibs"),
#'     param_grid = grid,
#'     folds   = 2,
#'     seed    = 123
#'   )
#'   head(res_cb)
#' }
#'
#' @export

tune_coxboost <- function(formula, data, times,
                          param_grid = expand.grid(
                            stepno = c(50, 100),
                            penalty = c(100, 200),
                            criterion = "pscore",
                            stepsize.factor = 1,
                            stringsAsFactors = FALSE
                          ),
                          metrics = c("cindex", "ibs"),
                          folds = 5,
                          seed = 123,
                          ncores = 1
                          ) {

  stopifnot(is.numeric(times), all(is.finite(times)), all(times >= 0))

  if (is.list(param_grid) && !is.data.frame(param_grid)) {
    param_grid <- do.call(expand.grid, c(param_grid, stringsAsFactors = FALSE))
  }

  .higher_is_better <- function(metric) metric %in% c("cindex", "auc", "accuracy")

  results <- .pmap_rbind_dt(
    param_grid,
    function(stepno, penalty, criterion, stepsize.factor) {

      params <- list(
        stepno = stepno,
        penalty = penalty,
        criterion = criterion,
        stepsize.factor = stepsize.factor
      )

      cv_results <- tryCatch({
        cv_survlearner(
          formula  = formula,
          data     = data,
          fit_fun  = fit_coxboost,
          pred_fun = predict_coxboost,
          times    = times,
          metrics  = metrics,
          folds    = folds,
          seed     = seed,
          ncores   = ncores,
          stepno   = stepno,
          penalty  = penalty,
          criterion = criterion,
          stepsize.factor = stepsize.factor
        )
      }, error = function(e) NULL)

      if (is.null(cv_results)) {
        return(do.call(data.table::data.table, c(params, list(failed = TRUE))))
      }

      smry <- cv_summary(cv_results)[, c("metric", "mean")]
      .wide_metric_row(c(params, list(failed = FALSE)), smry)
    }
  )

  primary <- metrics[1]
  if (primary %in% names(results)) {
    results <- .arrange_by_metric_dt(results, primary, .higher_is_better(primary))
  }

  results
}
