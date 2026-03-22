#' Fit a Predictor-Selection Cox Model (pec::selectCox, mlsurv_model-compatible)
#'
#' Fits a Cox proportional hazards model with automated predictor selection
#' using \pkg{pec}'s \code{selectCox()}, returning an `mlsurv_model` object
#' compatible with the `survalis` evaluation and cross‑validation pipeline.
#'
#' @section Engine:
#' Uses \strong{pec::selectCox}. Selection can be driven by Akaike's Information
#' Criterion (\code{rule = "aic"}) or by p‑value thresholds (\code{rule = "p"}),
#' as implemented by \pkg{pec}.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left‑hand side must be a `Surv()` object from the \pkg{survival} package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param rule Selection rule passed to \code{pec::selectCox()}, e.g., \code{"aic"}
#'   (default) or \code{"p"}.
#'
#' @return An object of class `mlsurv_model`, a named list with elements:
#' \describe{
#'   \item{model}{The fitted \code{pec::selectCox} model.}
#'   \item{learner}{Character scalar identifying the learner (`"selectcox"`).}
#'   \item{engine}{Character scalar naming the engine (`"selectCox"`).}
#'   \item{formula}{The original survival formula.}
#'   \item{data}{The training dataset (or a minimal subset needed for prediction).}
#'   \item{rule}{The selection rule used.}
#' }
#'
#' @details
#' This wrapper standardizes the return object to the `mlsurv_model` contract so
#' downstream prediction and evaluation behave consistently across learners.
#'
#' @seealso [predict_selectcox()], [tune_selectcox()], \code{pec::selectCox()}
#'
#' @references
#' Mogensen UB, Ishwaran H, Gerds TA (2012). Evaluating Random Forests for
#' Survival Analysis using \pkg{pec}.
#' Gerds TA, et al. \pkg{pec}: Prediction Error Curves for Survival Models.
#'
#' @examplesIf requireNamespace("pec", quietly = TRUE)
#' mod <- fit_selectcox(Surv(time, status) ~ age + celltype + karno, data = veteran)
#'
#' @export

fit_selectcox <- function(formula, data, rule = "aic") {
  stopifnot(requireNamespace("pec", quietly = TRUE))

  model <- pec::selectCox(formula = formula, data = data, rule = rule)

  structure(list(
    model = model,
    learner = "selectcox",
    formula = formula,
    data = data,
    rule = rule
  ), class = "mlsurv_model", engine = "selectCox")
  }


#' Predict Survival Probabilities with a Selected Cox Model
#'
#' Generates survival probabilities at specified time points using a fitted
#' selected‑predictor Cox model (`mlsurv_model`) via \code{pec::predictSurvProb()}.
#'
#' @param object A fitted `mlsurv_model` returned by [fit_selectcox()].
#' @param newdata A `data.frame` of new observations for prediction. Must include
#'   the same predictor variables used at training time.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   training survival time). Must be non‑negative and finite.
#'
#' @return A base `data.frame` with one row per observation in `newdata` and
#' columns named `t={time}` (character), containing numeric survival
#' probabilities in \[0, 1].
#'
#' @details
#' Internally calls \code{pec::predictSurvProb(object$model, newdata, times)}
#' and renames columns to the standard `t=...` convention used by `survalis`.
#'
#' @seealso [fit_selectcox()], [tune_selectcox()]
#'
#' @examplesIf requireNamespace("pec", quietly = TRUE)
#' mod <- fit_selectcox(Surv(time, status) ~ age + celltype + karno, data = veteran)
#' predict_selectcox(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export

predict_selectcox <- function(object, newdata, times) {


  if (!is.null(object$learner) && object$learner != "selectcox") {
    warning("Object passed to predict_selectcox() may not come from fit_selectcox().")
  }


  stopifnot(requireNamespace("pec", quietly = TRUE))

  pred <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  .finalize_survmat(as.matrix(pred), times = times)

  }


#' Tune SelectCox Rule (Cross-Validation)
#'
#' Cross‑validates \pkg{pec}'s \code{selectCox()} across one or more selection
#' rules and selects the best configuration by the primary metric. Optionally
#' refits the best rule on the full dataset.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left‑hand side must be a `Surv()` object from the \pkg{survival} package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   training survival time). Must be non‑negative and finite.
#' @param rules Character vector of selection rules to compare (e.g.,
#'   \code{c("aic","p")}).
#' @param metrics Character vector of metrics to evaluate/optimize
#'   (e.g., `"cindex"`, `"ibs"`, `"ise"`). The first entry is the primary
#'   selection metric.
#' @param folds Number of cross‑validation folds.
#' @param seed Integer random seed for reproducibility.
#' @param refit_best Logical; if `TRUE`, refits the best rule on the full data
#'   and returns it as an `mlsurv_model` (with tuning results attached).
#' @param ... Additional arguments forwarded to the underlying engine where applicable.
#'
#' @return If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) with a
#' row per rule and metric columns, ordered by the first metric. If
#' `refit_best = TRUE`, a fitted `mlsurv_model` from [fit_selectcox()] with
#' attribute `"tuning_results"` containing the full results.
#'
#' @details
#' \strong{Evaluation.} Internally calls `cv_survlearner()` with
#' `fit_selectcox()`/`predict_selectcox()` to ensure tuning uses the same code
#' paths as production. Rules typically include \code{"aic"} and/or \code{"p"}.
#'
#' @seealso [fit_selectcox()], [predict_selectcox()], \code{pec::selectCox()}
#'
#' @examplesIf requireNamespace("pec", quietly = TRUE)
#' res_selectcox <- tune_selectcox(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   rules = c("aic", "p"),
#'   metrics = c("cindex", "ibs", "ise"),
#'   folds = 3,
#'   refit_best = FALSE
#' )
#'
#' print(res_selectcox)
#' class(res_selectcox)
#'
#' mod_selectcox <- tune_selectcox(
#'   formula = Surv(time, status) ~ age + karno + celltype,
#'   data = veteran,
#'   times = c(100, 200, 300),
#'   rules = c("aic", "p"),
#'   metrics = c("cindex", "ibs", "ise"),
#'   folds = 3,
#'   refit_best = TRUE
#' )
#'
#' summary(mod_selectcox)
#'
#' @export

tune_selectcox <- function(formula, data, times,
                           rules = c("aic", "p"),
                           metrics = c("cindex", "ibs"),
                           folds = 5,
                           seed = 123,
                           refit_best = TRUE,
                           ...) {

  results <- purrr::map_dfr(rules, function(rule) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_selectcox,
      pred_fun = predict_selectcox,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      rule = rule,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(rule = rule) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    results <- dplyr::arrange(results, dplyr::desc(!!sym(metrics[1])))
  } else {
    results <- dplyr::arrange(results, !!sym(metrics[1]))
  }

  if (refit_best) {
    best_rule <- results$rule[1]
    model <- fit_selectcox(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      rule = best_rule,
      ...
    )
    attr(model, "tuning_results") <- results
    return(model)
    }

  class(results) <- c("tuned_surv", class(results))
  return(results)
  }
