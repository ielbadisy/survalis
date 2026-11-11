# R/catalogs.R

# internal helper
.exists_in_pkg <- function(fn) {
  !is.null(get0(fn, envir = asNamespace("survalis"), inherits = FALSE))
}

#' List the survival learners
#'
#' The ids accepted by the `model` argument of [fit()], [evaluate()],
#' [compare()], [tune()] and [estimate()].
#'
#' @param has_tune Logical (default `FALSE`). If `TRUE`, return only learners
#'   with a built-in tuning grid, so that [tune()] works without a `grid`.
#'
#' @return A data frame with columns `learner` (the id) and `tunable` (whether
#'   a built-in tuning grid exists).
#'
#' @examples
#' listlearners()
#' listlearners(has_tune = TRUE)
#' @export
listlearners <- function(has_tune = FALSE) {
  ids <- c(
    "coxph", "aalen", "glmnet", "selectcox", "aftgee", "flexsurvreg", "bnnsurv",
    "rpart", "bart", "xgboost", "coxboost", "fastgbm", "ranger", "rsf", "cforest", "blackboost",
    "survsvm", "densemlp", "orsf", "rpsurv", "nbsurv", "survmetalearner"
  )
  tunable <- vapply(paste0("tune_", ids), .exists_in_pkg, logical(1))
  tbl <- data.frame(learner = ids, tunable = unname(tunable),
                    stringsAsFactors = FALSE)
  if (isTRUE(has_tune)) tbl <- tbl[tbl$tunable, , drop = FALSE]
  rownames(tbl) <- NULL
  tbl
}


#' List the interpretability methods
#'
#' The ids accepted by the `method` argument of [interpret()].
#'
#' @return A data frame with columns `method` and `summary`.
#' @examples
#' listmethods()
#' @export
listmethods <- function() {
  data.frame(
    method = c("pdp", "ice", "ale", "shap", "shap_matrix", "varimp", "permute",
               "interaction", "counterfactual", "surrogate", "tree_surrogate",
               "calibration", "roc", "dca", "survcluster"),
    summary = c(
      "Partial dependence of predicted survival on one feature, over time.",
      "Individual conditional expectation curves for one feature.",
      "Accumulated local effects for one feature.",
      "Shapley values for one instance, aggregated over time.",
      "Shapley values for many instances (beeswarm plot).",
      "Permutation variable importance on the survival scale.",
      "Alias of varimp.",
      "Friedman H-statistic interaction strength per feature.",
      "Smallest feature changes that move predicted survival at a target time.",
      "Local linear surrogate for one instance at a target time.",
      "Global decision-tree surrogate of the predictions.",
      "Calibration curve at one time.",
      "Time-dependent ROC curve at one time.",
      "Decision curve (net benefit over risk thresholds) at one time.",
      "Clusters of subjects with similar predicted survival curves."
    ),
    stringsAsFactors = FALSE
  )
}



#' List the evaluation metrics
#'
#' The ids accepted by the `metrics` argument of [evaluate()] and [compare()]
#' and the `metric` argument of [tune()].
#'
#' Columns:
#' - `metric`: short metric name used throughout the package.
#' - `direction`: whether higher or lower values are better.
#' - `summary`: brief description.
#' - `range`: typical value range.
#'
#' `"range:<lower>-<upper>"` is a template: write the thresholds into the id,
#' for example `"range:0.05-0.20"`.
#'
#' @return A data.frame with one row per metric.
#' @examples
#' listmetrics()
#' @export
listmetrics <- function() {
  data.frame(
    metric    = c("cindex", "auc", "brier", "ibs", "iae", "ise", "ece",
                  "survacc", "range:<lower>-<upper>"),
    direction = c("maximize", "maximize", "minimize", "minimize", "minimize", "minimize", "minimize",
                  "maximize", "minimize"),
    summary   = c(
      "Harrell-style concordance index for survival predictions.",
      "Cumulative/dynamic time-dependent AUC at a selected evaluation time.",
      "Brier Score at specified evaluation time(s) (IPCW-weighted when needed).",
      "Integrated Brier Score over an evaluation time grid (IPCW-weighted).",
      "Integrated absolute error against the Kaplan-Meier curve.",
      "Integrated squared error against the Kaplan-Meier curve.",
      "Expected calibration error at a single evaluation time.",
      "SurvAcc: spherical decision-weighted score at a single time (IPCW), see dwscore().",
      "Decision-weighted loss over a risk-threshold range, e.g. \"range:0.05-0.20\", at a single time, see dwscore()."
    ),
    range     = c(
      "[0, 1] (higher is better)",
      "[0, 1] (higher is better)",
      "[0, 1] (lower is better)",
      "[0, 1] (lower is better)",
      "[0, Inf) (lower is better)",
      "[0, Inf) (lower is better)",
      "[0, 1] (lower is better)",
      "[0, 1] (higher is better)",
      "[0, 1] (lower is better)"
    ),
    stringsAsFactors = FALSE
  )
}
