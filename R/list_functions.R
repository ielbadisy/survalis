# R/catalogs.R

# internal helper
.exists_in_pkg <- function(fn) {
  !is.null(get0(fn, envir = asNamespace("survalis"), inherits = FALSE))
}

#' List survival learners available in survalis
#'
#' Returns a table of known survival learners, showing the expected
#' `fit_*`, `predict_*`, and (if present) `tune_*` functions, plus
#' booleans indicating which functions are available.
#'
#' @param has_tune Logical (default `FALSE`). If `TRUE`, return only learners
#'   that have a corresponding `tune_*` function. If `FALSE`, return all.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item \code{learner} – learner id (e.g., "ranger")
#'     \item \code{fit}, \code{predict}, \code{tune} – function names (tune may be \code{NA})
#'     \item \code{has_fit}, \code{has_predict}, \code{has_tune} – logical flags
#'     \item \code{available} – \code{has_fit & has_predict}
#'   }
#'
#' @examples
#' list_survlearners()                   # all learners (default)
#' list_survlearners(has_tune = TRUE)    # only tunable learners
#' @export
list_survlearners <- function(has_tune = FALSE) {
  # curated order
  ids <- c(
    "coxph", "aalen", "glmnet", "selectcox", "aftgee", "flexsurvreg", "stpm2", "bnnsurv",
    "rpart", "bart", "xgboost", "ranger", "rsf", "cforest", "blackboost", "survsvm",
    "survdnn", "orsf", "survmetalearner"
  )

  fun_exists <- function(x) isTRUE(exists(x, mode = "function", inherits = TRUE))

  tbl <- tibble::tibble(
    learner = ids,
    fit     = paste0("fit_", ids),
    predict = paste0("predict_", ids),
    tune    = paste0("tune_", ids)
  )

  tbl$has_fit     <- vapply(tbl$fit,     fun_exists, logical(1))
  tbl$has_predict <- vapply(tbl$predict, fun_exists, logical(1))
  tbl$has_tune    <- vapply(tbl$tune,    fun_exists, logical(1))
  tbl$tune[!tbl$has_tune] <- NA_character_
  tbl$available <- tbl$has_fit & tbl$has_predict

  if (isTRUE(has_tune)) {
    tbl <- tbl[tbl$has_tune, , drop = FALSE]
  }

  tbl
}


#' List interpretability methods available in survalis
#'
#' Returns a table mapping each \code{compute_*} function to its paired
#' \code{plot_*} helper (if any). Methods without a plot helper show \code{NA}.
#'
#' @return A tibble with columns \code{compute}, \code{plot},
#'   \code{has_compute}, and \code{has_plot}.
#' @examples
#' list_interpretability_methods()
#' subset(list_interpretability_methods(), is.na(plot))  # methods without plot
#' @export
list_interpretability_methods <- function() {
  compute_candidates <- c(
    "compute_shap",
    "compute_pdp",
    "compute_ale",
    "compute_surrogate",
    "compute_tree_surrogate",
    "compute_varimp",
    "compute_interactions",
    "compute_counterfactual"  # no plot helper (yet)
  )

  plot_candidates <- c(
    "plot_shap",
    "plot_pdp",
    "plot_ale",
    "plot_surrogate",
    "plot_tree_surrogate",
    "plot_varimp",
    "plot_interactions",
    NA_character_              # aligns with compute_counterfactual
  )

  stopifnot(length(compute_candidates) == length(plot_candidates))

  fun_exists <- function(x) {
    if (is.na(x)) return(FALSE)
    isTRUE(exists(x, mode = "function", inherits = TRUE))
  }

  tibble::tibble(
    compute     = compute_candidates,
    plot        = plot_candidates,
    has_compute = vapply(compute_candidates, fun_exists, logical(1)),
    has_plot    = vapply(plot_candidates,    fun_exists, logical(1))
  )
}



#' List tunable survival learners
#'
#' Filters learners to only those that provide a `tune_*` function in addition
#' to `fit_*` and `predict_*`.
#'
#' @return A base `data.frame` like `list_survlearners()` but containing only
#' rows where `tune` is not `NA`.
#'
#' @examples
#' list_tunable_survlearners()
#' @export
list_tunable_survlearners <- function() {
  all <- list_survlearners()
  if (is.null(all) || nrow(all) == 0L) return(all)
  all[!is.na(all$tune), , drop = FALSE]
}



#' List Available Evaluation Metrics
#'
#' Returns a data frame describing the evaluation metrics supported by
#' `survalis` (as used by helpers like `cv_survlearner()` and `score_survmodel()`).
#'
#' Columns:
#' - `metric`: short metric name used throughout the package.
#' - `direction`: whether higher or lower values are better.
#' - `summary`: brief description.
#' - `range`: typical value range.
#'
#' @return A tibble/data.frame with one row per metric.
#' @examples
#' list_metrics()
#' @export
list_metrics <- function() {
  tibble::tibble(
    metric    = c("cindex", "brier", "ibs"),
    direction = c("maximize", "minimize", "minimize"),
    summary   = c(
      "Harrell-style concordance index for survival predictions.",
      "Brier Score at specified evaluation time(s) (IPCW-weighted when needed).",
      "Integrated Brier Score over an evaluation time grid (IPCW-weighted)."
    ),
    range     = c(
      "[0, 1] (higher is better)",
      "[0, 1] (lower is better)",
      "[0, 1] (lower is better)"
    )
  )
}
