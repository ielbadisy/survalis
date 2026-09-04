#' Superseded functions in survalis 2.x
#'
#' survalis 2.0 introduced the seven-verb interface -- [fit()], [predict()],
#' [evaluate()], [tune()], [compare()], [interpret()], [estimate()] -- plus the
#' resampling constructors ([cv()], [holdout()], [group_cv()], [bootstrap()]).
#'
#' The granular functions listed below are the machinery those verbs are built
#' on. They remain **exported and fully supported for the whole 2.x series**, are
#' marked internal in the manual, and emit a one-time note per session. They are
#' scheduled to become internal (unexported) in survalis 3.0.
#'
#' @section Mapping:
#' \tabular{ll}{
#'   \strong{Superseded} \tab \strong{Use instead} \cr
#'   `fit_<learner>()` \tab `fit(formula, data, model = "<learner>")` \cr
#'   `predict_<learner>()` \tab `predict(fit, newdata, times, type)` \cr
#'   `tune_<learner>()` \tab `tune(formula, data, model = "<learner>", grid)` \cr
#'   `cv_survlearner()`, `cv_summary()`, `cv_plot()` \tab `evaluate(...)` \cr
#'   `score_survmodel()` \tab `evaluate(fit, resampling = ...)` \cr
#'   `benchmark()`, `benchmark_default_survlearners()`,
#'     `benchmark_tuned_survlearners()` \tab `compare(...)` \cr
#'   `summarise_benchmark()`, `summarize_benchmark_results()`,
#'     `best_survlearner()` \tab `summary()` / `print()` on a `survalis_compare` \cr
#'   `compute_*()` / `plot_*()` \tab `interpret(fit, method = ...)` then `plot()` \cr
#'   `estimate_rmst()` \tab `estimate(..., estimand = "RMST_diff")` \cr
#' }
#'
#' @name survalis-deprecated
#' @keywords internal
NULL

#' One-time session note that a legacy entry point was called directly
#'
#' A no-op while a verb is on the call stack (the verbs reuse this machinery).
#' @keywords internal
.legacy_notice <- function() {
  if (isTRUE(getOption("survalis.quiet_legacy", FALSE))) return(invisible())
  rlang::warn(
    paste0(
      "You are calling a granular survalis function directly. The seven-verb ",
      "interface (fit/predict/evaluate/tune/compare/interpret/estimate) is now ",
      "recommended; these helpers become internal in survalis 3.0. ",
      "See ?`survalis-deprecated`."
    ),
    .frequency = "once", .frequency_id = "survalis-legacy-api",
    class = "survalis_legacy_notice"
  )
  invisible()
}

#' Run an expression with the legacy note suppressed (used by every verb)
#' @keywords internal
.with_quiet_legacy <- function(expr) {
  old <- options(survalis.quiet_legacy = TRUE)
  on.exit(options(old), add = TRUE)
  force(expr)
}
