#' Low-level functions in survalis
#'
#' survalis 1.3.0 introduced the seven-verb interface -- [fit()], [predict()],
#' [evaluate()], [tune()], [compare()], [interpret()], [estimate()] -- plus the
#' resampling constructors ([cv()], [holdout()], [group_cv()], [bootstrap()]).
#'
#' The granular functions listed below are the machinery those verbs are built
#' on. They stay exported and supported for direct use, for example in
#' simulation code that calls one learner many times.
#'
#' @section Mapping:
#' \tabular{ll}{
#'   \strong{Low-level} \tab \strong{Verb} \cr
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
#' @name survalis-lowlevel
#' @keywords internal
NULL
