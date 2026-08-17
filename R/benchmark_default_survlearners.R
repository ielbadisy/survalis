
#' Benchmark Survival Learners (Meta Function)
#'
#' The single entry point for comparing multiple survival learners. With
#' \code{tune = FALSE} (default), each learner is run with fixed
#' hyperparameters via \code{\link{benchmark_default_survlearners}}. With
#' \code{tune = TRUE}, each learner is tuned internally per outer fold via
#' nested cross-validation using \code{\link{benchmark_tuned_survlearners}}.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ x1 + x2 + ...}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param learners Character vector of learner ids (without prefixes), e.g.
#'   \code{c("ranger", "coxph", "glmnet")}.
#' @param times Numeric vector of evaluation time points for survival predictions.
#' @param metrics Character vector of metrics to compute/optimize.
#' @param tune Logical; if \code{FALSE} (default), benchmarks learners with
#'   fixed hyperparameters (\code{folds}, \code{ncores} apply). If \code{TRUE},
#'   tunes each learner via nested CV (\code{outer_folds}, \code{inner_folds},
#'   \code{inner_ncores}, \code{learner_args}, \code{refit_final} apply).
#' @param folds Integer number of CV folds; used only when \code{tune = FALSE}.
#' @param outer_folds,inner_folds Integer numbers of outer/inner CV folds; used
#'   only when \code{tune = TRUE}.
#' @param seed Integer random seed.
#' @param ncores Integer CPU cores for fold-level mapping; used only when
#'   \code{tune = FALSE}.
#' @param inner_ncores Integer CPU cores for each learner's inner tuning CV;
#'   used only when \code{tune = TRUE}.
#' @param learner_args Optional named list of learner-specific tuning
#'   arguments; used only when \code{tune = TRUE}. See
#'   \code{\link{benchmark_tuned_survlearners}}.
#' @param refit_final Logical; if \code{TRUE} (and \code{tune = TRUE}), tunes
#'   and refits each learner on the full dataset in addition to the outer CV.
#' @param verbose Logical; print progress.
#' @param suppress_errors Logical; if \code{TRUE} (default), per-learner
#'   errors are caught and reported via \code{warning()} instead of stopping
#'   the whole benchmark.
#' @param ... Additional arguments forwarded to \code{benchmark_default_survlearners()}
#'   (fit arguments) when \code{tune = FALSE}, or to \code{benchmark_tuned_survlearners()}
#'   (common tuning arguments) when \code{tune = TRUE}.
#'
#' @return With \code{tune = FALSE}, a data frame of CV results (see
#'   \code{\link{benchmark_default_survlearners}}). With \code{tune = TRUE}, a
#'   list of class \code{"nested_surv_benchmark"} (see
#'   \code{\link{benchmark_tuned_survlearners}}).
#'
#' @examples
#' \donttest{
#' res <- benchmark(
#'   Surv(time, status) ~ age + karno + trt,
#'   data = veteran,
#'   learners = c("coxph", "rpart"),
#'   times = c(80, 160),
#'   metrics = c("cindex", "ibs"),
#'   folds = 2,
#'   seed = 1
#' )
#' head(res)
#' }
#'
#' @seealso [benchmark_default_survlearners()], [benchmark_tuned_survlearners()]
#' @export
benchmark <- function(formula, data, learners, times,
                      metrics = c("cindex", "ibs"),
                      tune = FALSE,
                      folds = 5,
                      outer_folds = 5,
                      inner_folds = 5,
                      seed = 123,
                      ncores = 1,
                      inner_ncores = 1,
                      learner_args = list(),
                      refit_final = FALSE,
                      verbose = FALSE,
                      suppress_errors = TRUE,
                      ...) {
  if (isTRUE(tune)) {
    benchmark_tuned_survlearners(
      formula = formula, data = data, learners = learners, times = times,
      metrics = metrics, outer_folds = outer_folds, inner_folds = inner_folds,
      seed = seed, inner_ncores = inner_ncores, learner_args = learner_args,
      refit_final = refit_final, verbose = verbose,
      suppress_errors = suppress_errors, ...
    )
  } else {
    benchmark_default_survlearners(
      formula = formula, data = data, learners = learners, times = times,
      metrics = metrics, folds = folds, seed = seed, ncores = ncores,
      verbose = verbose, suppress_errors = suppress_errors, ...
    )
  }
}

#' Benchmark Multiple Survival Learners (Cross-Validation Wrapper)
#'
#' Runs \code{cv_survlearner()} for a set of learner names (e.g., \code{"ranger"},
#' \code{"coxph"}) by dynamically dispatching \code{fit_<learner>} and
#' \code{predict_<learner>} functions. Returns the row‑bound CV results across
#' all requested learners.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ x1 + x2 + ...}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param learners Character vector of learner ids (without prefixes), e.g.
#'   \code{c("ranger","coxph","glmnet")}. For each id, \code{fit_<id>} and
#'   \code{predict_<id>} must exist.
#' @param times Numeric vector of evaluation time points for survival predictions.
#' @param metrics Character vector of metrics to compute in CV (e.g., \code{"cindex"},
#'   \code{"ibs"}, \code{"iae"}, \code{"ise"}). The metric semantics are those of
#'   \code{cv_survlearner()}.
#' @param folds Integer number of CV folds (default \code{5}).
#' @param seed Integer random seed for fold generation.
#' @param ncores Integer; number of CPU cores passed to \code{cv_survlearner()}
#'   for fold-level parallel mapping (default \code{1}).
#' @param verbose Logical; if \code{TRUE}, prints progress per learner.
#' @param suppress_errors Logical; if \code{TRUE} (default) errors from individual
#'   learners are caught and reported via \code{warning()} and benchmarking continues.
#'   If \code{FALSE}, the first error is thrown.
#' @param ... Additional arguments forwarded to each learner's \code{fit_*()}.
#'
#' @details
#' Learners are run independently using identical CV splits and scoring settings.
#' Any learner whose \code{fit_*()} or \code{predict_*()} function is missing will
#' be skipped with a warning. At least one learner must complete successfully or an
#' error is raised.
#'
#' @return A data frame of CV results (as returned by \code{cv_survlearner()})
#' with an extra column \code{learner} identifying the source learner.
#'
#' @examples
#' \donttest{
#' res <- benchmark_default_survlearners(
#'   Surv(time, status) ~ age + karno + trt,
#'   data = veteran,
#'   learners = c("coxph", "rpart"),
#'   times = c(80, 160),
#'   metrics = c("cindex", "ibs"),
#'   folds = 2,
#'   seed = 1
#' )
#' head(res)
#' }
#'
#' @seealso [cv_survlearner()], [plot_benchmark()], [summarise_benchmark()]
#' @export

benchmark_default_survlearners <- function(formula, data, learners, times,
                                   metrics = c("cindex", "ibs"),
                                   folds = 5, seed = 123,
                                   ncores = 1,
                                   verbose = FALSE,
                                   suppress_errors = TRUE,
                                   ...) {
  stopifnot(is.character(learners), length(learners) > 0)

  results_list <- lapply(learners, function(learner) {
    fit_fun_name  <- paste0("fit_", learner)
    pred_fun_name <- paste0("predict_", learner)

    if (!exists(fit_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", fit_fun_name))
      return(NULL)
    }
    if (!exists(pred_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", pred_fun_name))
      return(NULL)
    }

    fit_fun  <- get(fit_fun_name)
    pred_fun <- get(pred_fun_name)

    if (verbose) message("Benchmarking learner: ", learner)

    result <- try(
      cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_fun,
        pred_fun = pred_fun,
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        ncores = ncores,
        verbose = verbose,
        ...
      ),
      silent = suppress_errors
    )

    if (inherits(result, "try-error")) {
      if (!suppress_errors) stop(result)
      warning(sprintf("Learner '%s' failed: %s", learner, conditionMessage(attr(result, "condition"))))
      return(NULL)
    }

    result$learner <- learner
    result
  })

  # combine results
  results_combined <- data.table::rbindlist(results_list, fill = TRUE)
  if (nrow(results_combined) == 0) {
    stop("All learners failed or returned empty results.")
  }

  results_combined
}


.nested_surv_prepare_data <- function(formula, data, verbose = FALSE) {
  parsed_formula <- .parse_surv_formula(formula, data)
  status_col <- parsed_formula$status_col
  event_value <- parsed_formula$event_value
  recode_status <- parsed_formula$recode_status

  all_vars <- all.vars(formula)
  n_before <- nrow(data)
  data <- .complete_cases_df(data, all_vars)
  n_after <- nrow(data)

  if (verbose && n_after < n_before) {
    message("Dropped ", n_before - n_after, " rows with missing values (from ", n_before, " to ", n_after, ").")
  }

  if (recode_status) {
    data[[status_col]] <- as.integer(data[[status_col]] == event_value)
  }

  list(data = data, parsed_formula = parsed_formula)
}


.nested_surv_score_predictions <- function(formula, data, pred, times, metrics) {
  parsed_formula <- .parse_surv_formula(formula, data)
  time_col <- parsed_formula$time_col
  status_col <- parsed_formula$status_col
  event_value <- parsed_formula$event_value

  if (parsed_formula$recode_status) {
    status_vector <- as.integer(data[[status_col]] == event_value)
  } else {
    status_vector <- data[[status_col]]
  }

  surv_obj <- survival::Surv(time = data[[time_col]], event = status_vector)

  .score_metrics(surv_obj, pred, times, metrics)
}


.nested_surv_higher_is_better <- function(metric) {
  !metric %in% c("ibs", "brier", "iae", "ise", "ece")
}


.nested_surv_param_cols <- function(tuning_results) {
  metric_cols <- c("cindex", "auc", "brier", "ibs", "iae", "ise", "ece")
  setdiff(names(tuning_results), c(metric_cols, "failed"))
}


.nested_surv_row_to_list <- function(row_df) {
  if (!is.data.frame(row_df) || nrow(row_df) != 1L) {
    stop("Expected a one-row data frame of tuning parameters.")
  }

  out <- as.list(row_df)
  for (nm in names(out)) {
    if (is.list(out[[nm]]) && length(out[[nm]]) == 1L) {
      out[[nm]] <- out[[nm]][[1]]
    }
  }
  out
}


.nested_surv_sort_tuning_results <- function(tuning_results, primary_metric) {
  if (!primary_metric %in% names(tuning_results)) {
    stop("Primary metric column not found in tuning results: ", primary_metric)
  }

  tuning_results <- tuning_results[!is.na(tuning_results[[primary_metric]]), , drop = FALSE]
  if (nrow(tuning_results) == 0L) {
    stop("All tuning configurations failed or returned NA for the primary metric.")
  }

  .arrange_by_metric_dt(tuning_results, primary_metric, .nested_surv_higher_is_better(primary_metric))
}


.nested_surv_tune_once <- function(tune_fun, formula, data, times, metrics,
                                   inner_folds, seed, inner_ncores = 1,
                                   tune_args = list()) {
  args <- c(
    list(
      formula = formula,
      data = data,
      times = times,
      metrics = metrics,
      folds = inner_folds,
      seed = seed,
      ncores = inner_ncores
    ),
    tune_args
  )

  if ("refit_best" %in% names(formals(tune_fun))) {
    args$refit_best <- FALSE
  }

  tuning_results <- do.call(tune_fun, args)

  if (!is.data.frame(tuning_results) || nrow(tuning_results) == 0L) {
    stop("Tuning did not return any results.")
  }

  if ("failed" %in% names(tuning_results)) {
    failed <- tuning_results$failed
    failed[is.na(failed)] <- FALSE
    tuning_results <- tuning_results[!failed, , drop = FALSE]
    if (nrow(tuning_results) == 0L) {
      stop("All tuning configurations failed.")
    }
  }

  .nested_surv_sort_tuning_results(tuning_results, metrics[1])
}


#' Benchmark Tuned Survival Learners with Nested Cross-Validation
#'
#' Runs nested cross-validation for one or more learners that expose
#' `tune_<learner>()`, `fit_<learner>()`, and `predict_<learner>()`. The outer
#' folds estimate performance, while each inner training split is tuned using the
#' learner's existing `tune_*()` implementation only on the outer training data.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ x1 + x2 + ...`.
#' @param data A data frame containing the variables in `formula`.
#' @param learners Character vector of learner ids (without prefixes), e.g.
#'   `c("ranger", "glmnet")`.
#' @param times Numeric vector of evaluation time points for survival predictions.
#' @param metrics Character vector of metrics to optimize and report.
#' @param outer_folds Integer number of outer CV folds used for performance estimation.
#' @param inner_folds Integer number of inner CV folds used by each learner's
#'   tuning routine.
#' @param seed Integer random seed used for outer and inner resampling.
#' @param inner_ncores Integer; number of CPU cores passed to each learner's
#'   inner \code{cv_survlearner()} calls during tuning (default \code{1}).
#' @param learner_args Optional named list of learner-specific arguments. Each
#'   entry can be either a list of tuning arguments passed to `tune_<learner>()`,
#'   or a list with `tune` and/or `fit` components. The `fit` component is
#'   forwarded only when refitting the selected hyperparameters on an outer
#'   training split or on the full dataset.
#' @param refit_final Logical; if `TRUE`, tunes each learner on the full dataset
#'   and refits a final model with the selected hyperparameters.
#' @param verbose Logical; if `TRUE`, prints progress per learner and outer fold.
#' @param suppress_errors Logical; if `TRUE` (default), errors from individual
#'   learners or outer folds are caught and reported via `warning()`.
#' @param ... Additional arguments passed to each learner's `tune_*()` function.
#'
#' @return A list of class `"nested_surv_benchmark"` with components:
#' \describe{
#'   \item{outer_results}{Per-outer-fold metric values with columns `learner`,
#'   `id`, `outer_fold`, `metric`, and `value`.}
#'   \item{outer_summary}{Mean and dispersion summary of `outer_results` by learner
#'   and metric.}
#'   \item{selected_params}{One row per learner and outer fold containing the
#'   selected hyperparameters.}
#'   \item{final_models}{Named list of final refit models when `refit_final = TRUE`,
#'   otherwise `NULL`.}
#'   \item{settings}{List of benchmark settings.}
#' }
#'
#' @examples
#' \donttest{
#' res <- benchmark_tuned_survlearners(
#'   Surv(time, status) ~ age + karno + trt,
#'   data = veteran,
#'   learners = c("ranger", "glmnet"),
#'   times = c(100, 200),
#'   outer_folds = 3,
#'   inner_folds = 2,
#'   learner_args = list(
#'     ranger = list(tune = list(
#'       param_grid = expand.grid(
#'         num.trees = c(100, 200),
#'         mtry = c(1, 2),
#'         min.node.size = c(3, 5)
#'       )
#'     ))
#'   )
#' )
#' res$outer_summary
#' }
#'
#' @seealso \code{\link{benchmark_default_survlearners}},
#'   \code{\link{cv_survlearner}}
#' @export
benchmark_tuned_survlearners <- function(formula, data, learners, times,
                                         metrics = c("cindex", "ibs"),
                                         outer_folds = 5,
                                         inner_folds = 5,
                                         seed = 123,
                                         inner_ncores = 1,
                                         learner_args = list(),
                                         refit_final = FALSE,
                                         verbose = FALSE,
                                         suppress_errors = TRUE,
                                         ...) {
  stopifnot(is.character(learners), length(learners) > 0)
  stopifnot(is.list(learner_args))

  prepared <- .nested_surv_prepare_data(formula, data, verbose = verbose)
  data <- prepared$data
  status_col <- prepared$parsed_formula$status_col
  common_tune_args <- list(...)

  if ("ncores" %in% names(common_tune_args)) {
    warning("Ignoring `ncores` in `...`; use `inner_ncores` in `benchmark_tuned_survlearners()`.")
    common_tune_args$ncores <- NULL
  }

  set.seed(seed)
  outer_resamples <- rsample::vfold_cv(data, v = outer_folds, strata = !!rlang::sym(status_col))

  outer_results <- list()
  selected_params <- list()
  final_models <- list()

  for (learner in learners) {
    fit_fun_name <- paste0("fit_", learner)
    pred_fun_name <- paste0("predict_", learner)
    tune_fun_name <- paste0("tune_", learner)

    if (!exists(fit_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", fit_fun_name))
      next
    }
    if (!exists(pred_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", pred_fun_name))
      next
    }
    if (!exists(tune_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", tune_fun_name))
      next
    }

    fit_fun <- get(fit_fun_name)
    pred_fun <- get(pred_fun_name)
    tune_fun <- get(tune_fun_name)

    learner_spec <- learner_args[[learner]]
    if (is.null(learner_spec)) {
      learner_spec <- list()
    }
    if (!is.list(learner_spec)) {
      stop("Each entry in `learner_args` must be a list.")
    }

    has_named_tune_fit <- any(names(learner_spec) %in% c("tune", "fit"))
    tune_args <- if (has_named_tune_fit) learner_spec$tune else learner_spec
    fit_args_extra <- if (has_named_tune_fit) learner_spec$fit else list()

    if (is.null(tune_args)) {
      tune_args <- list()
    }
    if (is.null(fit_args_extra)) {
      fit_args_extra <- list()
    }

    if (!is.list(tune_args) || !is.list(fit_args_extra)) {
      stop("Learner-specific `tune` and `fit` entries must be lists.")
    }

    if ("ncores" %in% names(tune_args)) {
      warning(sprintf(
        "Ignoring learner-specific `ncores` for '%s'; use `inner_ncores` in `benchmark_tuned_survlearners()`.",
        learner
      ))
      tune_args$ncores <- NULL
    }

    tune_args <- c(common_tune_args, tune_args)

    if (verbose) {
      message("Nested CV learner: ", learner)
    }

    learner_results <- vector("list", length(outer_resamples$splits))
    learner_params <- vector("list", length(outer_resamples$splits))

    for (i in seq_along(outer_resamples$splits)) {
      split <- outer_resamples$splits[[i]]
      fold_id <- outer_resamples$id[[i]]
      train <- rsample::analysis(split)
      test <- rsample::assessment(split)

      if (verbose) {
        message("  Outer fold ", i, "/", length(outer_resamples$splits))
      }

      fold_result <- tryCatch({
        tuning_results <- .nested_surv_tune_once(
          tune_fun = tune_fun,
          formula = formula,
          data = train,
          times = times,
          metrics = metrics,
          inner_folds = inner_folds,
          seed = seed + i - 1L,
          inner_ncores = inner_ncores,
          tune_args = tune_args
        )

        best_param_df <- .select_cols(tuning_results[1, ], .nested_surv_param_cols(tuning_results))
        fit_args <- c(
          list(formula = formula, data = train),
          .nested_surv_row_to_list(best_param_df),
          fit_args_extra
        )

        model <- do.call(fit_fun, fit_args)
        pred <- pred_fun(model, newdata = test, times = times)
        pred <- .finalize_survmat(pred, times = times)

        scores <- .nested_surv_score_predictions(
          formula = formula,
          data = test,
          pred = pred,
          times = times,
          metrics = metrics
        )
        scores[, `:=`(learner = learner, id = fold_id, outer_fold = i)]
        data.table::setcolorder(scores, c("learner", "id", "outer_fold"))

        params <- data.table::as.data.table(best_param_df)
        params[, `:=`(learner = learner, id = fold_id, outer_fold = i)]
        data.table::setcolorder(params, c("learner", "id", "outer_fold"))

        list(scores = scores, params = params)
      }, error = function(e) {
        if (!suppress_errors) stop(e)
        warning(sprintf(
          "Learner '%s' failed on outer fold %s: %s",
          learner, i, conditionMessage(e)
        ))
        NULL
      })

      if (!is.null(fold_result)) {
        learner_results[[i]] <- fold_result$scores
        learner_params[[i]] <- fold_result$params
      }
    }

    outer_results[[learner]] <- data.table::rbindlist(learner_results, fill = TRUE)
    selected_params[[learner]] <- data.table::rbindlist(learner_params, fill = TRUE)

    if (isTRUE(refit_final)) {
      final_model <- tryCatch({
        final_tuning <- .nested_surv_tune_once(
          tune_fun = tune_fun,
          formula = formula,
          data = data,
          times = times,
          metrics = metrics,
          inner_folds = inner_folds,
          seed = seed + outer_folds,
          inner_ncores = inner_ncores,
          tune_args = tune_args
        )

        final_param_df <- .select_cols(final_tuning[1, ], .nested_surv_param_cols(final_tuning))
        fit_args <- c(
          list(formula = formula, data = data),
          .nested_surv_row_to_list(final_param_df),
          fit_args_extra
        )

        model <- do.call(fit_fun, fit_args)
        attr(model, "tuning_results") <- final_tuning
        model
      }, error = function(e) {
        if (!suppress_errors) stop(e)
        warning(sprintf(
          "Final refit for learner '%s' failed: %s",
          learner, conditionMessage(e)
        ))
        NULL
      })

      final_models[[learner]] <- final_model
    }
  }

  outer_results <- data.table::rbindlist(outer_results, fill = TRUE)
  selected_params <- data.table::rbindlist(selected_params, fill = TRUE)

  if (nrow(outer_results) == 0L) {
    stop("All learners failed or returned empty nested CV results.")
  }

  structure(
    list(
      outer_results = outer_results,
      outer_summary = summarise_benchmark(outer_results),
      selected_params = selected_params,
      final_models = if (isTRUE(refit_final)) final_models else NULL,
      settings = list(
        learners = learners,
        metrics = metrics,
        times = times,
        outer_folds = outer_folds,
        inner_folds = inner_folds,
        inner_ncores = inner_ncores,
        seed = seed
      )
    ),
    class = "nested_surv_benchmark"
  )
}



#' Summarise Benchmark Results (Mean  SD with Wald CI)
#'
#' Aggregates cross‑validated benchmark results by learner and metric, reporting
#' mean, standard deviation, standard error, and an approximate 95% Wald
#' confidence interval.
#'
#' @param benchmark_results A data frame produced by
#'   \code{benchmark_default_survlearners()}, containing at least
#'   \code{learner}, \code{metric}, and \code{value}.
#' @param digits Integer number of decimal places for \code{mean}, \code{sd},
#'   \code{se}, \code{lower}, \code{upper} (default \code{3}).
#'
#' @return A data.table with columns \code{learner}, \code{metric}, \code{mean},
#' \code{sd}, \code{n}, \code{se}, \code{lower}, \code{upper}.
#'
#' @examples
#' res <- data.frame(
#'   learner = c("coxph", "coxph", "rpart", "rpart"),
#'   metric = c("cindex", "ibs", "cindex", "ibs"),
#'   value = c(0.64, 0.19, 0.60, 0.23)
#' )
#' summarise_benchmark(res)
#'
#' @seealso [benchmark_default_survlearners()], [plot_benchmark()],
#'   [summarize_benchmark_results()]
#' @export

summarise_benchmark <- function(benchmark_results, digits = 3) {
  DT <- data.table::as.data.table(benchmark_results)
  out <- DT[, list(
    mean = mean(value, na.rm = TRUE),
    sd   = stats::sd(value, na.rm = TRUE),
    n    = .N
  ), by = list(learner, metric)]
  out[, se := sd / sqrt(n)]
  out[, lower := mean - 1.96 * se]
  out[, upper := mean + 1.96 * se]
  out[, c("mean", "sd", "se", "lower", "upper") :=
    lapply(.SD, round, digits = digits), .SDcols = c("mean", "sd", "se", "lower", "upper")]
  out[]
}



#' Plot Benchmark Distributions Across Learners
#'
#' Produces box‑and‑jitter plots of CV metric values per learner, faceted by
#' metric for quick visual comparison.
#'
#' @param benchmark_results A data frame from
#'   \code{benchmark_default_survlearners()} with columns \code{learner},
#'   \code{metric}, and \code{value}.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely (e.g., for journals
#'   requiring caption-only figures).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' res <- data.frame(
#'   learner = c("coxph", "coxph", "rpart", "rpart"),
#'   metric = c("cindex", "ibs", "cindex", "ibs"),
#'   value = c(0.64, 0.19, 0.60, 0.23)
#' )
#' plot_benchmark(res)
#'
#' @seealso [benchmark_default_survlearners()], [summarise_benchmark()]
#' @export

plot_benchmark <- function(benchmark_results, title) {
  if (missing(title)) title <- "Cross-Validated Performance of Survival Learners"
  p <- ggplot(benchmark_results, ggplot2::aes(x = learner, y = value)) +
    geom_boxplot(fill = .survalis_palette[1], alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      x = "Learner",
      y = "Metric Value"
    ) +
    theme_survalis()
  if (!is.null(title)) p <- p + labs(title = title)
  p
}



#' Compact Table of Mean  SD by Learner and Metric
#'
#' Creates a wide table summarizing each learner's performance as
#' \code{mean  sd} per metric, suitable for reporting.
#'
#' @param results A data frame with columns \code{learner}, \code{metric},
#'   and \code{value} (e.g., output of \code{benchmark_default_survlearners()}).
#' @param digits Integer number of decimal places in the formatted summary
#'   (default \code{3}).
#'
#' @return A wide tibble with one row per learner and one column per metric,
#'   containing formatted strings \code{"mean  sd"}.
#'
#' @examples
#' res <- data.frame(
#'   learner = c("coxph", "coxph", "rpart", "rpart"),
#'   metric = c("cindex", "ibs", "cindex", "ibs"),
#'   value = c(0.64, 0.19, 0.60, 0.23)
#' )
#' summarize_benchmark_results(res, digits = 2)
#'
#' @seealso [summarise_benchmark()], [benchmark_default_survlearners()]
#' @export

summarize_benchmark_results <- function(results, digits = 3) {

  stopifnot(is.data.frame(results), all(c("learner", "metric", "value") %in% colnames(results)))

  DT <- data.table::as.data.table(results)
  summarised <- DT[, list(
    mean = mean(value, na.rm = TRUE),
    sd   = stats::sd(value, na.rm = TRUE)
  ), by = list(learner, metric)]
  summarised[, summary := sprintf("%.*f  %.*f", digits, mean, digits, sd)]

  out <- data.table::dcast(summarised, learner ~ metric, value.var = "summary")
  data.table::setorderv(out, "learner")
  out[]
}





#' Select the Best Survival Learner by a Given Metric
#'
#' Extracts the top‑performing learner(s) under a chosen metric from benchmark
#' results, using the average value across folds.
#'
#' @param benchmark_results A data frame with columns \code{learner},
#'   \code{metric}, and \code{value}.
#' @param metric Character name of the metric to optimize (e.g., \code{"cindex"},
#'   \code{"ibs"}). Must exist in \code{benchmark_results$metric}.
#' @param maximize Logical; whether to maximize the metric. If \code{NULL}
#'   (default), the function maximizes for concordance‑like metrics
#'   (\code{"cindex"}) and minimizes for error‑like metrics
#'   (\code{"ibs"}, \code{"brier"}, \code{"iae"}, \code{"ise"}, \code{"ece"}).
#'
#' @return A tibble with columns \code{learner}, \code{metric}, and the selected
#'   average \code{value} for the best learner(s). Ties are returned as multiple rows.
#'
#' @examples
#' res <- data.frame(
#'   learner = c("coxph", "coxph", "rpart", "rpart"),
#'   metric = c("cindex", "ibs", "cindex", "ibs"),
#'   value = c(0.64, 0.19, 0.60, 0.23)
#' )
#' best_survlearner(res, metric = "cindex")
#' best_survlearner(res, metric = "ibs")
#'
#' @seealso [benchmark_default_survlearners()], [summarise_benchmark()]
#' @export

best_survlearner <- function(benchmark_results, metric, maximize = NULL) {
  if (!metric %in% benchmark_results$metric) {
    stop("Metric not found in results: ", metric)
  }

  if (is.null(maximize)) {
    # default: maximize cindex, minimize others
    maximize <- !(metric %in% c("ibs", "brier", "iae", "ise", "ece"))
  }

  target_metric <- metric
  DT <- data.table::as.data.table(benchmark_results)
  summary <- DT[metric == target_metric, list(value = mean(value, na.rm = TRUE)), by = list(learner, metric)]

  best <- if (maximize) summary[value == max(value)] else summary[value == min(value)]

  return(best[])
  }
