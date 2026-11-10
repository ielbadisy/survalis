#' Tune a survival learner over a hyperparameter grid
#'
#' Grid search over resampling folds. Learners listed by
#' `list_survlearners(has_tune = TRUE)` have a built-in default grid; any other
#' learner, including a custom one (see [fit()]), is tuned over the `grid` you
#' supply, with every configuration scored on the same splits. The best
#' configuration is refit on the full data and returned as a ready-to-use
#' `survalis_fit`.
#'
#' @param formula A survival formula.
#' @param data A data frame.
#' @param model A learner id (`list_survlearners()$learner`) or a custom
#'   learner.
#' @param grid A data frame of hyperparameter combinations (one row per
#'   configuration, for example from [expand.grid()]). `NULL` uses the learner's
#'   built-in default grid, if it has one. A `grid` need not name every tuned
#'   parameter: parameters it leaves out take the learner's default values
#'   (those of [fit()] with no `spec`).
#' @param times Evaluation grid; `NULL` picks one from the data quantiles.
#' @param resampling A resampling spec. The built-in tuner uses [cv()] only
#'   (its `v` and `seed` drive the inner cross-validation). It runs when `grid`
#'   is `NULL` or has exactly the columns of the built-in grid; any other
#'   `grid` is scored on shared splits and accepts any spec.
#' @param metric A single metric to optimise (default `"ibs"`).
#' @param ncores Cores for the inner CV (default `1`).
#' @param ... Fixed learner arguments: forwarded to the built-in tuner, or
#'   added to every configuration of a supplied `grid`.
#'
#' @return A `survalis_tune` object: `fit` (the refit `survalis_fit` at the best
#'   configuration), `results` (the scored grid), `best` (one-row best config),
#'   `times`, `resampling`, `metric`, `model`.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' tuned <- tune(Surv(time, status) ~ age + karno, veteran, "ranger",
#'               grid = expand.grid(num.trees = c(100, 300), mtry = 2, min.node.size = c(5, 15)),
#'               times = c(90, 180), resampling = cv(3), metric = "ibs")
#' tuned
#' predict(tuned, veteran[1:3, ], times = c(90, 180))
#' @export
tune <- function(formula, data, model, grid = NULL, times = NULL,
                 resampling = cv(), metric = "ibs", ncores = 1, ...) {
  lrn <- .resolve_learner(model)
  if (length(metric) != 1L) stop("`metric` must be a single value.", call. = FALSE)
  if (!.is_metric(metric)) {
    stop("Unknown metric ", sQuote(metric), ". See list_metrics().", call. = FALSE)
  }
  if (!inherits(resampling, "survalis_resampling")) {
    stop("`resampling` must come from cv(), holdout(), group_cv() or bootstrap().",
         call. = FALSE)
  }
  builtin <- !is.null(lrn$tune)
  if (!builtin && is.null(grid)) {
    stop("Learner ", sQuote(lrn$id), " has no built-in grid; supply `grid`.",
         call. = FALSE)
  }
  # The built-in tuners pass every tuned parameter by name, so they only take a
  # grid with exactly their own columns. Any other grid (one that leaves out a
  # tuned parameter, or adds an engine argument) goes through the generic
  # search, where parameters not in the grid keep the learner's defaults.
  if (builtin && !is.null(grid) &&
      !setequal(names(grid), .builtin_grid_names(lrn$tune))) {
    builtin <- FALSE
  }
  if (builtin && resampling$type != "cv") {
    stop("The built-in tuner of ", sQuote(lrn$id), " supports resampling = cv() only.",
         call. = FALSE)
  }


  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula))
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)
  .check_metrics(metric, times)

  tuned_names <- if (!is.null(grid)) names(grid) else if (builtin) .builtin_grid_names(lrn$tune)
  clash <- intersect(names(list(...)), tuned_names)
  if (length(clash)) {
    stop("Argument(s) ", paste(sQuote(clash), collapse = ", "), " are both tuned ",
         "in the grid and fixed in `...`; drop one of them.", call. = FALSE)
  }

  # The built-in tuners run one plain v-fold split. For repeated or
  # differently stratified CV, search their default grid on shared splits.
  if (builtin && (.null_default(resampling$params$repeats, 1L) > 1L ||
                  !identical(resampling$params$strata, "status"))) {
    if (is.null(grid)) {
      grid <- eval(formals(lrn$tune)$param_grid, envir = environment(lrn$tune))
    }
    builtin <- FALSE
  }

  if (!builtin) {
    return(.tune_grid(formula, data, model, lrn, grid, times, resampling, metric,
                      parsed, ncores, fixed = list(...)))
  }

  tuner <- lrn$tune
  tf <- names(formals(tuner))
  std <- .std_outcome(formula, data)
  args <- list(formula = std$formula, data = std$data, times = times,
               metrics = metric, folds = resampling$params$v,
               seed = .null_default(resampling$params$seed, 123L))
  if ("ncores" %in% tf) args$ncores <- ncores
  if (!is.null(grid)) args$param_grid <- grid
  dots <- list(...)
  if (!"..." %in% tf) dots <- list()
  if (length(dots)) args <- c(args, dots)

  # Ask the tuner for its scored grid only and refit here: several tuners
  # return just the refit model (without the grid) when refit_best = TRUE.
  if ("refit_best" %in% tf) args$refit_best <- FALSE
  results <- do.call(tuner, args)
  results <- if (is.null(results)) data.frame() else as.data.frame(results)
  class(results) <- "data.frame"
  if (!nrow(results)) {
    stop("Tuning ", sQuote(lrn$id), " returned no results: every configuration ",
         "failed during cross-validation. Fit one configuration with fit() to ",
         "see the underlying error.", call. = FALSE)
  }
  if ("failed" %in% names(results)) {
    failed <- as.logical(results$failed)
    if (all(failed)) {
      stop("All ", nrow(results), " tuning configurations of ", sQuote(lrn$id),
           " failed during cross-validation. Fit one configuration with fit() ",
           "to see the underlying error.", call. = FALSE)
    }
    # the tuner has already warned about each failed configuration
    results <- results[!failed, setdiff(names(results), "failed"), drop = FALSE]
  }
  if (!metric %in% names(results)) {
    stop("The tuner of ", sQuote(lrn$id), " did not report ", sQuote(metric), ".",
         call. = FALSE)
  }
  results <- as.data.frame(.arrange_by_metric_dt(results, metric, .metric_maximize(metric)))
  class(results) <- "data.frame"
  rownames(results) <- NULL
  param_cols0 <- setdiff(names(results), names(results)[.is_metric(names(results))])
  spec0 <- as.list(results[1, param_cols0, drop = FALSE])
  best_model <- do.call(lrn$fit, c(list(formula = formula, data = data), spec0, dots))

  metric_cols <- names(results)[.is_metric(names(results))]
  param_cols <- setdiff(names(results), c(metric_cols, "failed"))
  best_row <- results[1, , drop = FALSE]
  # keep the fixed engine arguments, so fit$spec reproduces the best model
  spec <- c(as.list(best_row[, param_cols, drop = FALSE]), dots)

  fit_obj <- .as_survalis_fit(best_model, model_id = lrn$id, spec = spec,
                              formula = formula, data = data, parsed = parsed)

  structure(
    list(fit = fit_obj, results = as.data.frame(results), best = best_row,
         times = times, resampling = resampling, metric = metric, model = model),
    class = "survalis_tune"
  )
}

#' Generic grid search on shared splits (learners without a built-in tuner)
#' @noRd
.tune_grid <- function(formula, data, model, lrn, grid, times, resampling,
                       metric, parsed, ncores = 1, fixed = list()) {
  grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  for (nm in names(grid)) if (is.factor(grid[[nm]])) grid[[nm]] <- as.character(grid[[nm]])
  if (!nrow(grid)) stop("`grid` has no rows.", call. = FALSE)
  splits <- .make_splits(resampling, data, parsed$status_col)
  vals <- vapply(seq_len(nrow(grid)), function(i) {
    spec <- c(as.list(grid[i, , drop = FALSE]), fixed)
    per <- tryCatch(.evaluate_splits(splits, formula, data, model, spec, times,
                                     metric, parsed, ncores),
                    error = function(e) .tune_config_failed(lrn$id, grid[i, , drop = FALSE], e))
    if (is.null(per)) NA_real_ else mean(per$value, na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(vals))) {
    stop("All ", nrow(grid), " tuning configurations failed. Fit one ",
         "configuration with fit() to see the underlying error.", call. = FALSE)
  }
  results <- grid
  results[[metric]] <- vals
  results <- results[!is.na(vals), , drop = FALSE]
  ord <- order(results[[metric]], decreasing = .metric_maximize(metric))
  results <- results[ord, , drop = FALSE]
  rownames(results) <- NULL
  spec <- c(as.list(results[1, names(grid), drop = FALSE]), fixed)
  best_model <- do.call(lrn$fit, c(list(formula = formula, data = data), spec))
  fit_obj <- .as_survalis_fit(best_model, model_id = lrn$id, spec = spec,
                              formula = formula, data = data, parsed = parsed)
  if (lrn$custom) fit_obj$learner_def <- model
  structure(
    list(fit = fit_obj, results = results, best = results[1, , drop = FALSE],
         times = times, resampling = resampling, metric = metric, model = lrn$id),
    class = "survalis_tune"
  )
}

#' Column names of a built-in tuner's default grid
#' @noRd
.builtin_grid_names <- function(tuner) {
  g <- eval(formals(tuner)$param_grid, envir = environment(tuner))
  names(g)
}

#' @export
print.survalis_tune <- function(x, ...) {
  cat(sprintf("<survalis_tune> model = %s  |  optimise %s over %d configs\n",
              x$model, x$metric, nrow(x$results)))
  print(x$resampling)
  cat("  best: ",
      paste(names(x$fit$spec),
            vapply(x$fit$spec, function(v) paste(deparse(v), collapse = ""), ""),
            sep = " = ", collapse = ", "), "\n", sep = "")
  if (x$metric %in% names(x$best)) {
    cat(sprintf("  %s at best: %.4f\n", x$metric, x$best[[x$metric]]))
  }
  invisible(x)
}

#' @export
summary.survalis_tune <- function(object, ...) object$results

#' @export
plot.survalis_tune <- function(x, title, ...) {
  res <- x$results
  if (!x$metric %in% names(res)) stop("Metric column not found in tuning results.")
  param_cols <- setdiff(names(res), names(res)[.is_metric(names(res))])
  varying <- param_cols[vapply(param_cols, function(c) length(unique(res[[c]])) > 1L, logical(1))]
  xcol <- if (length(varying)) varying[1] else param_cols[1]
  if (missing(title)) title <- sprintf("Tuning %s: %s vs %s", x$model, x$metric, xcol)
  p <- ggplot2::ggplot(res, ggplot2::aes(x = .data[[xcol]], y = .data[[x$metric]])) +
    ggplot2::geom_point(color = .survalis_palette[1]) +
    ggplot2::geom_line(color = .survalis_palette[1], alpha = 0.5) +
    ggplot2::labs(x = xcol, y = x$metric) +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}

#' @method predict survalis_tune
#' @export
predict.survalis_tune <- function(object, newdata = NULL, times = NULL,
                                  type = "survival", ...) {
  predict(object$fit, newdata = newdata, times = times, type = type, ...)
}

#' @export
evaluate.survalis_tune <- function(object, newdata = NULL, times = NULL,
                                   resampling = cv(),
                                   metrics = c("cindex", "ibs"), ncores = 1, ...) {
  evaluate(object$fit, newdata = newdata, times = times, resampling = resampling,
           metrics = metrics, ncores = ncores, ...)
}
