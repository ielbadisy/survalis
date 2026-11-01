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
#'   built-in default grid, if it has one.
#' @param times Evaluation grid; `NULL` picks one from the data quantiles.
#' @param resampling A resampling spec. Built-in grids use [cv()] only (its `v`
#'   and `seed` drive the inner cross-validation); a supplied `grid` for a
#'   learner without a built-in tuner accepts any spec.
#' @param metric A single metric to optimise (default `"ibs"`).
#' @param ncores Cores for the inner CV (default `1`).
#' @param ... Forwarded to the learner's built-in tuner.
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
  if (!metric %in% list_metrics()$metric) {
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
  if (builtin && resampling$type != "cv") {
    stop("The built-in tuner of ", sQuote(lrn$id), " supports resampling = cv() only.",
         call. = FALSE)
  }


  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula))
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)

  if (!builtin) {
    return(.tune_grid(formula, data, model, lrn, grid, times, resampling, metric,
                      parsed, ncores))
  }

  tuner <- lrn$tune
  tf <- names(formals(tuner))
  args <- list(formula = formula, data = data, times = times,
               metrics = metric, folds = resampling$params$v,
               seed = .null_default(resampling$params$seed, 123L))
  if ("ncores" %in% tf) args$ncores <- ncores
  if (!is.null(grid)) args$param_grid <- grid
  dots <- list(...)
  if (length(dots) && "..." %in% tf) args <- c(args, dots)

  if ("refit_best" %in% tf) {
    args$refit_best <- TRUE
    best_model <- do.call(tuner, args)
    results <- as.data.frame(attr(best_model, "tuning_results"))
  } else {
    results <- as.data.frame(do.call(tuner, args))
    if ("failed" %in% names(results) && all(as.logical(results$failed))) {
      stop("All ", nrow(results), " tuning configurations of ", sQuote(model),
           " failed during cross-validation. Fit one configuration with fit() ",
           "to see the underlying error.", call. = FALSE)
    }
    if ("failed" %in% names(results)) {
      keep <- !as.logical(results$failed)
      if (any(keep)) results <- results[keep, , drop = FALSE]
    }
    results <- as.data.frame(.arrange_by_metric_dt(results, metric, .metric_maximize(metric)))
    param_cols0 <- setdiff(names(results), c(intersect(names(results), list_metrics()$metric), "failed"))
    spec0 <- as.list(results[1, param_cols0, drop = FALSE])
    best_model <- do.call(lrn$fit, c(list(formula = formula, data = data), spec0))
  }

  metric_cols <- intersect(names(results), list_metrics()$metric)
  param_cols <- setdiff(names(results), c(metric_cols, "failed"))
  best_row <- results[1, , drop = FALSE]
  spec <- as.list(best_row[, param_cols, drop = FALSE])

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
                       metric, parsed, ncores = 1) {
  grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  for (nm in names(grid)) if (is.factor(grid[[nm]])) grid[[nm]] <- as.character(grid[[nm]])
  if (!nrow(grid)) stop("`grid` has no rows.", call. = FALSE)
  splits <- .make_splits(resampling, data, parsed$status_col)
  vals <- vapply(seq_len(nrow(grid)), function(i) {
    spec <- as.list(grid[i, , drop = FALSE])
    per <- tryCatch(.evaluate_splits(splits, formula, data, model, spec, times,
                                     metric, parsed, ncores),
                    error = function(e) NULL)
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
  spec <- as.list(results[1, names(grid), drop = FALSE])
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
  param_cols <- setdiff(names(res), intersect(names(res), list_metrics()$metric))
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
