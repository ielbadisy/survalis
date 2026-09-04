#' Tune a survival learner over a hyperparameter grid
#'
#' Grid search with cross-validation for learners that expose a `tune_<model>()`
#' routine (`list_survlearners(has_tune = TRUE)`). The best configuration is
#' refit on the full data and returned as a ready-to-use `survalis_fit`.
#'
#' @param formula A survival formula.
#' @param data A data frame.
#' @param model A tunable learner id.
#' @param grid A data frame / list of hyperparameter combinations. `NULL` uses
#'   the learner's built-in default grid.
#' @param times Evaluation grid; `NULL` uses [default_times()].
#' @param resampling A [cv()] spec (only `cv()` is supported here); its `v` and
#'   `seed` drive the inner cross-validation.
#' @param metric A single metric to optimise (default `"ibs"`).
#' @param ncores Cores for the inner CV (default `1`).
#' @param ... Forwarded to `tune_<model>()` / `fit_<model>()`.
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
  tunable <- list_survlearners(has_tune = TRUE)$learner
  if (length(model) != 1L || !model %in% tunable) {
    stop("Model ", sQuote(model %||% "NULL"), " is not tunable. See ",
         "list_survlearners(has_tune = TRUE).", call. = FALSE)
  }
  if (!inherits(resampling, "survalis_resampling") || resampling$type != "cv") {
    stop("tune() currently supports resampling = cv() only.", call. = FALSE)
  }
  if (length(metric) != 1L) stop("`metric` must be a single value.", call. = FALSE)

  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula))
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)

  tuner <- get(paste0("tune_", model), envir = asNamespace("survalis"))
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
    if ("failed" %in% names(results)) {
      keep <- !as.logical(results$failed)
      if (any(keep)) results <- results[keep, , drop = FALSE]
    }
    results <- as.data.frame(.arrange_by_metric_dt(results, metric, .metric_maximize(metric)))
    param_cols0 <- setdiff(names(results), c(intersect(names(results), list_metrics()$metric), "failed"))
    spec0 <- as.list(results[1, param_cols0, drop = FALSE])
    fit_fun <- get(paste0("fit_", model), envir = asNamespace("survalis"))
    best_model <- do.call(fit_fun, c(list(formula = formula, data = data), spec0))
  }

  metric_cols <- intersect(names(results), list_metrics()$metric)
  param_cols <- setdiff(names(results), c(metric_cols, "failed"))
  best_row <- results[1, , drop = FALSE]
  spec <- as.list(best_row[, param_cols, drop = FALSE])

  fit_obj <- .as_survalis_fit(best_model, model_id = model, spec = spec,
                              formula = formula, data = data, parsed = parsed)

  structure(
    list(fit = fit_obj, results = as.data.frame(results), best = best_row,
         times = times, resampling = resampling, metric = metric, model = model),
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
evaluate.survalis_tune <- function(object, times = NULL, resampling = cv(),
                                   metrics = c("cindex", "ibs"), ncores = 1, ...) {
  evaluate(object$fit, times = times, resampling = resampling,
           metrics = metrics, ncores = ncores, ...)
}
