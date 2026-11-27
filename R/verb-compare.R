#' Compare several survival learners on shared resampling folds
#'
#' Evaluates each learner in `models` on the *same* set of resampling splits, so
#' the per-fold values are paired across models. Returns a leaderboard plus the
#' full per-split table.
#'
#' @param formula A survival formula.
#' @param data A data frame.
#' @param models Character vector of learner ids, or a list mixing ids and
#'   custom learners (see [fit()]). List names, when given, label the learners.
#' @param specs Named list of engine-argument lists, keyed by learner label
#'   (`specs[["ranger"]] <- list(num.trees = 500)`).
#' @param times Evaluation grid; `NULL` picks up to 50 event-time quantiles (5th to 95th percentile).
#' @param resampling A [resampling] spec (default `cv()`); the split set is drawn
#'   once and reused for every model.
#' @param metrics Character vector from `listmetrics()$metric`.
#' @param tune Logical. When `TRUE`, each tunable model is tuned once on the full
#'   data (its default grid) and the winning configuration is then evaluated on
#'   the shared folds. This is not nested CV, so the tuned results are optimistic;
#'   a message is emitted. Default `FALSE`.
#' @param ncores Number of parallel workers for the split loop (default `1`). Workers are
#'   fresh R sessions (PSOCK), each with one OpenMP thread, so parallel runs are safe after an
#'   OpenMP engine (xgboost, rsf) ran in the session. A custom learner must therefore not rely
#'   on objects of the global environment other than its own arguments.
#' @param ... Unused.
#'
#' @return A `survalis_compare` object: `leaderboard` (model x metric, mean and
#'   se), `per_split` (long), `times`, `resampling`, `metrics`, `models`, `tuned`.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' cmp <- compare(Surv(time, status) ~ age + karno + celltype, veteran,
#'                models = c("coxph", "ranger"), times = c(90, 180, 365),
#'                resampling = cv(3))
#' cmp
#' plot(cmp)
#' @export
compare <- function(formula, data, models, specs = list(), times = NULL,
                    resampling = cv(), metrics = c("cindex", "ibs"),
                    tune = FALSE, ncores = 1, ...) {
  if (.is_custom_learner(models)) models <- list(models)
  models <- as.list(models)
  if (length(models) < 2L) stop("`models` must name at least two learners.", call. = FALSE)
  for (m in models) .resolve_learner(m)
  labels <- vapply(models, .learner_label, character(1))
  if (!is.null(names(models))) {
    nm <- names(models)
    labels[nzchar(nm)] <- nm[nzchar(nm)]
  }
  if (anyDuplicated(labels)) {
    stop("Learner labels must be unique; got ", paste(labels, collapse = ", "),
         ". Name the list elements or give custom learners a `name`.", call. = FALSE)
  }


  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula))
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)
  .check_metrics(metrics, times, multi = TRUE)
  if (.start_workers(ncores)) on.exit(.stop_workers(), add = TRUE)

  splits <- .make_splits(resampling, data, parsed$status_col)
  tunable <- listlearners(has_tune = TRUE)$learner

  per <- list()
  for (k in seq_along(models)) {
    mdl <- models[[k]]
    lab <- labels[k]
    spec <- as.list(.null_default(specs[[lab]], list()))
    if (isTRUE(tune) && is.character(mdl) && mdl %in% tunable) {
      tuned <- tune(formula, data, mdl, times = times,
                    resampling = if (resampling$type == "cv") resampling else cv(seed = 1),
                    metric = metrics[1], ncores = ncores)
      spec <- tuned$fit$spec
      message("compare(): '", mdl, "' tuned on the full data (non-nested); ",
              "its resampled scores are optimistic.")
    }
    d <- .evaluate_splits(splits, formula, data, mdl, spec, times, metrics, parsed, ncores)
    d$model <- lab
    per[[lab]] <- d
  }
  per <- basetable::rbindfill(unname(per), fill = TRUE)

  mean_tbl <- basetable::aggregate(per, by = c("model", "metric"), value = "value", fun = mean)
  se_tbl <- basetable::aggregate(per, by = c("model", "metric"), value = "value",
                                 fun = function(v) stats::sd(v) / sqrt(length(v)))
  names(mean_tbl)[names(mean_tbl) == "value"] <- "mean"
  names(se_tbl)[names(se_tbl) == "value"] <- "se"
  leaderboard <- basetable::merge(mean_tbl, se_tbl, by = c("model", "metric"))
  leaderboard <- basetable::orderrows(leaderboard, by = c("metric", "model"))
  if (any(grepl("@", leaderboard$metric))) leaderboard <- .order_bytime(leaderboard)
  rownames(leaderboard) <- NULL

  structure(
    list(leaderboard = leaderboard, per_split = per, times = times,
         resampling = resampling, metrics = metrics, models = labels,
         tuned = isTRUE(tune), formula = formula),
    class = "survalis_compare"
  )
}

#' @keywords internal
.metric_maximize <- function(metric) {
  if (grepl("^range:", metric)) return(FALSE)
  m <- listmetrics()
  isTRUE(m$direction[match(sub("@.*", "", metric), m$metric)] == "maximize")
}

#' @export
print.survalis_compare <- function(x, ...) {
  cat(sprintf("<survalis_compare> %d models%s\n", length(x$models),
              if (x$tuned) " (tuned, non-nested)" else ""))
  print(x$resampling)
  for (mt in x$metrics) {
    sub <- x$leaderboard[x$leaderboard$metric == mt, , drop = FALSE]
    ord <- order(sub$mean, decreasing = .metric_maximize(mt))
    sub <- sub[ord, , drop = FALSE]
    cat(sprintf("\n  %s (%s is better):\n", mt,
                if (.metric_maximize(mt)) "higher" else "lower"))
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("   %s%-16s %.4f  (se %.4f)\n",
                  if (i == 1L) "* " else "  ", sub$model[i], sub$mean[i], sub$se[i]))
    }
  }
  invisible(x)
}

#' @export
summary.survalis_compare <- function(object, ...) {
  basetable::towide(object$leaderboard[, c("model", "metric", "mean")],
                    names = "metric", values = "mean", idcols = "model")
}

#' Plot a model comparison
#'
#' One panel per metric, with the per-split values of each model as a box plot
#' and points. Each panel has its own value axis, and its strip says whether
#' higher or lower is better. Models are ordered by the first metric, best on top.
#'
#' @param x A `survalis_compare` object.
#' @param title,subtitle,caption Plot title, subtitle and caption; `NULL`
#'   (default) draws none.
#' @param scales Passed to [ggplot2::facet_wrap()]: `"free_x"` (default) gives
#'   each metric its own value axis, `"fixed"` one shared axis.
#' @param ... Unused.
#' @return A ggplot.
#' @export
plot.survalis_compare <- function(x, title = NULL, scales = "free_x", subtitle = NULL,
                                  caption = NULL, ...) {
  d <- x$per_split
  mts <- unique(d$metric)
  d$metric <- factor(d$metric, mts, .metric_strip(mts))
  lb <- x$leaderboard[x$leaderboard$metric == mts[1], ]
  # ggplot draws the first factor level at the bottom, so the best model goes last
  d$model <- factor(d$model, lb$model[order(lb$mean, decreasing = !.metric_maximize(mts[1]))])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = model)) +
    ggplot2::geom_boxplot(fill = "grey90", colour = "grey40", outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0, height = 0.12, alpha = 0.6, size = 1.4, colour = "grey20") +
    ggplot2::facet_wrap(~ metric, scales = scales) +
    ggplot2::labs(x = "Value over resampling splits", y = NULL, title = title,
                  subtitle = subtitle, caption = caption) +
    theme_survalis()
  if (scales %in% c("free", "free_x")) p <- p + .value_axis()
  p
}
