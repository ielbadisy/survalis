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
#' @param metrics Character vector from `list_metrics()$metric`.
#' @param tune Logical. When `TRUE`, each tunable model is tuned once on the full
#'   data (its default grid) and the winning configuration is then evaluated on
#'   the shared folds. This is not nested CV, so the tuned results are optimistic;
#'   a message is emitted. Default `FALSE`.
#' @param ncores Cores for the split loop (default `1`).
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

  splits <- .make_splits(resampling, data, parsed$status_col)
  tunable <- list_survlearners(has_tune = TRUE)$learner

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
  per <- do.call(rbind, per)

  agg <- stats::aggregate(value ~ model + metric, data = per,
                          FUN = function(v) c(mean = mean(v), se = stats::sd(v) / sqrt(length(v))))
  leaderboard <- data.frame(
    model = agg$model, metric = agg$metric,
    mean = round(agg$value[, "mean"], 4), se = round(agg$value[, "se"], 4)
  )

  structure(
    list(leaderboard = leaderboard, per_split = per, times = times,
         resampling = resampling, metrics = metrics, models = labels,
         tuned = isTRUE(tune), formula = formula),
    class = "survalis_compare"
  )
}

#' @keywords internal
.metric_maximize <- function(metric) {
  m <- list_metrics()
  isTRUE(m$direction[match(metric, m$metric)] == "maximize")
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
  stats::reshape(object$leaderboard[, c("model", "metric", "mean")],
                 idvar = "model", timevar = "metric", direction = "wide")
}

#' @export
plot.survalis_compare <- function(x, title, ...) {
  if (missing(title)) title <- "Model comparison"
  p <- ggplot2::ggplot(x$per_split,
                       ggplot2::aes(x = stats::reorder(model, value), y = value)) +
    ggplot2::geom_boxplot(fill = .survalis_palette[1], alpha = 0.3, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.5) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(x = "Model", y = "Value") +
    ggplot2::coord_flip() +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}
