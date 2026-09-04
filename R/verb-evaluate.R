#' Resampled performance of a survival learner
#'
#' Refits `model` on each training split of `resampling` and scores its
#' predictions on the held-out split, returning per-split metric values plus a
#' summary. A `survalis_fit` may be passed directly, in which case its stored
#' formula, data and `spec` are reused.
#'
#' @param object A `survalis_fit` from [fit()], or a survival `formula`.
#' @param data When `object` is a formula: the data frame.
#' @param model When `object` is a formula: a learner id
#'   (`list_survlearners()$learner`).
#' @param spec When `object` is a formula: engine arguments (named list).
#' @param times Evaluation grid; `NULL` uses [default_times()] and records it.
#' @param resampling A [resampling] spec (default `cv()`).
#' @param metrics Character vector from `list_metrics()$metric`. `"brier"` and
#'   `"ece"` require a single `times` value.
#' @param ncores Cores for the split loop (default `1`).
#' @param ... Unused.
#'
#' @return A `survalis_eval` object: `per_split` (long data frame), `summary`
#'   (mean / sd / se / 95% interval per metric), `times`, `resampling`,
#'   `metrics`, `model`.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
#'          spec = list(num.trees = 100))
#' evaluate(m, times = c(90, 180, 365), resampling = cv(3))
#' @export
evaluate <- function(object, ...) UseMethod("evaluate")

#' @rdname evaluate
#' @export
evaluate.survalis_fit <- function(object, times = NULL, resampling = cv(),
                                  metrics = c("cindex", "ibs"), ncores = 1, ...) {
  .evaluate_core(object$formula, object$data, object$model_id, object$spec,
                 times = times, resampling = resampling, metrics = metrics,
                 ncores = ncores)
}

#' @rdname evaluate
#' @export
evaluate.formula <- function(object, data, model, spec = list(), times = NULL,
                             resampling = cv(), metrics = c("cindex", "ibs"),
                             ncores = 1, ...) {
  .evaluate_core(object, data, model, as.list(spec),
                 times = times, resampling = resampling, metrics = metrics,
                 ncores = ncores)
}

#' @keywords internal
.recode_status_vec <- function(parsed, df) {
  s <- df[[parsed$status_col]]
  if (isTRUE(parsed$recode_status)) as.integer(s == parsed$event_value) else as.integer(s != 0)
}

#' Shared split-loop used by evaluate() and compare()
#' @keywords internal
.evaluate_splits <- function(splits, formula, data, model, spec, times, metrics,
                             parsed, ncores = 1) {
  fit_fun  <- get(paste0("fit_", model), envir = asNamespace("survalis"))
  pred_fun <- get(paste0("predict_", model), envir = asNamespace("survalis"))

  one <- function(sp) {
    tr <- data[sp$analysis, , drop = FALSE]
    te <- data[sp$assessment, , drop = FALSE]
    m  <- do.call(fit_fun, c(list(formula = formula, data = tr), spec))
    p  <- .finalize_survmat(pred_fun(m, newdata = te, times = times), times = times)
    surv <- survival::Surv(te[[parsed$time_col]], .recode_status_vec(parsed, te))
    sc <- .score_metrics(surv, p, times, metrics)
    data.frame(id = sp$id, metric = sc$metric, value = sc$value,
               stringsAsFactors = FALSE)
  }

  parts <- if (ncores > 1) {
    functionals::fmapn(list(sp = splits), function(sp) one(sp), ncores = ncores, pb = FALSE)
  } else {
    lapply(splits, one)
  }
  do.call(rbind, parts)
}

#' @keywords internal
.evaluate_core <- function(formula, data, model, spec, times, resampling,
                           metrics, ncores = 1) {
  op <- options(survalis.quiet_legacy = TRUE)
  on.exit(options(op), add = TRUE)
  known <- list_survlearners()$learner
  if (!model %in% known) stop("Unknown model ", sQuote(model), ".", call. = FALSE)
  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula))
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)

  splits <- .make_splits(resampling, data, parsed$status_col)
  per <- .evaluate_splits(splits, formula, data, model, as.list(spec), times,
                          metrics, parsed, ncores)
  summ <- cv_summary(per)

  structure(
    list(per_split = per, summary = as.data.frame(summ), times = times,
         resampling = resampling, metrics = metrics, model = model,
         formula = formula),
    class = "survalis_eval"
  )
}

#' @export
print.survalis_eval <- function(x, ...) {
  cat(sprintf("<survalis_eval> model = %s\n", x$model))
  print(x$resampling)
  cat(sprintf("  grid: %d times in [%s, %s]\n",
              length(x$times), format(min(x$times)), format(max(x$times))))
  s <- x$summary
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-7s %.4f  (sd %.4f, 95%% [%.4f, %.4f], n = %d)\n",
                s$metric[i], s$mean[i], s$sd[i], s$lower[i], s$upper[i], s$n[i]))
  }
  invisible(x)
}

#' @export
summary.survalis_eval <- function(object, ...) object$summary

#' @export
plot.survalis_eval <- function(x, title, ...) {
  if (missing(title)) title <- sprintf("Resampled performance: %s", x$model)
  cv_plot(x$per_split, title = title)
}
