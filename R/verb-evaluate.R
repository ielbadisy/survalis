#' Performance of a survival learner
#'
#' Scores survival predictions in one of three ways:
#' \describe{
#'   \item{Resampled}{`evaluate(fit)` or `evaluate(formula, data, model)` refits
#'     the learner on each training split of `resampling` and scores it on the
#'     held-out split (internal validation).}
#'   \item{New data}{`evaluate(fit, newdata = test)` scores a fitted learner on
#'     data it was not trained on (external or temporal validation).}
#'   \item{Supplied predictions}{`evaluate(S, y = Surv(time, status))` scores an
#'     `n x length(times)` matrix of survival probabilities produced elsewhere,
#'     for example by another package.}
#' }
#'
#' @param object A `survalis_fit` from [fit()], a survival `formula`, or a
#'   numeric matrix of survival probabilities (rows = subjects, columns = times).
#' @param data When `object` is a formula: the data frame.
#' @param model When `object` is a formula: a learner id
#'   (`listlearners()$learner`) or a custom learner (see [fit()]).
#' @param spec When `object` is a formula: engine arguments (named list).
#' @param newdata When `object` is a `survalis_fit`: optional data frame with the
#'   outcome columns. If given, the fit is scored on it once and `resampling` is
#'   ignored.
#' @param y When `object` is a matrix: a [survival::Surv()] object with one
#'   entry per row of `object`.
#' @param times Evaluation grid. `NULL` picks a grid from the data quantiles and
#'   records it; for a matrix, `NULL` reads the `"times"` attribute or column
#'   names of the form `"t=<time>"`.
#' @param resampling A resampling spec from [cv()], [holdout()], [groupcv()] or
#'   [bootstrap()] (default `cv()`).
#' @param metrics Character vector from `listmetrics()$metric`. `"brier"`,
#'   `"ece"`, `"survacc"` and range scores such as `"range:0.05-0.20"` require
#'   a single `times` value (see [dwscore()]).
#' @param ncores Cores for the split loop (default `1`).
#' @param ... Unused.
#'
#' @return A `survalis_eval` object: `per_split` (long data frame), `summary`
#'   (mean / sd / se / 95% interval per metric), `times`, `resampling` (`NULL`
#'   when scored once), `metrics`, `model`.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
#'          spec = list(num.trees = 100))
#' evaluate(m, times = c(90, 180, 365), resampling = cv(3))
#'
#' # train / test split
#' idx <- seq_len(100)
#' m_tr <- fit(Surv(time, status) ~ age + karno, veteran[idx, ], "coxph")
#' evaluate(m_tr, newdata = veteran[-idx, ], times = c(90, 180))
#'
#' # predictions from anywhere
#' S <- predict(m_tr, veteran[-idx, ], times = c(90, 180))
#' evaluate(S, y = Surv(veteran$time[-idx], veteran$status[-idx]))
#' @export
evaluate <- function(object, ...) UseMethod("evaluate")

#' @rdname evaluate
#' @export
evaluate.survalis_fit <- function(object, newdata = NULL, times = NULL,
                                  resampling = cv(),
                                  metrics = c("cindex", "ibs"), ncores = 1, ...) {
  model <- .null_default(object$learner_def, object$model_id)
  if (!is.null(newdata)) {
    parsed <- attr(object, "surv_parsed")
    if (is.null(parsed)) parsed <- .parse_surv_formula(object$formula, object$data)
    need <- c(parsed$time_col, parsed$status_col)
    if (!all(need %in% names(newdata))) {
      stop("`newdata` must contain the outcome columns ",
           paste(sQuote(need), collapse = ", "), ".", call. = FALSE)
    }
    newdata <- .complete_cases_df(newdata, all.vars(object$formula))
    if (!is.null(times)) .check_metrics(metrics, times)
    S <- predict(object, newdata = newdata, times = times)
    y <- survival::Surv(newdata[[parsed$time_col]], .recode_status_vec(parsed, newdata))
    out <- evaluate.matrix(S, y = y, times = attr(S, "times"), metrics = metrics)
    out$model <- object$model_id
    out$formula <- object$formula
    return(out)
  }
  .evaluate_core(object$formula, object$data, model, object$spec,
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

#' @rdname evaluate
#' @export
evaluate.matrix <- function(object, y, times = NULL,
                            metrics = c("cindex", "ibs"), ...) {
  if (!inherits(y, "Surv")) stop("`y` must be a Surv() object.", call. = FALSE)
  inp <- .survmat_input(object, NROW(y), times)
  P <- inp$S; times <- inp$times
  .check_metrics(metrics, times)
  colnames(P) <- paste0("t=", times)
  sc <- .score_metrics(y, as.data.frame(P), times, metrics)
  per <- data.frame(id = "supplied", metric = sc$metric, value = sc$value,
                    stringsAsFactors = FALSE)
  summ <- data.frame(metric = sc$metric, mean = sc$value, sd = NA_real_,
                     n = 1L, se = NA_real_, lower = NA_real_, upper = NA_real_)
  structure(
    list(per_split = per, summary = summ, times = times, resampling = NULL,
         metrics = metrics, model = "supplied predictions", formula = NULL,
         n = nrow(P)),
    class = "survalis_eval"
  )
}

#' @keywords internal
.recode_status_vec <- function(parsed, df) {
  s <- df[[parsed$status_col]]
  if (isTRUE(parsed$recode_status)) return(as.integer(s == parsed$event_value))
  # survival's own coding rules: 0/1, 1/2 (1 = censored), or logical
  as.integer(survival::Surv(df[[parsed$time_col]], s)[, 2])
}

#' Shared split-loop used by evaluate() and compare()
#' @keywords internal
.evaluate_splits <- function(splits, formula, data, model, spec, times, metrics,
                             parsed, ncores = 1) {
  lrn <- .resolve_learner(model)
  fit_fun  <- lrn$fit
  pred_fun <- lrn$predict

  # One seed per split, drawn from the caller's RNG stream, so stochastic
  # learners give the same results with any `ncores`.
  seeds <- sample.int(.Machine$integer.max, length(splits))
  one <- function(sp, seed) {
    set.seed(seed)
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
    functionals::fmapn(list(sp = splits, seed = seeds), function(sp, seed) one(sp, seed),
                       ncores = ncores, pb = FALSE)
  } else {
    Map(one, splits, seeds)
  }
  .abort_on_worker_error(parts)
  basetable::rbindfill(unname(parts), fill = TRUE)
}

#' @keywords internal
.evaluate_core <- function(formula, data, model, spec, times, resampling,
                           metrics, ncores = 1) {
  .resolve_learner(model)
  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula))
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)
  .check_metrics(metrics, times)

  splits <- .make_splits(resampling, data, parsed$status_col)
  per <- .evaluate_splits(splits, formula, data, model, as.list(spec), times,
                          metrics, parsed, ncores)
  summ <- cv_summary(per)

  structure(
    list(per_split = per, summary = as.data.frame(summ), times = times,
         resampling = resampling, metrics = metrics, model = .learner_label(model),
         formula = formula),
    class = "survalis_eval"
  )
}

#' @export
print.survalis_eval <- function(x, ...) {
  cat(sprintf("<survalis_eval> model = %s\n", x$model))
  if (is.null(x$resampling)) {
    cat(sprintf("  scored once on %d subjects (no resampling)\n", x$n))
  } else {
    print(x$resampling)
  }
  cat(sprintf("  grid: %d times in [%s, %s]\n",
              length(x$times), format(min(x$times)), format(max(x$times))))
  s <- x$summary
  lab <- formatC(s$metric, width = -max(7L, nchar(s$metric)))
  for (i in seq_len(nrow(s))) {
    if (is.null(x$resampling)) {
      cat(sprintf("  %s %.4f\n", lab[i], s$mean[i]))
    } else {
      cat(sprintf("  %s %.4f  (sd %.4f, 95%% [%.4f, %.4f], n = %d)\n",
                  lab[i], s$mean[i], s$sd[i], s$lower[i], s$upper[i], s$n[i]))
    }
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
