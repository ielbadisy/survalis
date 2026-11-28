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
#'   Rows with a missing value in a model variable are dropped with a warning.
#' @param model When `object` is a formula: a learner id
#'   (`listlearners()$learner`) or a custom learner (see [fit()]).
#' @param spec When `object` is a formula: engine arguments (named list).
#' @param newdata When `object` is a `survalis_fit`: optional data frame with the
#'   outcome columns. If given, the fit is scored on it once and `resampling` is
#'   ignored.
#' @param y When `object` is a matrix: a [survival::Surv()] object with one
#'   entry per row of `object`.
#' @param times Evaluation grid. `NULL` picks a grid from the data quantiles and
#'   records it (for resampling, from the full data before splitting, so every split
#'   shares one grid); for a matrix, `NULL` reads the `"times"` attribute or column
#'   names of the form `"t=<time>"`.
#' @param resampling A resampling spec from [cv()], [holdout()], [groupcv()] or
#'   [bootstrap()] (default `cv()`).
#' @param metrics Character vector from `listmetrics()$metric`. `"ece"`,
#'   `"survacc"` and range scores such as `"range:0.05-0.20"` require a single
#'   `times` value (see [dwscore()]). `"brier"` with several `times` gives one
#'   metric per time, `"brier@<time>"`; with one time it is `"brier"`. `"cindex"` and `"auc"`
#'   use the last time of `times` only, so their value depends on the grid. `"iae"` and `"ise"`
#'   compare the mean predicted curve with the Kaplan-Meier curve (calibration in the large),
#'   and `"ece"` bins subjects by deciles of predicted survival (see [listmetrics()]).
#' @param ncores Number of parallel workers for the split loop (default `1`). Workers are
#'   fresh R sessions (PSOCK), each with one OpenMP thread, so parallel runs are safe after an
#'   OpenMP engine (xgboost, rsf) ran in the session. A custom learner must therefore not rely
#'   on objects of the global environment other than its own arguments.
#' @param ... Unused.
#'
#' @return A `survalis_eval` object: `per_split` (long data frame), `summary`
#'   (mean / sd / se / 95% interval per metric), `times`, `resampling` (`NULL`
#'   when scored once), `metrics`, `model`. The `se` is the split standard
#'   deviation over the square root of the number of splits and the interval is
#'   mean +/- 1.96 se; this treats the splits as independent, which they are not
#'   (training sets overlap), so the interval is too narrow. With
#'   `bootstrap(estimator = ".632", ".632+" or "optimism")` the `mean` is that
#'   estimate, `sd`, `se` and the interval are `NA`, and a `bootstrap` element holds,
#'   per metric, the apparent, out-of-bag and no-information scores, the relative
#'   overfitting rate `R`, the .632+ weight `w` and the three corrected estimates
#'   (see [bootstrap()]); `per_split` then also has each fit's score on all rows
#'   (`on_full`) and on its own resample (`on_inbag`).
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
    newdata <- .complete_cases_df(newdata, all.vars(object$formula), warn = TRUE)
    if (!is.null(times)) .check_metrics(metrics, times, multi = TRUE)
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
  .check_metrics(metrics, times, multi = TRUE)
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
                             parsed, ncores = 1, full = FALSE) {
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
    out <- data.frame(id = sp$id, metric = sc$metric, value = sc$value,
                      stringsAsFactors = FALSE)
    if (full) {
      # scores of this fit on all rows and on its own (in-bag) rows, for the optimism bootstrap
      pf <- .finalize_survmat(pred_fun(m, newdata = data, times = times), times = times)
      out$on_full <- .score_rows(pf, data, parsed, times, metrics)
      out$on_inbag <- .score_rows(pf, data, parsed, times, metrics, sp$analysis)
    }
    out
  }

  parts <- .pmap(list(splits, seeds), one, ncores)
  .abort_on_worker_error(parts)
  basetable::rbindfill(unname(parts), fill = TRUE)
}

#' @keywords internal
.evaluate_core <- function(formula, data, model, spec, times, resampling,
                           metrics, ncores = 1) {
  .resolve_learner(model)
  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, all.vars(formula), warn = TRUE)
  status <- .recode_status_vec(parsed, data)
  times <- .resolve_times(times, data[[parsed$time_col]], status)
  .check_metrics(metrics, times, multi = TRUE)

  splits <- .make_splits(resampling, data, parsed$status_col)
  est <- if (resampling$type == "bootstrap") .null_default(resampling$params$estimator, "oob") else "oob"
  per <- .evaluate_splits(splits, formula, data, model, as.list(spec), times,
                          metrics, parsed, ncores, full = est != "oob")
  summ <- as.data.frame(.order_bytime(cv_summary(per[c("id", "metric", "value")])))
  boot <- NULL
  if (est != "oob") {
    boot <- .boot_estimates(per, formula, data, model, as.list(spec), times, metrics, parsed,
                            .null_default(resampling$params$nperm, 20L))
    i <- match(summ$metric, boot$metric)
    summ$mean <- boot[[est]][i]
    # the resamples give no simple standard error for a corrected estimate
    summ[c("sd", "se", "lower", "upper")] <- NA_real_
  }

  structure(
    list(per_split = per, summary = summ, times = times,
         resampling = resampling, metrics = metrics, model = .learner_label(model),
         formula = formula, bootstrap = boot),
    class = "survalis_eval"
  )
}

# Metric values of prediction matrix P scored on rows `i` of data (all rows by default).
.score_rows <- function(P, data, parsed, times, metrics, i = seq_len(nrow(data))) {
  d <- data[i, , drop = FALSE]
  y <- survival::Surv(d[[parsed$time_col]], .recode_status_vec(parsed, d))
  .score_metrics(y, P[i, , drop = FALSE], times, metrics)$value
}

# Apparent, out-of-bag, no-information, .632, .632+ and optimism estimates per metric
# (see ?bootstrap). `per` holds the out-of-bag scores and the scores on all and on in-bag rows.
.boot_estimates <- function(per, formula, data, model, spec, times, metrics, parsed, nperm) {
  lrn <- .resolve_learner(model)
  set.seed(sample.int(.Machine$integer.max, 1L))
  m <- do.call(lrn$fit, c(list(formula = formula, data = data), spec))
  S <- .finalize_survmat(lrn$predict(m, newdata = data, times = times), times = times)
  y <- survival::Surv(data[[parsed$time_col]], .recode_status_vec(parsed, data))
  app <- .score_metrics(y, S, times, metrics)
  ni <- .noinf_scores(y, as.matrix(S), times, metrics, app$metric, nperm)
  agg <- function(v) tapply(v, factor(per$metric, app$metric), mean, na.rm = TRUE)[app$metric]
  e1 <- unname(agg(per$value))
  sg <- ifelse(vapply(app$metric, .metric_maximize, logical(1)), -1, 1)
  # on the loss scale (sign sg): R = 0 without overfitting or information, 1 at no information
  a <- sg * app$value; b <- sg * e1; g <- sg * ni
  R <- ifelse(b <= a | g <= a, 0, pmin((b - a) / (g - a), 1))
  w <- 0.632 / (1 - 0.368 * R)
  data.frame(metric = app$metric, apparent = app$value, oob = e1, noinf = ni, R = R, w = w,
             `.632` = 0.368 * app$value + 0.632 * e1, `.632+` = (1 - w) * app$value + w * e1,
             optimism = app$value + unname(agg(per$on_full - per$on_inbag)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

# No-information score of each metric: predictions and outcomes paired at random over all n^2
# pairs. The IPCW Brier score is quadratic in a common prediction c, B(c) = k0 - 2 k1 c + k2 c^2
# with weights that follow the outcome, so its n^2-pair mean is k0 - 2 k1 mean(s) + k2 mean(s^2),
# exact; other metrics average `nperm` random row permutations.
.noinf_scores <- function(y, S, times, metrics, labels, nperm) {
  n <- nrow(S)
  bk <- function(k) {
    B <- vapply(c(0, 0.5, 1), function(c) unname(brier(y, rep(c, n), times[k])), 0)
    k2 <- 2 * (B[3] + B[1] - 2 * B[2]); k1 <- (B[1] - B[3] + k2) / 2
    B[1] - 2 * k1 * mean(S[, k]) + k2 * mean(S[, k]^2)
  }
  nb <- NULL
  perm <- NULL
  out <- lapply(metrics, function(mt) {
    if (mt %in% c("brier", "ibs")) {
      if (is.null(nb)) nb <<- vapply(seq_along(times), bk, 0)
      if (mt == "ibs" && length(times) > 1) {
        o <- order(times)
        return(sum(nb[o][-length(times)] * diff(times[o])) / diff(range(times)))
      }
      return(nb)
    }
    if (is.null(perm)) perm <<- replicate(nperm, sample.int(n))
    mean(apply(perm, 2, function(p) .score_one(mt, y, S[p, , drop = FALSE], times)))
  })
  unlist(out)
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
    if (is.null(x$resampling) || !is.null(x$bootstrap)) {
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

#' Plot an evaluation
#'
#' `type = "metrics"` draws one panel per metric with the per-split values as a
#' box plot and points; each panel has its own value axis and its strip says
#' whether higher or lower is better. `type = "time"` draws the time-specific
#' metrics (`"brier"` evaluated at several `times`, stored as `"brier@<time>"`)
#' against time: the mean over the splits with its 95% interval, or the single
#' value when the predictions were scored once.
#'
#' @param x A `survalis_eval` object.
#' @param type `"metrics"` (default) or `"time"`.
#' @param reference For `type = "time"`: an optional second `survalis_eval` on
#'   the same times, drawn dashed in grey, for example the Kaplan-Meier curve
#'   given to everyone (`evaluate()` on a matrix).
#' @param labels For `type = "time"`: legend labels of `x` and `reference`.
#'   `NULL` uses their model names.
#' @param xlab For `type = "time"`: the time axis label.
#' @param time.labels For `type = "time"`: `NULL` labels the evaluation times as
#'   they are; a function maps them to labels, e.g. `function(t) t / 365`.
#' @param title,subtitle,caption Plot title, subtitle and caption; `NULL`
#'   (default) draws none.
#' @param ... Unused.
#' @return A ggplot.
#' @examples
#' m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
#' ev <- evaluate(m, times = c(30, 90, 180), resampling = cv(3, seed = 1),
#'                metrics = c("cindex", "brier"))
#' plot(ev)
#' km <- summary(survival::survfit(Surv(time, status) ~ 1, veteran), times = c(30, 90, 180))$surv
#' ref <- evaluate(matrix(km, nrow(veteran), 3, byrow = TRUE),
#'                 y = Surv(veteran$time, veteran$status), times = c(30, 90, 180),
#'                 metrics = "brier")
#' plot(ev, type = "time", reference = ref, labels = c("Cox model", "Kaplan-Meier"),
#'      xlab = "Days")
#' @export
plot.survalis_eval <- function(x, type = c("metrics", "time"), reference = NULL, labels = NULL,
                               xlab = "Time", time.labels = NULL, title = NULL,
                               subtitle = NULL, caption = NULL, ...) {
  type <- match.arg(type)
  if (type == "metrics") {
    return(cv_plot(x$per_split, title = title) + ggplot2::labs(subtitle = subtitle, caption = caption))
  }
  d <- .eval_bytime(x)
  if (!nrow(d)) stop("`type = \"time\"` needs a metric scored at several times, such as ",
                     "`metrics = \"brier\"` with several `times`.", call. = FALSE)
  if (is.null(labels)) labels <- c(x$model, reference$model)
  if (!is.null(reference)) {
    if (labels[1] == labels[2]) labels <- c("Model", "Reference")
    r <- .eval_bytime(reference)
    d <- rbind(cbind(d, curve = labels[1]), cbind(r, curve = labels[2]))
  } else d$curve <- labels[1]
  labels <- labels[seq_len(1 + !is.null(reference))]
  d$curve <- factor(d$curve, labels)
  col <- stats::setNames(c("black", "grey55")[seq_along(labels)], labels)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$time, .data$mean, colour = .data$curve,
                                       linetype = .data$curve, shape = .data$curve))
  if (any(is.finite(d$lower))) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                                    width = 0, linetype = "solid", na.rm = TRUE)
  }
  p + ggplot2::geom_line() + ggplot2::geom_point(size = 2, fill = "white") +
    ggplot2::scale_colour_manual(values = col) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22")[seq_along(labels)], labels)) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(19, 21)[seq_along(labels)], labels)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_x_continuous(breaks = unique(d$time),
                                labels = if (is.null(time.labels)) ggplot2::waiver() else time.labels) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(x = xlab, y = if (length(unique(d$metric)) == 1) .metric_name(d$metric[1]) else NULL,
                  colour = NULL, linetype = NULL, shape = NULL,
                  title = title, subtitle = subtitle, caption = caption) +
    theme_survalis() +
    ggplot2::theme(strip.text = if (length(unique(d$metric)) == 1) ggplot2::element_blank() else
      ggplot2::element_text(face = "bold"),
      legend.position = if (is.null(reference)) "none" else "bottom")
}

# "brier@365" before "brier@1095": time-specific rows in time order, the rest as before.
.order_bytime <- function(s) {
  tt <- suppressWarnings(as.numeric(sub("^[^@]*@?", "", s$metric)))
  s <- s[order(sub("@.*", "", s$metric), tt), , drop = FALSE]
  rownames(s) <- NULL
  s
}

# Time-specific metrics ("brier@365") of an evaluation: metric, time, mean, lower, upper.
.eval_bytime <- function(x) {
  s <- x$summary[grepl("@", x$summary$metric), , drop = FALSE]
  data.frame(metric = sub("@.*", "", s$metric), time = as.numeric(sub(".*@", "", s$metric)),
             mean = s$mean, lower = s$lower, upper = s$upper)
}

.metric_name <- function(m) switch(m, brier = "Brier score", m)
