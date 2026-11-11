#' Decision-weighted scores for survival predictions
#'
#' Scores predicted risks at one or more horizons with a decision-weighted
#' proper score: the loss of the rule "act if predicted risk exceeds `c`",
#' averaged over thresholds `c` with a weight `w(c)`. The weight states which
#' risk thresholds matter. Every such score equals a weighted area under the
#' net benefit curve, with weight `(1 - c) w(c)`, so a threshold range
#' restricts model comparison to the decisions clinicians actually make.
#'
#' @param object A numeric matrix of survival probabilities (rows = subjects,
#'   columns = times), for example from [predict()].
#' @param y A [survival::Surv()] object with one entry per row of `object`.
#' @param times Evaluation times. `NULL` reads the `"times"` attribute or
#'   column names of the form `"t=<time>"`.
#' @param weight Threshold weight:
#'   \describe{
#'     \item{`"range"`}{uniform on `thresholds` (default); only decisions with
#'       a threshold in that range count.}
#'     \item{`"survacc"`}{spherical weight \eqn{(c^2 + (1 - c)^2)^{-3/2}}; the
#'       loss is `1 - SurvAcc`.}
#'     \item{`"brier"`}{constant weight 2; the loss is the Brier score.}
#'     \item{`"log"`}{weight \eqn{1 / (c (1 - c))}; the loss is the
#'       logarithmic score.}
#'     \item{`"beta"`}{Beta density with parameters `shape`.}
#'   }
#' @param thresholds Length-2 numeric range `c(lower, upper)` in `[0, 1]`
#'   for `weight = "range"` (default `c(0.05, 0.20)`).
#' @param shape Length-2 numeric Beta parameters for `weight = "beta"`
#'   (default `c(2, 18)`, mean threshold 0.10).
#'
#' @details
#' With risk \eqn{r = 1 - S(t \mid X)} and event status \eqn{D = I(T \le t)},
#' the loss is \eqn{L_w(r, 0) = \int_0^r w(c)\, c\, dc} for a subject event
#' free at `t` and \eqn{L_w(r, 1) = \int_r^1 w(c) (1 - c)\, dc} for a subject
#' with an event by `t`. Right censoring is handled by inverse probability of
#' censoring weights from the Kaplan-Meier estimate of the censoring
#' distribution: subjects event free at `t` are weighted by \eqn{1 / G(t)},
#' subjects with an observed event at \eqn{U \le t} by \eqn{1 / G(U-)}, and
#' subjects censored before `t` get weight zero. This assumes censoring is
#' independent of the covariates. Range and Beta weights integrate to one, so
#' those losses lie on the scale of an average regret per decision.
#'
#' In [evaluate()], [compare()] and [tune()], use the metric ids `"survacc"`
#' (SurvAcc, higher is better) and `"range:<lower>-<upper>"`, for example
#' `"range:0.05-0.20"` (loss, lower is better). Both need a single time.
#'
#' @return A data frame with columns `time`, `weight` and `loss` (lower is
#'   better for every weight).
#'
#' @examples
#' m <- fit(Surv(time, status) ~ age + karno, veteran[1:100, ], "coxph")
#' S <- predict(m, veteran[-(1:100), ], times = c(90, 180))
#' y <- Surv(veteran$time[-(1:100)], veteran$status[-(1:100)])
#' dwscore(S, y)
#' dwscore(S, y, weight = "brier")
#' dwscore(S, y, weight = "range", thresholds = c(0.10, 0.30))
#'
#' # the same scores inside the verbs (single time)
#' S180 <- predict(m, veteran[-(1:100), ], times = 180)
#' evaluate(S180, y = y, metrics = c("brier", "survacc", "range:0.05-0.20"))
#' @export
dwscore <- function(object, y, times = NULL,
                    weight = c("range", "survacc", "brier", "log", "beta"),
                    thresholds = c(0.05, 0.20), shape = c(2, 18)) {
  weight <- match.arg(weight)
  if (!inherits(y, "Surv")) stop("`y` must be a Surv() object.", call. = FALSE)
  P <- .survmat_input(object, NROW(y), times)
  wt <- .dw_weight(weight, thresholds, shape)
  loss <- vapply(seq_along(P$times), function(k) {
    .dw_ipcw(y, 1 - P$S[, k], P$times[k], wt)
  }, numeric(1))
  data.frame(time = P$times, weight = wt$name, loss = loss,
             stringsAsFactors = FALSE)
}

## Check a survival matrix and read its evaluation times.
#' @keywords internal
#' @noRd
.survmat_input <- function(object, n, times = NULL) {
  S <- unclass(as.matrix(object))
  if (!is.numeric(S)) stop("`object` must be a numeric matrix.", call. = FALSE)
  if (nrow(S) != n) {
    stop("`object` has ", nrow(S), " rows but `y` has ", n, " entries.",
         call. = FALSE)
  }
  if (is.null(times)) times <- attr(object, "times")
  if (is.null(times) && !is.null(colnames(S))) {
    times <- suppressWarnings(as.numeric(sub("^t=", "", colnames(S))))
    if (anyNA(times)) times <- NULL
  }
  if (is.null(times)) {
    stop("Give `times`, or name the columns \"t=<time>\".", call. = FALSE)
  }
  times <- as.numeric(times)
  if (length(times) != ncol(S)) {
    stop("`times` has length ", length(times), " but `object` has ", ncol(S),
         " columns.", call. = FALSE)
  }
  attributes(S) <- list(dim = dim(S))
  list(S = S, times = times)
}

## Closed-form losses l0(r) = L_w(r, 0) and l1(r) = L_w(r, 1) for each weight.
#' @keywords internal
#' @noRd
.dw_weight <- function(weight, thresholds = c(0.05, 0.20), shape = c(2, 18),
                       eps = 1e-12) {
  switch(weight,
    brier = list(name = "brier", l0 = function(r) r^2, l1 = function(r) (1 - r)^2),
    survacc = {
      nrm <- function(r) sqrt(r^2 + (1 - r)^2)
      list(name = "survacc", l0 = function(r) 1 - (1 - r) / nrm(r),
           l1 = function(r) 1 - r / nrm(r))
    },
    log = list(name = "log", l0 = function(r) -log(pmax(1 - r, eps)),
               l1 = function(r) -log(pmax(r, eps))),
    beta = {
      if (length(shape) != 2L || any(!is.finite(shape)) || any(shape <= 0)) {
        stop("`shape` must be two positive numbers.", call. = FALSE)
      }
      a <- shape[1]; b <- shape[2]
      list(name = sprintf("beta(%g,%g)", a, b),
           l0 = function(r) a / (a + b) * stats::pbeta(r, a + 1, b),
           l1 = function(r) b / (a + b) * stats::pbeta(r, a, b + 1, lower.tail = FALSE))
    },
    range = {
      if (length(thresholds) != 2L || any(!is.finite(thresholds)) ||
          thresholds[1] < 0 || thresholds[2] > 1 || thresholds[1] >= thresholds[2]) {
        stop("`thresholds` must be c(lower, upper) with 0 <= lower < upper <= 1.",
             call. = FALSE)
      }
      cl <- thresholds[1]; cu <- thresholds[2]; h <- 1 / (cu - cl)
      list(name = sprintf("range(%g,%g)", cl, cu),
           l0 = function(r) { x <- pmin(pmax(r, cl), cu); h * (x^2 - cl^2) / 2 },
           l1 = function(r) { x <- pmin(pmax(r, cl), cu); h * ((cu - x) - (cu^2 - x^2) / 2) })
    },
    stop("Unknown weight: ", weight, call. = FALSE)
  )
}

## IPCW weights at horizon t from the Kaplan-Meier censoring estimate:
## w0 = I(U > t) / G(t) (event free at t), w1 = I(U <= t, Delta = 1) / G(U-)
## (event by t). NULL if G(t) is NA or zero. Shared by the scores and the
## decision curve, so that loss = int w(c) (1 - c) (pibar - NB(c)) dc holds
## exactly in every sample.
#' @keywords internal
#' @noRd
.ipcw_weights <- function(y, t) {
  U <- y[, 1]; Delta <- y[, 2]
  km <- survival::survfit(survival::Surv(U, 1 - Delta) ~ 1)
  Gt <- summary(km, times = t, extend = TRUE)$surv
  if (!length(Gt) || !is.finite(Gt) || Gt <= 0) return(NULL)
  ev <- U <= t & Delta == 1
  idx <- findInterval(U[ev], km$time, left.open = TRUE)
  G_ev <- rep(1, length(idx)); G_ev[idx > 0L] <- km$surv[idx[idx > 0L]]
  w1 <- numeric(length(U)); w1[ev] <- 1 / G_ev
  list(w0 = as.numeric(U > t) / Gt, w1 = w1)
}

## IPCW decision-weighted loss at horizon t.
#' @keywords internal
#' @noRd
.dw_ipcw <- function(y, risk, t, wt) {
  w <- .ipcw_weights(y, t)
  if (is.null(w)) {
    warning("dwscore: G(t) is NA or zero at t = ", t, "; returning NA.", call. = FALSE)
    return(NA_real_)
  }
  s <- numeric(length(risk))
  surv <- w$w0 > 0; ev <- w$w1 > 0
  s[surv] <- w$w0[surv] * wt$l0(risk[surv])
  s[ev] <- w$w1[ev] * wt$l1(risk[ev])
  mean(s)
}

## Parse a range metric id such as "range:0.05-0.20"; NULL if not one.
#' @keywords internal
#' @noRd
.parse_range_metric <- function(metric) {
  m <- regmatches(metric, regexec("^range:([0-9.]+)-([0-9.]+)$", metric))[[1]]
  if (!length(m)) return(NULL)
  th <- suppressWarnings(as.numeric(m[2:3]))
  if (anyNA(th) || th[1] < 0 || th[2] > 1 || th[1] >= th[2]) {
    stop("Metric ", sQuote(metric), ": use \"range:<lower>-<upper>\" with ",
         "0 <= lower < upper <= 1.", call. = FALSE)
  }
  th
}

## TRUE for metric ids the scorers accept (fixed ids and range ids).
#' @keywords internal
#' @noRd
.is_metric <- function(metric) {
  metric %in% listmetrics()$metric | grepl("^range:", metric)
}

## Validate metric ids against the time grid before any model is fitted, so a
## bad request fails at once rather than after a full resampling run.
#' @keywords internal
#' @noRd
.check_metrics <- function(metrics, times) {
  bad <- metrics[!.is_metric(metrics)]
  if (length(bad)) {
    stop("Unknown metric(s): ", paste(sQuote(bad), collapse = ", "),
         ". See listmetrics().", call. = FALSE)
  }
  for (m in metrics[grepl("^range:", metrics)]) .parse_range_metric(m)
  single <- metrics[metrics %in% c("brier", "ece", "survacc") |
                      grepl("^range:", metrics)]
  if (length(single) && length(times) != 1L) {
    stop("Metric(s) ", paste(sQuote(single), collapse = ", "),
         " need a single evaluation time, but `times` has ", length(times),
         " values. Use `times = t` or drop them; \"auc\" and \"cindex\" use ",
         "max(times).", call. = FALSE)
  }
  invisible(TRUE)
}
