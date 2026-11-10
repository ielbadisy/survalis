
#' Concordance Index from a Survival-Probability Matrix
#'
#' Computes Harrell's concordance index using predicted survival probabilities
#' at a specified time point (or the last column if \code{t_star} is \code{NULL}).
#'
#' @param object A \code{Surv} object of length \eqn{n}.
#' @param predicted An \code{n x k} matrix or data frame of survival probabilities
#'   with columns named \code{"t=<time>"}.
#' @param t_star Optional numeric time at which to evaluate the c-index; if omitted,
#'   the rightmost column of \code{predicted} is used.
#'
#' @details
#' Risk scores are defined as \code{1 - S(t)} at the chosen time. Ties receive
#' partial credit (0.5). Pairs not comparable due to censoring are excluded.
#'
#' @return A named numeric scalar: \code{"C index"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' sp <- matrix(
#'   stats::plogis(scale(veteran$karno)),
#'   ncol = 1,
#'   dimnames = list(NULL, "t=100")
#' )
#' cindex_survmat(y, predicted = sp, t_star = 100)
#' @noRd

cindex_survmat <- function(object, predicted, t_star = NULL) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")
  time <- object[, 1]
  status <- object[, 2]

  # Extract column for the given time point
  if (!is.null(t_star)) {
    t_name <- paste0("t=", t_star)
    if (!(t_name %in% colnames(predicted))) {
      stop("t_star = ", t_star, " not found in predicted survival matrix.")
    }
    surv_prob <- if (is.matrix(predicted)) predicted[, t_name] else predicted[[t_name]]
  } else {
    surv_prob <- if (is.matrix(predicted)) predicted[, ncol(predicted)] else predicted[[ncol(predicted)]]
  }

  # Harrell's C on the risk 1 - S(t). Pairs tied on risk count 1/2; pairs of
  # events tied in time are not comparable; an event tied in time with a
  # censored subject is comparable (the censored subject was still at risk).
  # This is the definition of survival::concordance(), which also runs in
  # O(n log n).
  risk_score <- 1 - surv_prob
  ok <- is.finite(time) & !is.na(status) & is.finite(risk_score)
  if (!any(status[ok] == 1)) return(stats::setNames(NA_real_, "C index"))
  cc <- survival::concordance(survival::Surv(time[ok], status[ok]) ~ risk_score[ok],
                              reverse = TRUE)
  C_index <- unname(cc$concordance)
  names(C_index) <- "C index"
  return(C_index)
}

#' Time-Dependent AUC from a Survival-Probability Matrix
#'
#' Computes a cumulative/dynamic time-dependent AUC using predicted survival
#' probabilities at a specified time point (or the last column if \code{t_star}
#' is \code{NULL}), following the Uno et al. (2007) / Heagerty and Zheng (2005)
#' estimator as implemented in the \pkg{timeROC} package. Cases are subjects
#' with an observed event strictly before \code{t_star}; controls are subjects
#' known to survive beyond \code{t_star}. Subjects censored before \code{t_star}
#' are handled through IPCW weighting.
#'
#' @param object A \code{Surv} object of length \eqn{n}.
#' @param predicted An \code{n x k} matrix or data frame of survival probabilities
#'   with columns named \code{"t=<time>"}.
#' @param t_star Optional numeric time at which to evaluate AUC; if omitted,
#'   the rightmost column of \code{predicted} is used.
#'
#' @details
#' Risk scores are defined as \code{1 - S(t)} at the chosen time. The AUC is
#' computed over case-control pairs using inverse-probability-of-censoring
#' weights for cases and partial credit (0.5) for ties.
#'
#' @return A named numeric scalar: \code{"auc"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' sp <- matrix(
#'   stats::plogis(scale(veteran$karno)),
#'   ncol = 1,
#'   dimnames = list(NULL, "t=100")
#' )
#' auc_survmat(y, predicted = sp, t_star = 100)
#' @noRd
auc_survmat <- function(object, predicted, t_star = NULL) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")

  time <- object[, 1]
  status <- object[, 2]

  if (!is.null(t_star)) {
    t_name <- paste0("t=", t_star)
    if (!(t_name %in% colnames(predicted))) {
      stop("t_star = ", t_star, " not found in predicted survival matrix.")
    }
    surv_prob <- if (is.matrix(predicted)) predicted[, t_name] else predicted[[t_name]]
  } else {
    if (is.matrix(predicted)) {
      surv_prob <- predicted[, ncol(predicted)]
      t_star <- as.numeric(sub("^t=", "", colnames(predicted)[ncol(predicted)]))
    } else {
      surv_prob <- predicted[[ncol(predicted)]]
      t_star <- as.numeric(sub("^t=", "", names(predicted)[ncol(predicted)]))
    }
  }

  risk_score <- 1 - surv_prob
  cases <- which(time < t_star & status == 1)
  controls <- which(time > t_star)

  if (!length(cases) || !length(controls)) {
    warning("AUC is undefined because no comparable case/control pairs are available.")
    return(stats::setNames(NA_real_, "auc"))
  }

  # IPCW with the marginal Kaplan-Meier of the censoring distribution:
  # cases are weighted by G(T_i-), controls by G(t), as in timeROC.
  km_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  Gt <- summary(km_fit, times = t_star, extend = TRUE)$surv
  idx <- findInterval(time[cases], km_fit$time, left.open = TRUE)
  G_case <- ifelse(idx > 0, km_fit$surv[pmax(idx, 1L)], 1)
  valid_cases <- which(is.finite(G_case) & G_case > 0)

  if (!length(valid_cases) || !is.finite(Gt) || Gt <= 0) {
    warning("AUC is undefined because censoring weights are not available at t_star.")
    return(stats::setNames(NA_real_, "auc"))
  }

  cases <- cases[valid_cases]
  case_weights <- 1 / G_case[valid_cases]
  risk_cases <- risk_score[cases]
  risk_controls <- risk_score[controls]

  comparisons <- outer(risk_cases, risk_controls, FUN = "-")
  concordant <- (comparisons > 0) + 0.5 * (comparisons == 0)

  auc_value <- sum(rowSums(concordant) * case_weights) /
    (length(controls) * sum(case_weights))

  names(auc_value) <- "auc"
  auc_value
}

#' Time-Dependent ROC/AUC Curve from a Survival-Probability Matrix
#'
#' Computes the cumulative/dynamic time-dependent AUC (Uno et al., 2007;
#' Heagerty and Zheng, 2005) at each time in \code{times}, using marginal
#' (Kaplan-Meier) inverse-probability-of-censoring weighting for cases. This
#' is a thin vectorized wrapper around \code{auc_survmat} that returns
#' a full \eqn{AUC(t)} curve instead of a single time point, matching the
#' quantity reported by \code{timeROC::timeROC()} with
#' \code{weighting = "marginal"}.
#'
#' @param object A \code{Surv} object of length \eqn{n}.
#' @param predicted An \code{n x k} matrix or data frame of survival probabilities
#'   with columns named \code{"t=<time>"}, one column per entry of \code{times}.
#' @param times Numeric vector of evaluation time points; each must have a
#'   matching \code{"t=<time>"} column in \code{predicted}.
#'
#' @return A \code{data.frame} (class \code{"timeroc_survmat"}) with columns
#'   \code{time} and \code{auc}, one row per entry of \code{times}.
#'
#' @seealso `auc_survmat()`
#'
#' @references
#' Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for
#' t-year survivors with censored regression models. \doi{10.1198/016214507000000149}
#'
#' Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC
#' curves. \doi{10.1111/j.0006-341X.2005.030814.x}
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' times <- c(30, 90, 180)
#' sp <- sapply(times, function(t) stats::plogis(-scale(veteran$karno) + log(t) / 5))
#' colnames(sp) <- paste0("t=", times)
#' timeroc_survmat(y, predicted = sp, times = times)
#' @noRd
timeroc_survmat <- function(object, predicted, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")
  if (!length(times)) stop("times must be a non-empty numeric vector.")

  auc_vals <- vapply(times, function(t_star) {
    unname(auc_survmat(object, predicted = predicted, t_star = t_star))
  }, numeric(1))

  out <- data.frame(time = times, auc = auc_vals)
  class(out) <- c("timeroc_survmat", class(out))
  out
}

#' Brier Score with IPCW for a Single Time Point
#'
#' Computes the inverse-probability-of-censoring weighted (IPCW) Brier score
#' at a single time \code{t_star}.
#'
#' @param object A \code{Surv} object.
#' @param pre_sp Numeric vector of predicted survival probabilities \eqn{S(t^{*}\mid x_i)}.
#' @param t_star Numeric evaluation time.
#'
#' @details
#' The censoring distribution \eqn{G(t)} is estimated via Kaplan-Meier on
#' \code{1 - status}. Observed events at or before \code{t_star} contribute
#' \eqn{S(t_i)^2 / G(t_i-)}; subjects event free after \code{t_star} contribute
#' \eqn{(1 - S(t^{*}))^2 / G(t^{*})} (Graf et al. 1999). This is the
#' \code{"brier"} member of [dwscore()]. Returns \code{NA} if \eqn{G(t^{*})} is
#' undefined or zero.
#'
#' @return A named numeric scalar: \code{"brier"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' pre_sp <- stats::plogis(scale(veteran$karno))
#' brier(y, pre_sp = pre_sp, t_star = 100)
#' @noRd

brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")
  bs_value <- suppressWarnings(.dw_ipcw(object, 1 - pre_sp, t_star, .dw_weight("brier")))
  if (is.na(bs_value)) {
    warning("Brier: Gt(t_star) is NA or zero at t = ", t_star, "; returning NA.")
    return(NA_real_)
  }
  names(bs_value) <- "brier"
  bs_value
}

#' Integrated Brier Score (Discrete Integration)
#'
#' Computes the Integrated Brier Score (IBS) over a vector of times, using
#' discrete left Riemann integration of IPCW Brier scores.
#'
#' @param object A \code{Surv} object.
#' @param sp_matrix Matrix/data frame of survival probabilities with one column
#'   per time in \code{times}.
#' @param times Numeric vector of strictly increasing times (length must equal
#'   \code{ncol(sp_matrix)}).
#'
#' @return A named numeric scalar: \code{"ibs"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' times <- c(60, 120)
#' lp <- stats::plogis(scale(veteran$karno))
#' sp <- cbind("t=60" = pmin(1, lp + 0.05), "t=120" = pmax(0, lp - 0.05))
#' colnames(sp) <- c("t=60", "t=120")
#' ibs_survmat(y, sp_matrix = sp, times = times)
#' @noRd

## Columns named "t=<time>" (as predict() returns them, in increasing time)
## are matched to `times` by name, so the order of `times` does not matter.
#' @keywords internal
#' @noRd
.align_survmat <- function(sp_matrix, times) {
  cn <- colnames(sp_matrix)
  if (is.null(cn) || !all(grepl("^t=", cn))) return(sp_matrix)
  ct <- suppressWarnings(as.numeric(sub("^t=", "", cn)))
  idx <- match(times, ct)
  if (anyNA(idx)) {
    idx <- vapply(times, function(t) {
      j <- which(abs(ct - t) <= 1e-8 * max(1, abs(t)))
      if (length(j)) j[1] else NA_integer_
    }, integer(1))
  }
  if (anyNA(idx)) return(sp_matrix)
  sp_matrix[, idx, drop = FALSE]
}

ibs_survmat <- function(object, sp_matrix, times) {
  sp_matrix <- .align_survmat(sp_matrix, times)
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")
  if (length(times) == 1) {
    return(brier(object, pre_sp = sp_matrix[, 1], t_star = times[1]))
  }

  o <- order(times)
  times <- times[o]
  sp_matrix <- sp_matrix[, o, drop = FALSE]
  brier_vec <- vapply(seq_along(times), function(i) {
    brier(object, pre_sp = sp_matrix[, i], t_star = times[i])
  }, numeric(1))
  if (anyNA(brier_vec)) return(c(ibs = NA_real_))

  dt <- diff(times)
  integral <- sum(brier_vec[-length(brier_vec)] * dt)
  ibs_value <- integral / (max(times) - min(times))
  names(ibs_value) <- "ibs"
  return(ibs_value)
}

#' Integrated Absolute Error Against Kaplan-Meier
#'
#' Computes \eqn{\int \mid \bar{S}(t) - \hat{S}_{KM}(t) \mid \, dt} where
#' \eqn{\bar{S}(t)} is the mean predicted survival across subjects and
#' \eqn{\hat{S}_{KM}(t)} is the Kaplan-Meier estimate. Integration is carried
#' out over KM event times using a left Riemann sum.
#'
#' @param object A \code{Surv} object.
#' @param sp_matrix Matrix/data frame of survival probabilities (rows = subjects,
#'   columns aligned with \code{times}).
#' @param times Numeric vector of times corresponding to the columns of \code{sp_matrix}.
#'
#' @return A named numeric scalar: \code{"iae"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' times <- c(60, 120)
#' lp <- stats::plogis(scale(veteran$karno))
#' sp <- cbind("t=60" = pmin(1, lp + 0.05), "t=120" = pmax(0, lp - 0.05))
#' colnames(sp) <- c("t=60", "t=120")
#' iae_survmat(y, sp_matrix = sp, times = times)
#' @noRd

iae_survmat <- function(object, sp_matrix, times) {
  sp_matrix <- .align_survmat(sp_matrix, times)
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  o <- order(times)
  times <- times[o]
  mean_pred_surv <- colMeans(sp_matrix)[o]
  km_fit <- survfit(object ~ 1)
  ## Integrate over [min(times), max(times)] only: both curves are step
  ## functions, so evaluate them on the union of the evaluation times and the
  ## Kaplan-Meier jump times inside that window.
  jumps <- km_fit$time[km_fit$n.event > 0]
  g <- sort(unique(c(times, jumps[jumps > min(times) & jumps < max(times)])))
  if (length(g) < 2L) return(c(iae = NA_real_))
  km_g <- summary(km_fit, times = g, extend = TRUE)$surv
  p_g <- approx(x = times, y = mean_pred_surv, xout = g,
                method = "constant", rule = 2)$y
  iae <- sum((abs(km_g - p_g))[-length(g)] * diff(g))
  names(iae) <- "iae"
  return(iae)
}

#' Integrated Squared Error Against Kaplan-Meier
#'
#' Computes \eqn{\int ( \bar{S}(t) - \hat{S}_{KM}(t) )^2 \, dt} where
#' \eqn{\bar{S}(t)} is the mean predicted survival and \eqn{\hat{S}_{KM}(t)}
#' is the Kaplan-Meier curve. Integration uses KM event times and a left
#' Riemann sum.
#'
#' @param object A \code{Surv} object.
#' @param sp_matrix Matrix/data frame of survival probabilities (rows = subjects,
#'   columns aligned with \code{times}).
#' @param times Numeric vector of times corresponding to the columns of \code{sp_matrix}.
#'
#' @return A named numeric scalar: \code{"ise"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' times <- c(60, 120)
#' lp <- stats::plogis(scale(veteran$karno))
#' sp <- cbind("t=60" = pmin(1, lp + 0.05), "t=120" = pmax(0, lp - 0.05))
#' colnames(sp) <- c("t=60", "t=120")
#' ise_survmat(y, sp_matrix = sp, times = times)
#' @noRd

ise_survmat <- function(object, sp_matrix, times) {
  sp_matrix <- .align_survmat(sp_matrix, times)
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  o <- order(times)
  times <- times[o]
  mean_pred_surv <- colMeans(sp_matrix)[o]
  km_fit <- survfit(object ~ 1)
  ## Integrate over [min(times), max(times)] only: both curves are step
  ## functions, so evaluate them on the union of the evaluation times and the
  ## Kaplan-Meier jump times inside that window.
  jumps <- km_fit$time[km_fit$n.event > 0]
  g <- sort(unique(c(times, jumps[jumps > min(times) & jumps < max(times)])))
  if (length(g) < 2L) return(c(ise = NA_real_))
  km_g <- summary(km_fit, times = g, extend = TRUE)$surv
  p_g <- approx(x = times, y = mean_pred_surv, xout = g,
                method = "constant", rule = 2)$y
  ise <- sum(((km_g - p_g)^2)[-length(g)] * diff(g))
  names(ise) <- "ise"
  return(ise)
}




#' Expected Calibration Error (ECE) at a Single Time Point
#'
#' Computes a binned Expected Calibration Error at time \code{t_star}, using
#' quantile-based bins of predicted survival probabilities.
#'
#' @param object A \code{Surv} object.
#' @param sp_matrix Matrix or data frame of survival probabilities with rows =
#'   subjects and (optionally) columns named \code{"t=<time>"}.
#' @param t_star Numeric evaluation time (must be a single value).
#' @param n_bins Integer number of quantile bins (default \code{10}).
#' @param p Power for the L^p error (default \code{1} gives absolute error).
#' @param weighted Logical; if \code{TRUE} (default) weights bin errors by
#'   \code{n_k / N} where \code{n_k} is the bin size.
#'
#' @details
#' Let \eqn{\hat{S}_k(t^*)} be the mean predicted survival in bin \eqn{k} and
#' \eqn{\tilde{S}_k(t^*)} the Kaplan-Meier survival at \eqn{t^*} in that bin.
#' The (weighted) ECE is
#' \deqn{
#'   \mathrm{ECE}(t^*) = \left(
#'     \sum_k w_k \, |\hat{S}_k(t^*) - \tilde{S}_k(t^*)|^p
#'   \right)^{1/p}
#' }
#' with \eqn{w_k = n_k / N}.
#'
#' @return A named numeric scalar: \code{"ece"}.
#'
#' @seealso The calibration *curve* is available as
#'   \code{interpret(fit, method = "calibration")} (via
#'   \code{compute_calibration()}); \code{ece_survmat()} is the scalar summary of
#'   the same idea used by \code{evaluate(metrics = "ece")}.
#'
#' @examples
#' y <- survival::Surv(
#'   time = c(1, 2, 3, 4, 6, 7, 8, 9),
#'   event = c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' sp <- cbind("t=5" = c(0.15, 0.20, 0.35, 0.40, 0.55, 0.65, 0.75, 0.80))
#' ece_survmat(y, sp_matrix = sp, t_star = 5, n_bins = 4)
#' @noRd
ece_survmat <- function(object,
  sp_matrix,
  t_star,
  n_bins   = 10L,
  p        = 1,
  weighted = TRUE) {
if (!inherits(object, "Surv")) {
stop("object must be a survival object (from Surv()).")
}
if (length(t_star) != 1L) {
stop("t_star must be a single numeric value for ECE.")
}

# coerce to matrix
if (is.data.frame(sp_matrix)) {
sp_mat <- as.matrix(sp_matrix)
} else if (is.matrix(sp_matrix)) {
sp_mat <- sp_matrix
} else {
stop("sp_matrix must be a matrix or data frame of survival probabilities.")
}

if (nrow(sp_mat) != nrow(object)) {
stop("Number of rows in sp_matrix must match length of Surv object.")
}

# Extract column at t_star (or the only column if unlabelled)
t_name <- paste0("t=", t_star)
if (!is.null(colnames(sp_mat)) && t_name %in% colnames(sp_mat)) {
pre_sp <- sp_mat[, t_name]
} else if (ncol(sp_mat) == 1L) {
pre_sp <- sp_mat[, 1L]
} else {
stop("Column for t_star = ", t_star,
" not found in sp_matrix (expected colname '", t_name, "').")
}

time   <- object[, 1]
status <- object[, 2]

# quantile-based bins on predicted survival
qs <- stats::quantile(
pre_sp,
probs = seq(0, 1, length.out = n_bins + 1L),
na.rm = TRUE
)
qs <- unique(qs)

if (length(qs) <= 1L) {
# all predictions identical: single bin
bins <- rep.int(1L, length(pre_sp))
} else {
bins <- cut(pre_sp, breaks = qs, include.lowest = TRUE, labels = FALSE)
}

uniq_bins <- sort(unique(bins[!is.na(bins)]))

n_vec    <- numeric(length(uniq_bins))
mean_vec <- numeric(length(uniq_bins))
obs_vec  <- numeric(length(uniq_bins))

for (i in seq_along(uniq_bins)) {
b   <- uniq_bins[i]
idx <- which(bins == b)
if (!length(idx)) next

n_vec[i]    <- length(idx)
mean_vec[i] <- mean(pre_sp[idx], na.rm = TRUE)

# KM in bin
fit_b   <- survival::survfit(survival::Surv(time[idx], status[idx]) ~ 1)
summ_b  <- summary(fit_b, times = t_star, extend = TRUE)
obs_vec[i] <- if (length(summ_b$surv) == 0L) NA_real_ else summ_b$surv[1L]
}

valid <- which(!is.na(mean_vec) & !is.na(obs_vec) & n_vec > 0)
if (!length(valid)) {
warning("ECE: no valid bins (all NA or empty); returning NA.")
return(NA_real_)
}

n_vec    <- n_vec[valid]
mean_vec <- mean_vec[valid]
obs_vec  <- obs_vec[valid]

diff_p <- abs(mean_vec - obs_vec)^p

if (weighted) {
w <- n_vec / sum(n_vec)
ece <- (sum(w * diff_p))^(1 / p)
} else {
ece <- (mean(diff_p))^(1 / p)
}

names(ece) <- "ece"
ece
}

.score_metrics <- function(surv_obj, pred, times, metrics, digits = NULL) {
  values <- vapply(metrics, function(m) {
    th <- .parse_range_metric(m)
    if (!is.null(th) || m == "survacc") {
      if (length(times) != 1) stop(m, " requires a single time point.")
      risk <- 1 - as.matrix(pred)[, 1]
      wt <- .dw_weight(if (is.null(th)) "survacc" else "range", thresholds = th)
      loss <- .dw_ipcw(surv_obj, risk, times, wt)
      return(if (is.null(th)) 1 - loss else loss)
    }
    unname(switch(m,
      "cindex" = cindex_survmat(surv_obj, predicted = pred, t_star = max(times)),
      "auc"    = auc_survmat(surv_obj, predicted = pred, t_star = max(times)),
      "brier"  = {
        if (length(times) != 1) stop("Brier requires a single time point.")
        brier(surv_obj, pre_sp = pred[, 1], t_star = times)
      },
      "ibs" = ibs_survmat(surv_obj, sp_matrix = pred, times = times),
      "iae" = iae_survmat(surv_obj, sp_matrix = pred, times = times),
      "ise" = ise_survmat(surv_obj, sp_matrix = pred, times = times),
      "ece" = {
        if (length(times) != 1) stop("ECE requires a single time point.")
        ece_survmat(surv_obj, sp_matrix = pred, t_star = times)
      },
      stop("Unknown metric: ", m)
    ))
  }, numeric(1))

  data.frame(metric = metrics, value = if (is.null(digits)) values else round(values, digits), stringsAsFactors = FALSE)
}



#' Cross-Validate a Survival Learner (fold-mapped with \code{fmapn})
#'
#' Runs k-fold cross-validation for any pair of \code{fit_fun}/\code{pred_fun}
#' that follow the package's learner contracts, and returns tidy per-fold metric values.
#' Fold iteration is handled by \code{functionals::fmapn()} with optional parallel
#' execution and a progress bar.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing all variables in \code{formula}.
#' @param fit_fun Function with signature \code{fit_fun(formula, data, ...)} that
#'   returns an \code{mlsurv_model}.
#' @param pred_fun Function with signature \code{pred_fun(object, newdata, times, ...)}
#'   returning a survival-probability matrix/data frame (columns named \code{"t=<time>"}).
#' @param times Numeric vector of evaluation times (passed to \code{pred_fun}).
#' @param metrics Character vector of metrics to compute. Supported:
#'   \code{"cindex"}, \code{"brier"} (single time), \code{"ibs"}, \code{"iae"}, \code{"ise"}, \code{"ece"} (single time).
#' @param folds Integer; number of folds (default \code{5}).
#' @param seed Integer random seed for reproducibility (default \code{123}).
#' @param verbose Logical; print row-dropping due to missingness (default \code{FALSE}).
#' @param ncores Integer; number of CPU cores for \code{fmapn} parallel mapping
#'   (default \code{1}). Set \code{> 1} to enable parallelism.
#' @param pb Logical; show a progress bar during fold mapping (default
#'   \code{interactive()}).
#' @param ... Additional arguments forwarded to \code{fit_fun}.
#'
#' @details
#' The routine:
#' \enumerate{
#' \item Validates \code{Surv(...)} on the LHS and warns against using \code{.} in formulas.
#' \item Drops rows with missing values in any variables referenced by \code{formula}.
#' \item Supports \code{Surv(time, status == k)} by recoding the status to 0/1.
#' \item Builds stratified v-folds on the status indicator (\pkg{rsample}).
#' \item For each fold: fits on the analysis set, predicts on the assessment set, and computes metrics.
#' }
#'
#' Fold iteration is performed via \code{functionals::fmapn()}, which preserves
#' per-fold identifiers (\code{id}, \code{fold}) and returns a list ready for
#' base R's `rbind()`.
#'
#' @return A \code{data.frame} with columns: \code{splits} (rsample split object),
#'   \code{id}, \code{fold}, \code{metric}, and \code{value}.
#'
#' @examples
#' \donttest{
#' cv_results <- cv_survlearner(
#'   formula = Surv(time, status) ~ age + karno,
#'   data = veteran,
#'   fit_fun = fit_coxph,
#'   pred_fun = predict_coxph,
#'   times = c(80, 160),
#'   metrics = c("cindex", "ibs"),
#'   folds = 2,
#'   seed = 1
#' )
#' }
#'
#' @seealso \code{fmapn}
#' @importFrom functionals fmapn
#' @keywords survival cross-validation fmapn parallel
#' @keywords internal
#' @noRd

cv_survlearner <- function(formula, data,
  fit_fun, pred_fun,
  times,
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  verbose = FALSE,
  ncores = 1,
  pb = interactive(),
  ...) {

if ("." %in% all.vars(update(formula, . ~ 0))) {
warning("Please avoid using '.' in the formula. Specify all predictors explicitly.")
}

parsed_formula <- .parse_surv_formula(formula, data)
time_col <- parsed_formula$time_col
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

set.seed(seed)
folds_data <- rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col))

# cross-validation loop via fmapn
results <- functionals::fmapn(
list(split = folds_data$splits, id = folds_data$id, fold = seq_along(folds_data$splits)),
function(split, id, fold) {
train <- rsample::analysis(split)
test  <- rsample::assessment(split)

model <- fit_fun(formula = formula, data = train, ...)
pred  <- pred_fun(model, newdata = test, times = times)
pred  <- .finalize_survmat(pred, times = times)

parsed_formula <- .parse_surv_formula(formula, test)
surv_obj <- survival::Surv(time = test[[parsed_formula$time_col]],
                           event = .recode_status_vec(parsed_formula, test))

scored <- .score_metrics(surv_obj, pred, times, metrics)
data.frame(
  splits = I(replicate(nrow(scored), split, simplify = FALSE)),
  id = id,
  fold = fold,
  metric = scored$metric,
  value = scored$value,
  stringsAsFactors = FALSE
)}, ncores = ncores, pb = pb)
  .abort_on_worker_error(results)
  .rbind_fill_dt(results)
}



#' Summarize Cross-Validation Results
#'
#' Produces mean, standard deviation, standard error, and 95\% Wald intervals
#' for each metric returned by \code{cv_survlearner}.
#'
#' @param cv_results A tibble/data frame as returned by \code{cv_survlearner()}.
#' @param digits Optional number of decimal places for \code{mean}, \code{sd},
#'   \code{se}, \code{lower}, \code{upper}. \code{NULL} (default) keeps full
#'   precision, so that tuners and leaderboards do not rank on rounded values.
#'
#' @return A data.frame with columns: \code{metric}, \code{mean}, \code{sd}, \code{n},
#'   \code{se}, \code{lower}, \code{upper}.
#'
#' @examples
#' cv_results <- data.frame(
#'   metric = c("cindex", "cindex", "ibs", "ibs"),
#'   value = c(0.62, 0.66, 0.19, 0.21)
#' )
#' cv_summary(cv_results)
#' @keywords internal
#' @noRd

cv_summary <- function(cv_results, digits = NULL) {
  df <- as.data.frame(cv_results)
  mean_tbl <- basetable::aggregate(df, by = "metric", value = "value", fun = "mean", na.rm = TRUE)
  names(mean_tbl)[names(mean_tbl) == "value"] <- "mean"
  sd_tbl <- basetable::aggregate(df, by = "metric", value = "value", fun = "sd", na.rm = TRUE)
  names(sd_tbl)[names(sd_tbl) == "value"] <- "sd"
  df$.ok <- as.integer(!is.na(df$value))
  n_tbl <- basetable::aggregate(df, by = "metric", value = ".ok", fun = "sum")
  names(n_tbl)[names(n_tbl) == ".ok"] <- "n"

  out <- basetable::merge(mean_tbl, sd_tbl, by = "metric", sort = TRUE)
  out <- basetable::merge(out, n_tbl, by = "metric", sort = TRUE)
  out$se <- out$sd / sqrt(out$n)
  out$lower <- out$mean - 1.96 * out$se
  out$upper <- out$mean + 1.96 * out$se
  if (!is.null(digits)) {
    for (col in c("mean", "sd", "se", "lower", "upper")) {
      out[[col]] <- round(out[[col]], digits)
    }
  }
  out
}


#' Boxplot of Cross-Validation Metric Distributions
#'
#' Visualizes per-fold metric values from \code{cv_survlearner} using a
#' boxplot with jittered points.
#'
#' @param cv_results A tibble/data frame as returned by \code{cv_survlearner()}.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely (e.g., for journals
#'   requiring caption-only figures).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' cv_results <- data.frame(
#'   metric = c("cindex", "cindex", "ibs", "ibs"),
#'   value = c(0.62, 0.66, 0.19, 0.21)
#' )
#' cv_plot(cv_results)
#' @keywords internal
#' @noRd

cv_plot <- function(cv_results, title) {
  if (missing(title)) title <- "Cross-Validation Performance"
  p <- ggplot(cv_results, aes(x = metric, y = value)) +
    geom_boxplot(fill = .survalis_palette[1], outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    labs(x = "Metric", y = "Value") +
    theme_survalis()
  if (!is.null(title)) p <- p + labs(title = title)
  p
}


#' Score a Fitted Survival Model on Its Training Data
#'
#' Computes one or more performance metrics for a fitted \code{mlsurv_model},
#' predicting on the training data with the model's corresponding \code{predict_*}
#' function.
#'
#' @param model An object of class \code{"mlsurv_model"}.
#' @param times Numeric vector of evaluation times. For \code{"brier"}, must be a single time.
#' @param metrics Character vector of metrics to compute. Supported:
#'   \code{"cindex"}, \code{"auc"}, \code{"ibs"}, \code{"brier"},
#'   \code{"iae"}, \code{"ise"}, \code{"ece"}.
#'
#' @details
#' The function constructs the appropriate \code{predict_*} function name from
#' \code{model$learner}, predicts survival probabilities on \code{model$data},
#' builds a \code{Surv} object from \code{model$formula}, and computes the metrics.
#' If \code{"brier"} is requested with multiple \code{times}, an error is thrown.
#'
#' @return A data.frame with columns \code{metric} and \code{value}.
#'
#' @examples
#' fitted_model <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
#' score_survmodel(
#'   fitted_model,
#'   times = c(80, 160),
#'   metrics = c("cindex", "ibs")
#' )
#' @keywords internal
#' @noRd

score_survmodel <- function(model, times, metrics = c("cindex", "ibs", "brier", "iae", "ise")) {
  stopifnot(inherits(model, "mlsurv_model"))
  stopifnot(!missing(times))

  # check brier with multiple times
  if ("brier" %in% metrics && length(times) != 1) {
    stop("The Brier score requires a single evaluation time.\n",
         "Please use `times = t` with a single value or remove 'brier' from the metric list.")
  }

  # check ece with multiple times
  if ("ece" %in% metrics && length(times) != 1) {
  stop("The ECE requires a single evaluation time.\n",
       "Please use `times = t` with a single value or remove 'ece' from the metric list.")
  }

  # extract model fields
  learner <- model$learner
  formula <- model$formula
  data    <- model$data

  if (is.null(learner)) stop("Model object must have a 'learner' field.")

  pred_fun <- .predict_fun(model)
  data <- .complete_cases_df(data, all.vars(formula))
  sp_matrix <- .finalize_survmat(pred_fun(model, newdata = data, times = times), times = times)
  parsed <- .parse_surv_formula(formula, data)
  surv_obj <- survival::Surv(data[[parsed$time_col]], .recode_status_vec(parsed, data))

  # score each metric
  .score_metrics(surv_obj, sp_matrix, times, metrics)
}
