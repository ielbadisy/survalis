
#' Concordance Index from a Survival-Probability Matrix
#'
#' Computes Harrell's concordance index using predicted survival probabilities
#' at a specified time point (or the last column if \code{t_star} is \code{NULL}).
#'
#' @param object A \code{\link[survival]{Surv}} object of length \eqn{n}.
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
#' @export

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

  risk_score <- 1 - surv_prob
  permissible <- 0
  concord <- 0
  par_concord <- 0
  n <- length(time)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((time[i] < time[j] & status[i] == 0) | (time[j] < time[i] & status[j] == 0)) next
      if (time[i] == time[j] & status[i] == 0 & status[j] == 0) next
      permissible <- permissible + 1
      if (time[i] != time[j]) {
        if ((time[i] < time[j] & risk_score[i] > risk_score[j]) |
            (time[j] < time[i] & risk_score[j] > risk_score[i])) {
          concord <- concord + 1
        } else if (risk_score[i] == risk_score[j]) {
          par_concord <- par_concord + 0.5
        }
      } else {
        if (status[i] + status[j] > 0) {
          if (risk_score[i] == risk_score[j]) concord <- concord + 1
          else par_concord <- par_concord + 0.5
        }
      }
    }
  }

  C_index <- (concord + par_concord) / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}

#' Time-Dependent AUC from a Survival-Probability Matrix
#'
#' Computes a cumulative/dynamic time-dependent AUC using predicted survival
#' probabilities at a specified time point (or the last column if \code{t_star}
#' is \code{NULL}). Cases are subjects with an observed event by \code{t_star};
#' controls are subjects known to survive beyond \code{t_star}. Subjects
#' censored before \code{t_star} are handled through IPCW weighting.
#'
#' @param object A \code{\link[survival]{Surv}} object of length \eqn{n}.
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
#' @export
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
  cases <- which(time <= t_star & status == 1)
  controls <- which(time > t_star)

  if (!length(cases) || !length(controls)) {
    warning("AUC is undefined because no comparable case/control pairs are available.")
    return(stats::setNames(NA_real_, "auc"))
  }

  km_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  eval_times <- sort(unique(c(t_star, time[cases])))
  G_all <- summary(km_fit, times = eval_times, extend = TRUE)$surv
  names(G_all) <- as.character(eval_times)

  Gt <- G_all[as.character(t_star)]
  G_case <- G_all[as.character(time[cases])]
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
  round(auc_value, 6)
}

#' Brier Score with IPCW for a Single Time Point
#'
#' Computes the inverse-probability-of-censoring weighted (IPCW) Brier score
#' at a single time \code{t_star}.
#'
#' @param object A \code{\link[survival]{Surv}} object.
#' @param pre_sp Numeric vector of predicted survival probabilities \eqn{S(t^{*}\mid x_i)}.
#' @param t_star Numeric evaluation time.
#'
#' @details
#' The censoring distribution \eqn{G(t)} is estimated via Kaplan-Meier on
#' \code{1 - status}. Observed events before \code{t_star} contribute
#' \eqn{S(t_i)^2 / G(t_i)}; those at risk at \code{t_star} contribute
#' \eqn{(1 - S(t^{*}))^2 / G(t^{*})}. Returns \code{NA} if \eqn{G(t^{*})} is
#' undefined or zero.
#'
#' @return A named numeric scalar: \code{"brier"}.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' pre_sp <- stats::plogis(scale(veteran$karno))
#' brier(y, pre_sp = pre_sp, t_star = 100)
#' @export

brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  time <- object[, 1]
  status <- object[, 2]
  km_fit <- survfit(Surv(time, 1 - status) ~ 1)

  all_times <- sort(unique(c(t_star, time)))
  G_all <- summary(km_fit, times = all_times, extend = TRUE)$surv
  names(G_all) <- as.character(all_times)

  Gt <- G_all[as.character(t_star)]

  if (is.na(Gt) || Gt == 0) {
    warning("Brier: Gt(t_star) is NA or zero at t = ", t_star, "; returning NA.")
    return(NA_real_)
  }

  early_event <- which(time < t_star & status == 1)
  later_risk  <- which(time >= t_star)
  score_vec <- numeric(length(time))

  if (length(early_event) > 0) {
    Gti <- G_all[as.character(time[early_event])]
    valid <- which(!is.na(Gti) & Gti > 0)
    score_vec[early_event[valid]] <- (pre_sp[early_event[valid]]^2) / Gti[valid]
  }

  if (Gt > 0 && length(later_risk) > 0) {
    score_vec[later_risk] <- ((1 - pre_sp[later_risk])^2) / Gt
  }

  bs_value <- mean(score_vec, na.rm = TRUE)
  names(bs_value) <- "brier"
  return(round(bs_value, 6))
}

#' Integrated Brier Score (Discrete Integration)
#'
#' Computes the Integrated Brier Score (IBS) over a vector of times, using
#' discrete left Riemann integration of IPCW Brier scores.
#'
#' @param object A \code{\link[survival]{Surv}} object.
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
#' @export

ibs_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")
  if (length(times) == 1) {
    return(brier(object, pre_sp = sp_matrix[, 1], t_star = times[1]))
  }

  times <- sort(times)
  brier_vec <- vapply(seq_along(times), function(i) {
    brier(object, pre_sp = sp_matrix[, i], t_star = times[i])
  }, numeric(1))

  dt <- diff(times)
  integral <- sum(brier_vec[-length(brier_vec)] * dt)
  ibs_value <- integral / (max(times) - min(times))
  names(ibs_value) <- "ibs"
  return(round(ibs_value, 6))
}

#' Integrated Absolute Error Against Kaplan-Meier
#'
#' Computes \eqn{\int \mid \bar{S}(t) - \hat{S}_{KM}(t) \mid \, dt} where
#' \eqn{\bar{S}(t)} is the mean predicted survival across subjects and
#' \eqn{\hat{S}_{KM}(t)} is the Kaplan-Meier estimate. Integration is carried
#' out over KM event times using a left Riemann sum.
#'
#' @param object A \code{\link[survival]{Surv}} object.
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
#' @export

iae_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  pred_at_km <- approx(x = times, y = mean_pred_surv, xout = km_data$time,
                       method = "constant", rule = 2)$y

  dt <- diff(km_data$time)
  t_IAE <- abs(km_data$surv - pred_at_km)
  iae <- sum(t_IAE[-length(t_IAE)] * dt)
  names(iae) <- "iae"
  return(round(iae, 4))
}

#' Integrated Squared Error Against Kaplan-Meier
#'
#' Computes \eqn{\int ( \bar{S}(t) - \hat{S}_{KM}(t) )^2 \, dt} where
#' \eqn{\bar{S}(t)} is the mean predicted survival and \eqn{\hat{S}_{KM}(t)}
#' is the Kaplan-Meier curve. Integration uses KM event times and a left
#' Riemann sum.
#'
#' @param object A \code{\link[survival]{Surv}} object.
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
#' @export

ise_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  pred_at_km <- approx(x = times, y = mean_pred_surv, xout = km_data$time,
                       method = "constant", rule = 2)$y

  dt <- diff(km_data$time)
  t_ISE <- (km_data$surv - pred_at_km)^2
  ise <- sum(t_ISE[-length(t_ISE)] * dt)
  names(ise) <- "ise"
  return(round(ise, 4))
}




#' Expected Calibration Error (ECE) at a Single Time Point
#'
#' Computes a binned Expected Calibration Error at time \code{t_star}, using
#' quantile-based bins of predicted survival probabilities.
#'
#' @param object A \code{\link[survival]{Surv}} object.
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
#' @examples
#' y <- survival::Surv(
#'   time = c(1, 2, 3, 4, 6, 7, 8, 9),
#'   event = c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' sp <- cbind("t=5" = c(0.15, 0.20, 0.35, 0.40, 0.55, 0.65, 0.75, 0.80))
#' ece_survmat(y, sp_matrix = sp, t_star = 5, n_bins = 4)
#' @export
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
round(ece, 6)
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
#' \code{dplyr::bind_rows()}.
#'
#' @return A tibble with columns: \code{splits} (rsample split object),
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
#' @seealso \code{\link[functionals]{fmapn}}
#' @importFrom functionals fmapn
#' @keywords survival cross-validation fmapn parallel
#' @export



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
data <- tidyr::drop_na(data, dplyr::all_of(all_vars))
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
time_col <- parsed_formula$time_col
status_col <- parsed_formula$status_col
event_value <- parsed_formula$event_value

if (parsed_formula$recode_status) {
status_vector <- as.integer(test[[status_col]] == event_value)
} else {
status_vector <- test[[status_col]]
}

surv_obj <- survival::Surv(time = test[[time_col]], event = status_vector)

tibble::tibble(metric = metrics) |>
      dplyr::mutate(value = lapply(metric, function(metric) {
        switch(metric,
          "cindex" = cindex_survmat(surv_obj, predicted = pred, t_star = max(times)),
          "auc" = auc_survmat(surv_obj, predicted = pred, t_star = max(times)),
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
      stop("Unknown metric: ", metric)
    )
  })) |> tidyr::unnest(cols = value) |>
  dplyr::mutate(splits = list(split), id = id, fold = fold) |>
  dplyr::relocate(splits, id, fold)}, ncores = ncores, pb = pb)
  dplyr::bind_rows(results)
}



#' Summarize Cross-Validation Results
#'
#' Produces mean, standard deviation, standard error, and 95\% Wald intervals
#' for each metric returned by \code{\link{cv_survlearner}}.
#'
#' @param cv_results A tibble/data frame as returned by \code{cv_survlearner()}.
#'
#' @return A tibble with columns: \code{metric}, \code{mean}, \code{sd}, \code{n},
#'   \code{se}, \code{lower}, \code{upper}.
#'
#' @examples
#' cv_results <- tibble::tibble(
#'   metric = c("cindex", "cindex", "ibs", "ibs"),
#'   value = c(0.62, 0.66, 0.19, 0.21)
#' )
#' cv_summary(cv_results)
#' @export

cv_summary <- function(cv_results) {
  cv_results |>
    group_by(metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      n    = dplyr::n(),
      se   = sd/sqrt(n),
      lower = (mean - 1.96 * se),
      upper = (mean + 1.96 * se),
      .groups = "drop"
    )
}


#' Boxplot of Cross-Validation Metric Distributions
#'
#' Visualizes per-fold metric values from \code{\link{cv_survlearner}} using a
#' boxplot with jittered points.
#'
#' @param cv_results A tibble/data frame as returned by \code{cv_survlearner()}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' cv_results <- tibble::tibble(
#'   metric = c("cindex", "cindex", "ibs", "ibs"),
#'   value = c(0.62, 0.66, 0.19, 0.21)
#' )
#' cv_plot(cv_results)
#' @export

cv_plot <- function(cv_results) {
  ggplot(cv_results, aes(x = metric, y = value)) +
    geom_boxplot(fill = "skyblue", outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    labs(title = "Cross-Validation Performance", x = "Metric", y = "Value") +
    theme_minimal()
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
#' @return A tibble with columns \code{metric} and \code{value}.
#'
#' @examples
#' fitted_model <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
#' score_survmodel(
#'   fitted_model,
#'   times = c(80, 160),
#'   metrics = c("cindex", "ibs")
#' )
#' @export

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

  # get prediction function
  pred_fun_name <- paste0("predict_", learner)
  if (!exists(pred_fun_name, mode = "function")) {
    stop("Prediction function not found: ", pred_fun_name)
  }
  pred_fun <- get(pred_fun_name, mode = "function")
  sp_matrix <- pred_fun(model, newdata = data, times = times)
  sp_matrix <- .finalize_survmat(sp_matrix, times = times)

  # extract Surv object
  parsed_formula <- .parse_surv_formula(formula, data)
  time_col <- parsed_formula$time_col
  status_col <- parsed_formula$status_col
  event_value <- parsed_formula$event_value

  if (parsed_formula$recode_status) {
    status_vector <- as.integer(data[[status_col]] == event_value)
  } else {
    status_vector <- data[[status_col]]
  }

  surv_obj <- survival::Surv(time = data[[time_col]], event = status_vector)

  # score each metric
  tibble::tibble(metric = metrics) |>
    dplyr::mutate(value = purrr::map(metric, function(metric) {
      switch(metric,
        "cindex" = cindex_survmat(surv_obj, predicted = sp_matrix, t_star = max(times)),
        "auc"    = auc_survmat(surv_obj, predicted = sp_matrix, t_star = max(times)),
        "brier"  = brier(surv_obj, pre_sp = sp_matrix[, 1], t_star = times),
        "ibs"    = ibs_survmat(surv_obj, sp_matrix, times),
        "iae"    = iae_survmat(surv_obj, sp_matrix, times),
        "ise"    = ise_survmat(surv_obj, sp_matrix, times),
        "ece"    = ece_survmat(surv_obj, sp_matrix = sp_matrix, t_star = times),
        stop("Unknown metric: ", metric)
      )      
    })) |>

    tidyr::unnest(cols = value)
  }
